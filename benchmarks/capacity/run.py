"""Opt-in, bounded Linux report-capacity benchmark. See README.md."""
import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def tree_rss(pid):
    """Sum live descendant RSS (shared pages can be counted more than once)."""
    total = 0
    pending = [pid]
    seen = set()
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        try:
            fields = Path(f'/proc/{current}/statm').read_text().split()
            total += int(fields[1]) * os.sysconf('SC_PAGE_SIZE')
            # A process may fork from any thread, including Chromium utility threads.
            for child_list in Path(f'/proc/{current}/task').glob('*/children'):
                try:
                    pending.extend(map(int, child_list.read_text().split()))
                except (FileNotFoundError, ProcessLookupError):
                    pass
        except (FileNotFoundError, ProcessLookupError):
            pass
    return total


def bounded(command, log, timeout, memory_mib):
    started = time.monotonic()
    peak = 0
    status = 'passed'
    with log.open('w') as output:
        process = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while process.poll() is None:
                peak = max(peak, tree_rss(process.pid))
                if peak > memory_mib * 1024 ** 2:
                    status = 'memory_limit'
                    break
                if time.monotonic() - started > timeout:
                    status = 'timeout'
                    break
                time.sleep(.01)
        finally:
            # Include Chromium descendants on failure, timeout, or interruption.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
    if status == 'passed' and process.returncode:
        status = 'failed'
    return {'status': status, 'wall_seconds': time.monotonic() - started,
            'peak_tree_rss_bytes': peak, 'exit_code': process.returncode, 'log': str(log)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--sizes', type=int, nargs='+', default=[1000, 10000, 50000, 100000])
    parser.add_argument('--stage-seconds', type=float, default=120)
    parser.add_argument('--browser-seconds', type=float, default=90)
    parser.add_argument('--navigation-seconds', type=float, default=30)
    parser.add_argument('--interaction-seconds', type=float, default=10)
    parser.add_argument('--python-memory-mib', type=int, default=2048)
    parser.add_argument('--browser-memory-mib', type=int, default=4096)
    parser.add_argument('--reference-bases', type=int, default=1000000)
    parser.add_argument('--skip-browser', action='store_true')
    args = parser.parse_args()
    if sys.platform != 'linux':
        parser.error('Linux /proc is required for the process-tree memory watchdog')
    if any(n < 20 or n % 4 for n in args.sizes):
        parser.error('sizes must be multiples of four, at least 20')
    if args.reference_bases < 1000 or min(args.stage_seconds, args.browser_seconds,
            args.navigation_seconds, args.interaction_seconds, args.python_memory_mib, args.browser_memory_mib) <= 0:
        parser.error('limits must be positive; reference must have at least 1000 bases')
    if tree_rss(os.getpid()) == 0:
        parser.error('real /proc RSS counters are required; this sandbox masks process memory')
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    settings = vars(args).copy()
    settings['output'] = str(args.output)
    import sjsift
    import pysam
    metadata = {'timestamp_utc': datetime.now(timezone.utc).isoformat(), 'python': sys.version,
                'python_executable': sys.executable, 'sjsift': sjsift.__version__, 'sjsift_path': sjsift.__file__,
                'pysam': pysam.__version__, 'htslib': pysam.version.__htslib_version__,
                'os': platform.platform(), 'cpu_count': os.cpu_count(),
                'cpu_model': next((line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                                   if line.startswith('model name')), 'unknown'),
                'meminfo': Path('/proc/meminfo').read_text(), 'affinity_cpus': len(os.sched_getaffinity(0)),
                'dependencies': {name: importlib.metadata.version(name) for name in ['setuptools', 'pytest']},
                'git_revision': subprocess.check_output(['git', '-C', str(HERE), 'rev-parse', 'HEAD'], text=True).strip(),
                'git_dirty': bool(subprocess.check_output(['git', '-C', str(HERE), 'status', '--porcelain'], text=True))}
    result = {'settings': settings, 'environment': metadata, 'cases': [], 'largest_complete_records': None}
    def save():
        (args.output / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    save()  # Persist limits before executing any workload.
    for count in args.sizes:
        for scenario in ['single', 'multiple']:
            directory = args.output / f'{scenario}-{count}'
            directory.mkdir()
            def worker(action, *arguments):
                name = '-'.join([action, *map(str, arguments[:1])])
                metric = directory / f'{name}.json'
                stage = bounded([sys.executable, str(HERE / 'worker.py'), action, str(metric), str(directory),
                                 *map(str, arguments)], directory / f'{name}.log',
                                args.stage_seconds, args.python_memory_mib)
                if metric.exists():
                    stage['metrics'] = json.loads(metric.read_text())
                return stage
            generation = worker('generate', count, scenario, args.reference_bases)
            for format_name in ['bam', 'cram']:
                case = {'records': count, 'scenario': scenario, 'format': format_name, 'generation': generation}
                result['cases'].append(case)
                print(f'{scenario} {count} {format_name}', flush=True)
                for name in ['report', 'validation', 'browser']:
                    case[name] = {'status': 'incomplete'}
                save()
                if generation['status'] != 'passed':
                    continue
                case['report'] = worker('report', format_name)
                save()
                if case['report']['status'] != 'passed':
                    continue
                validation = worker('validate', format_name, count, scenario)
                metrics = validation.pop('metrics', {})
                case['validation'] = {**validation, **metrics}
                save()
                if validation['status'] != 'passed':
                    continue
                if args.skip_browser:
                    case['browser'] = {'status': 'skipped'}
                else:
                    browser_file = directory / f'{format_name}.browser.json'
                    case['browser'] = bounded(['node', str(HERE / 'browser.cjs'),
                        str(directory / f'{format_name}.html'), str(browser_file), str(count), scenario,
                        str(args.interaction_seconds * 1000), str(args.navigation_seconds * 1000)], directory / f'{format_name}.browser.log',
                        args.browser_seconds, args.browser_memory_mib)
                    if browser_file.exists():
                        case['browser']['metrics'] = json.loads(browser_file.read_text())
                save()
    completed = [count for count in args.sizes if all(
        case[stage]['status'] == 'passed' for case in result['cases'] if case['records'] == count
        for stage in ['generation', 'report', 'validation', 'browser'])]
    result['largest_complete_records'] = max(completed, default=None)
    save()
    return 0 if args.skip_browser and all(case['validation']['status'] == 'passed' for case in result['cases']) else (
        0 if len(completed) == len(args.sizes) else 1)


if __name__ == '__main__':
    sys.exit(main())
