"""Run one bounded readiness stage; reuse the published capacity watchdog."""
import argparse
import json
from pathlib import Path
import runpy
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('result', type=Path)
    parser.add_argument('--seconds', type=float, default=120)
    parser.add_argument('--memory-mib', type=int, default=2048)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command
    if command[:1] == ['--']:
        command = command[1:]
    if not command or min(args.seconds, args.memory_mib) <= 0:
        parser.error('a command and positive limits are required')
    if args.result.exists():
        parser.error('result already exists')
    watchdog = runpy.run_path(str(Path(__file__).parents[1] / 'capacity' / 'run.py'))
    import os
    if not watchdog['tree_rss'](os.getpid()):
        parser.error('real Linux /proc memory counters are required')
    record = {'status': 'incomplete', 'seconds_limit': args.seconds,
              'memory_limit_mib': args.memory_mib, 'command': command}
    args.result.write_text(json.dumps(record, indent=2) + '\n')
    record.update(watchdog['bounded'](command, args.result.with_suffix('.log'),
                                     args.seconds, args.memory_mib))
    args.result.write_text(json.dumps(record, indent=2) + '\n')
    return 0 if record['status'] == 'passed' else 1


if __name__ == '__main__':
    sys.exit(main())
