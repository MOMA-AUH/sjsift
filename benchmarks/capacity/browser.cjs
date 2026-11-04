/* Opt-in report capacity measurements, using the pinned test Playwright tools. */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const {pathToFileURL} = require('node:url');
const {performance} = require('node:perf_hooks');
const {chromium} = require('playwright');
const [file, output, countText, scenario, timeoutText, navigationText] = process.argv.slice(2);
const count = Number(countText);
const timeout = Number(timeoutText);
const result = {status: 'incomplete', node: process.version,
  playwright: require('playwright/package.json').version, interactions_ms: {}, errors: []};
const save = () => fs.writeFileSync(output, JSON.stringify(result, null, 2) + '\n');
save();

(async () => {
  const browser = await chromium.launch({headless: true});
  result.chromium = browser.version();
  save();
  try {
    const context = await browser.newContext({viewport: {width: 1400, height: 1000}, serviceWorkers: 'block'});
    const page = await context.newPage();
    page.setDefaultTimeout(timeout);
    page.setDefaultNavigationTimeout(Number(navigationText));
    const url = pathToFileURL(file).href;
    page.on('pageerror', error => result.errors.push(String(error)));
    page.on('crash', () => { result.errors.push('page crashed'); save(); });
    page.on('console', message => { if (message.type() === 'error') result.errors.push(message.text()); });
    page.on('request', request => {
      if (request.url().split('#')[0] !== url) result.errors.push(`Unexpected request: ${request.url()}`);
    });
    await context.route(/^(https?|wss?):/, route => route.abort());
    const started = performance.now();
    await page.goto(url);
    await page.locator('#overview').waitFor({state: 'visible'});
    result.load_ms = performance.now() - started;
    assert.ok(result.load_ms <= Number(navigationText), 'Local opening exceeded navigation limit');
    save();
    const cdp = await context.newCDPSession(page);
    await cdp.send('Performance.enable');
    async function memory(label) {
      const metrics = await cdp.send('Performance.getMetrics');
      result[label] = Object.fromEntries(metrics.metrics.filter(m =>
        ['JSHeapUsedSize', 'JSHeapTotalSize', 'Nodes', 'Documents'].includes(m.name)).map(m => [m.name, m.value]));
      result[label].dom_counters = await cdp.send('Memory.getDOMCounters');
      save();
    }
    async function snapshot() {
      return page.evaluate(() => {
        let memberships = 0;
        let hash = 2166136261;
        let rows = 0;
        let rare = 0;
        for (const row of document.querySelectorAll('.read-row')) {
          rows++;
          memberships += row.dataset.support.split(' ').length;
          const name = row.querySelector('summary').title;
          if (name.startsWith('rare;')) rare++;
          const text = `${row.closest('.variant-detail').id}|${row.dataset.recordId}|${row.dataset.support}|${name}`;
          // Order-independent fingerprint: pairing reparents rows but preserves them.
          let recordHash = 2166136261;
          for (let i = 0; i < text.length; i++) recordHash = Math.imul(recordHash ^ text.charCodeAt(i), 16777619);
          hash = (hash + recordHash) >>> 0;
        }
        return {rows, memberships, hash, rare,
          support: [...document.querySelectorAll('.variant-detail table')].map(node => node.textContent),
          overview: [...document.querySelectorAll('#overview-table tbody tr')].map(node => node.textContent)};
      });
    }
    const before = await snapshot();
    assert.equal(before.rows, scenario === 'single' ? count : 2 * count - 2);
    assert.equal(before.memberships, scenario === 'single' ? count : 3 * count - 2);
    assert.equal(before.rare, scenario === 'single' ? 0 : 2);
    result.inventory = before;
    result.dom_elements = await page.locator('*').count();
    await memory('memory_loaded');
    async function interaction(name, action, check) {
      result.current_check = name;
      save();
      const start = performance.now();
      await action();
      await check();
      // Include a completed paint, not just event dispatch.
      await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
      const elapsed = performance.now() - start;
      result.interactions_ms[name] = elapsed;
      save();
      assert.ok(elapsed <= timeout, `${name} exceeded ${timeout} ms`);
    }
    await interaction('select_abundant', () => page.getByRole('link', {name: 'abundant', exact: true}).click(),
      () => page.locator('#variant-0').waitFor({state: 'visible'}));
    const group = page.locator('#pileup-v0');
    await interaction('expand', () => group.locator('.read-expand').click(),
      () => page.waitForFunction(() => document.querySelector('#pileup-v0 .read-expand').getAttribute('aria-expanded') === 'true'));
    const first = await group.locator('.read-row:visible').first().getAttribute('data-record-id');
    await interaction('next_page', () => group.locator('.read-next').click(), async () => {
      assert.notEqual(await group.locator('.read-row:visible').first().getAttribute('data-record-id'), first);
    });
    await interaction('previous_page', () => group.locator('.read-previous').click(), async () => {
      assert.equal(await group.locator('.read-row:visible').first().getAttribute('data-record-id'), first);
    });
    await interaction('mapping_unique', () => group.getByRole('combobox', {name: 'Mapping class'}).selectOption('unique'), async () => {
      assert.ok(await group.locator('.read-row:visible').evaluateAll(rows => rows.length > 0 && rows.every(row => row.dataset.mapping === 'unique')));
    });
    await interaction('mapping_all', () => group.getByRole('combobox', {name: 'Mapping class'}).selectOption('all'), async () => {
      assert.equal(await group.getByRole('combobox', {name: 'Mapping class'}).inputValue(), 'all');
    });
    await interaction('support_filter', () => group.getByRole('combobox', {name: 'Junction support'}).selectOption(scenario === 'single' ? '0' : '1'), async () => {
      const bucket = scenario === 'single' ? '0' : '1';
      assert.ok(await group.locator('.read-row:visible').evaluateAll((rows, value) => rows.length > 0 &&
        rows.every(row => row.dataset.support.split(' ').includes(value)), bucket));
    });
    await interaction('support_all', () => group.getByRole('combobox', {name: 'Junction support'}).selectOption('all'), async () => {
      assert.equal(await group.getByRole('combobox', {name: 'Junction support'}).inputValue(), 'all');
    });
    await interaction('pairing', () => group.getByRole('combobox', {name: 'Read layout'}).selectOption('pairs'), async () => {
      assert.ok(await group.locator('.read-pair:visible').count() > 0);
    });
    await interaction('single_layout', () => group.getByRole('combobox', {name: 'Read layout'}).selectOption('single'), async () => {
      assert.equal(await group.locator('.read-pair:visible').count(), 0);
    });
    if (scenario === 'multiple') {
      await interaction('select_rare', () => page.locator('#variant-picker').selectOption('1'),
        () => page.locator('#variant-1').waitFor({state: 'visible'}));
      const rare = page.locator('#pileup-v1');
      await interaction('rare_defining_support', () => rare.getByRole('combobox', {name: 'Junction support'}).selectOption('0'), async () => {
        assert.equal(await rare.locator('.read-row:visible').count(), 2);
        assert.ok(await rare.locator('.read-row:visible summary').evaluateAll(nodes => nodes.every(node => node.title.startsWith('rare;'))));
      });
    }
    assert.deepEqual(await snapshot(), before, 'Viewing changed support or alignment evidence');
    await memory('memory_after_interactions');
    assert.deepEqual(result.errors, []);
    result.status = 'passed';
    delete result.current_check;
    save();
  } catch (error) {
    result.status = 'failed';
    result.error = error.stack || String(error);
    save();
    throw error;
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
