/* Inspect a real report through file: with all external requests blocked. */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const crypto = require('node:crypto');
const {pathToFileURL} = require('node:url');
const {performance} = require('node:perf_hooks');
const {chromium} = require('playwright');
const [file, output] = process.argv.slice(2);
const result = {status: 'incomplete', node: process.version,
  playwright: require('playwright/package.json').version, interactions_ms: {}, errors: []};
const save = () => fs.writeFileSync(output, JSON.stringify(result, null, 2) + '\n');
save();

(async () => {
  const browser = await chromium.launch({headless: true});
  result.chromium = browser.version();
  const url = pathToFileURL(file).href;
  try {
    async function newPage(javaScriptEnabled) {
      const context = await browser.newContext({javaScriptEnabled, serviceWorkers: 'block',
        viewport: {width: 1400, height: 1000}});
      context.on('request', request => {
        if (request.url().split('#')[0] !== url || !request.isNavigationRequest())
          result.errors.push('Unexpected resource request');
      });
      await context.route(/^(https?|wss?):/, route => route.abort());
      await context.setOffline(true);
      const page = await context.newPage();
      page.setDefaultTimeout(10000);
      page.setDefaultNavigationTimeout(30000);
      page.on('pageerror', error => result.errors.push(String(error)));
      page.on('websocket', () => result.errors.push('Unexpected WebSocket'));
      page.on('console', message => { if (message.type() === 'error') result.errors.push(message.text()); });
      return page;
    }
    const page = await newPage(true);
    const start = performance.now();
    await page.goto(url);
    await page.locator('#overview').waitFor({state: 'visible'});
    result.load_ms = performance.now() - start;
    assert.ok(result.load_ms <= 30000);
    const evidence = async () => crypto.createHash('sha256').update(await page.evaluate(() => JSON.stringify({
      rows: [...document.querySelectorAll('.read-row')].map(row => [row.closest('.variant-detail').id,
        row.dataset.recordId, row.dataset.mapping, row.dataset.support, row.querySelector('summary').title]).sort(),
      counts: [...document.querySelectorAll('.variant-detail table')].map(table => table.textContent),
    }))).digest('hex');
    const before = await evidence();
    async function action(name, perform, inspect) {
      result.current_check = name;
      save();
      const started = performance.now();
      await perform();
      await inspect();
      await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
      result.interactions_ms[name] = performance.now() - started;
      save();
      assert.ok(result.interactions_ms[name] <= 10000, name);
    }
    result.variants = await page.locator('.variant-detail').count();
    assert.equal(result.variants, 17);
    await action('open_detail', () => page.locator('#overview-table a[href="#variant-0"]').click(),
      () => page.locator('#variant-0').waitFor({state: 'visible'}));
    // Exercise every navigation destination, including variants with zero observations.
    for (let i = 0; i < 17; i++) {
      await action(`navigate_${i}`, () => page.locator('#variant-picker').selectOption(String(i)),
        () => page.locator(`#variant-${i}`).waitFor({state: 'visible'}));
    }
    const choices = await page.locator('.variant-detail').evaluateAll(nodes => nodes.map((node, index) => ({
      index, rows: node.querySelectorAll('.read-row').length,
      paired: node.querySelectorAll('.read-row[data-pair]:not([data-pair=""])').length,
    })).sort((a, b) => b.paired - a.paired || b.rows - a.rows));
    const choice = choices[0];
    assert.ok(choice.rows > 0, 'No selected public alignment evidence');
    await action('select_evidence', () => page.locator('#variant-picker').selectOption(String(choice.index)),
      () => page.locator(`#variant-${choice.index}`).waitFor({state: 'visible'}));
    const group = page.locator(`#pileup-v${choice.index}`);
    await action('expand', () => group.locator('.read-expand').click(), async () => {
      assert.equal(await group.locator('.read-expand').getAttribute('aria-expanded'), 'true');
    });
    const mapping = group.getByRole('combobox', {name: 'Mapping class'});
    for (const value of await mapping.locator('option').evaluateAll(nodes => nodes.map(node => node.value))) {
      if (value === 'all') continue;
      await action(`mapping_${value}`, () => mapping.selectOption(value), async () => {
        assert.ok(await group.locator('.read-row:visible').evaluateAll((rows, expected) =>
          rows.length > 0 && rows.every(row => row.dataset.mapping === expected), value));
      });
    }
    await action('mapping_reset', () => mapping.selectOption('all'), async () => {
      assert.ok(await group.locator('.read-row:visible').count() > 0);
    });
    const support = group.getByRole('combobox', {name: 'Junction support'});
    for (const value of await support.locator('option').evaluateAll(nodes => nodes.map(node => node.value))) {
      if (value === 'all') continue;
      await action(`support_${value}`, () => support.selectOption(value), async () => {
        assert.ok(await group.locator('.read-row:visible').evaluateAll((rows, expected) =>
          rows.every(row => row.dataset.support.split(' ').includes(expected)), value));
        const eligible = await group.locator(`.read-row[data-support~="${value}"]`).count();
        assert.equal((await group.locator('.read-row:visible').count()) > 0, eligible > 0);
      });
    }
    await support.selectOption('all');
    const first = await group.locator('.read-row:visible').first().getAttribute('data-record-id');
    if (await group.locator('.read-next').isEnabled()) {
      await action('next_page', () => group.locator('.read-next').click(), async () => {
        assert.notEqual(await group.locator('.read-row:visible').first().getAttribute('data-record-id'), first);
      });
      await action('previous_page', () => group.locator('.read-previous').click(), async () => {
        assert.equal(await group.locator('.read-row:visible').first().getAttribute('data-record-id'), first);
      });
    } else result.pagination = 'not applicable: one page';
    await action('pairs', () => group.getByRole('combobox', {name: 'Read layout'}).selectOption('pairs'), async () => {
      result.visible_pairs = await group.locator('.read-pair:visible').count();
      assert.ok(result.visible_pairs > 0, 'No safely associated mate pair in selected evidence');
    });
    await action('read_inspection', () => group.locator('.read-row:visible summary').first().click(), async () => {
      const metadata = await group.locator('.read-inspection[open] .read-metadata').first().innerText();
      assert.match(metadata, /CIGAR:/);
      assert.match(metadata, /NH:/);
      assert.match(metadata, /Flags:/);
    });
    await action('single', () => group.getByRole('combobox', {name: 'Read layout'}).selectOption('single'), async () => {
      assert.equal(await group.locator('.read-pair:visible').count(), 0);
    });
    result.multimapping_rows = await page.locator('.read-row[data-mapping="multimapping"]').count();
    if (result.multimapping_rows) {
      const multiIndex = await page.locator('.read-row[data-mapping="multimapping"]').first()
        .evaluate(row => row.closest('.variant-detail').id.split('-')[1]);
      await action('select_multimapping_view', () => page.locator('#variant-picker').selectOption(multiIndex),
        () => page.locator(`#variant-${multiIndex}`).waitFor({state: 'visible'}));
      const multiGroup = page.locator(`#pileup-v${multiIndex}`);
      if (await multiGroup.locator('.read-expand').getAttribute('aria-expanded') !== 'true') {
        await action('expand_multimapping_view', () => multiGroup.locator('.read-expand').click(), async () => {
          assert.equal(await multiGroup.locator('.read-expand').getAttribute('aria-expanded'), 'true');
        });
      }
      const multiMapping = multiGroup.getByRole('combobox', {name: 'Mapping class'});
      const multiRows = multiGroup.locator('.read-row[data-mapping="multimapping"]');
      await action('show_real_multimapping', () => multiMapping.selectOption('all'), async () => {
        assert.ok(await multiRows.filter({visible: true}).count() > 0);
      });
      await action('exclude_real_multimapping', () => multiMapping.selectOption('unique'), async () => {
        assert.equal(await multiRows.filter({visible: true}).count(), 0);
        assert.ok(await multiGroup.locator('.read-row:visible').evaluateAll(rows => rows.length > 0 &&
          rows.every(row => row.dataset.mapping === 'unique')));
      });
      await action('restore_real_multimapping', () => multiMapping.selectOption('all'), async () => {
        assert.ok(await multiRows.filter({visible: true}).count() > 0);
      });
      await action('inspect_real_multimapping', () => multiRows.filter({visible: true}).first().locator('summary').click(), async () => {
        const metadata = await multiRows.locator('.read-inspection[open] .read-metadata').first().innerText();
        const nh = metadata.match(/NH: (\d+)/);
        assert.ok(nh && Number(nh[1]) > 1);
      });
      result.multimapping_filter = 'passed';
    } else result.multimapping_filter = 'not applicable: no selected multimapping records';
    assert.equal(await evidence(), before, 'Viewing changed counts or retained evidence');
    result.evidence_sha256 = before;
    result.embedded_rows = await page.locator('.read-row').count();
    const staticPage = await newPage(false);
    await staticPage.goto(url + `#variant-${choice.index}`);
    assert.equal(await staticPage.locator('.variant-detail:visible').count(), 17);
    await staticPage.locator(`#pileup-v${choice.index} .read-row:visible summary`).first().click();
    assert.ok(await staticPage.locator('.read-inspection[open] .read-metadata').first().isVisible());
    result.javascript_disabled = 'passed';
    assert.deepEqual(result.errors, []);
    result.status = 'passed';
    delete result.current_check;
    save();
  } catch (error) {
    result.status = 'failed';
    result.error = error.stack || String(error);
    save();
    throw error;
  } finally { await browser.close(); }
})().catch(error => { console.error(error); process.exitCode = 1; });
