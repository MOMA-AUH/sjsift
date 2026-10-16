/* Optional real-browser validation; see paired_view_browser.py. */
const assert = require("node:assert/strict");
const {pathToFileURL} = require("node:url");
const {chromium} = require("playwright");

(async () => {
  const browser = await chromium.launch({headless: true});
  try {
    const page = await browser.newPage({viewport: {width: 1400, height: 1000}});
    const errors = [];
    const requests = [];
    page.on("pageerror", error => errors.push(error.message));
    page.on("request", request => requests.push(request.url()));
    await page.goto(pathToFileURL(process.argv[2]).href + "#variant-0");
    const group = page.locator(".alignment-group");
    const support = await page.locator(".variant-detail table").innerText();
    await group.locator(".read-expand").click();
    assert.equal(await group.locator(".read-row:visible").count(), 51);
    const geometryBefore = await group.locator(".read-row").first().locator(".cigar-op").evaluateAll(nodes =>
      nodes.map(node => ({left: node.getBoundingClientRect().left, width: node.getBoundingClientRect().width})));
    await group.locator(".read-pairs").check();
    const geometryAfter = await group.locator(".read-row").first().locator(".cigar-op").evaluateAll(nodes =>
      nodes.map(node => ({left: node.getBoundingClientRect().left, width: node.getBoundingClientRect().width})));
    assert.deepEqual(geometryAfter, geometryBefore);
    assert.equal(await group.locator(".read-pair:visible").count(), 50);
    assert.equal(await group.locator(".read-row:visible").count(), 100);
    assert.match(await group.locator(".read-page").innerText(), /100 of 129.*records.*50 pairs/);
    const firstPage = await group.locator(".read-row:visible").evaluateAll(rows => rows.map(row => row.dataset.recordId));
    await group.locator(".read-next").click();
    assert.equal(await group.locator(".read-row:visible").count(), 29);
    assert.equal(await group.locator(".read-pair:visible").count(), 13);
    const secondPage = await group.locator(".read-row:visible").evaluateAll(rows => rows.map(row => row.dataset.recordId));
    assert.equal(new Set([...firstPage, ...secondPage]).size, 129);
    assert.equal(await group.locator(".read-next").isDisabled(), true);
    assert.match(await group.locator(".read-page").innerText(), /End of embedded reads/);

    const separated = group.locator(".read-pair:visible:not(.pair-overlap)");
    assert.equal(await separated.count(), 1);
    const positions = await separated.locator(".read-geometry").evaluateAll(nodes => nodes.map(node => node.getBoundingClientRect().top));
    assert.equal(positions[0], positions[1]);
    assert.equal(await separated.locator(".mate-connector").evaluate(node => getComputedStyle(node).borderTopStyle), "dotted");
    const overlapping = group.locator(".pair-overlap:visible").first();
    const overlapY = await overlapping.locator(".read-geometry").evaluateAll(nodes => nodes.map(node => node.getBoundingClientRect().top));
    assert.ok(overlapY[1] > overlapY[0]);
    const multi = group.locator(".read-hatched .op-M").first();
    assert.match(await multi.evaluate(node => getComputedStyle(node).backgroundImage), /repeating-linear-gradient/);
    const gray = await group.locator(".read-muted .cigar-op").evaluateAll(nodes => nodes.map(node => getComputedStyle(node).backgroundColor));
    assert.ok(gray.length > 0 && gray.every(color => color === "rgb(146, 152, 158)"));
    const secondary = group.locator('.read-hit[title*="Flags: secondary"]');
    assert.match(await secondary.getAttribute("title"), /multimapping.*NH: 2/);
    await secondary.focus();
    await page.keyboard.press("Enter");
    assert.equal(await group.locator(".read-inspection[open]").count(), 1);
    // A separated mate remains clickable despite sharing a display lane.
    await separated.locator(".read-row").first().locator(".op-M").first().click();
    assert.equal(await separated.locator(".read-inspection[open]").count(), 1);

    await group.locator(".support-filter").selectOption("1");
    assert.equal(await group.locator(".read-row:visible").count(), 1);
    assert.equal(await group.locator(".read-pair:visible").count(), 0);
    await group.locator(".support-filter").selectOption("all");
    assert.match(await group.locator(".read-page").innerText(), /Page 1 of 2/);
    await group.locator(".read-filter").selectOption("unique");
    assert.match(await group.locator(".read-page").innerText(), /of 126 matching records/);
    await group.locator(".read-expand").click();
    assert.equal(await group.locator(".read-pair:visible").count(), 10);
    assert.equal(await group.locator(".read-row:visible").count(), 20);
    await group.locator(".read-pairs").uncheck();
    assert.equal(await group.locator(".read-pair").count(), 0);
    assert.equal(await group.locator(".read-row:visible").count(), 11);
    assert.equal(await group.locator(".read-row").count(), 129);
    assert.equal(await page.locator(".variant-detail table").innerText(), support);
    await group.locator(".read-filter").selectOption("all");
    await group.locator(".read-pairs").check();
    await group.locator(".read-expand").click();
    await group.locator(".read-next").click();
    await page.screenshot({path: process.argv[2] + ".png", fullPage: true});
    assert.ok(requests.every(url => url.startsWith("file:")));
    assert.deepEqual(errors, []);
    console.log("Paired-view browser checks passed.");
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
