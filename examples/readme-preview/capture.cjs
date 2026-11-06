/* Capture the real report without changing its content or styles. */
const assert = require("node:assert/strict");
const path = require("node:path");
const {pathToFileURL} = require("node:url");
const {createRequire} = require("node:module");
const requireBrowser = createRequire(path.resolve(__dirname, "../../tests/package.json"));
const {chromium} = requireBrowser("playwright");

(async () => {
  if (process.argv.length !== 4) {
    throw new Error("Usage: node capture.cjs REPORT.html SCREENSHOT.png");
  }
  const browser = await chromium.launch({headless: true});
  try {
    const page = await browser.newPage({
      viewport: {width: 1200, height: 1400}, deviceScaleFactor: 2,
      colorScheme: "light", reducedMotion: "reduce", locale: "en-US",
    });
    const errors = [];
    page.on("pageerror", error => errors.push(error.message));
    page.on("console", message => {
      if (message.type() === "error") errors.push(message.text());
    });
    await page.route(/^(https?|wss?):/, route => {
      errors.push(`Unexpected external request: ${route.request().url()}`);
      return route.abort();
    });
    await page.goto(pathToFileURL(path.resolve(process.argv[2])).href);
    await page.getByRole("link", {name: "METex14", exact: true}).click();
    const detail = page.locator(".variant-detail:visible");
    await detail.locator(".junction-arc").first().waitFor();
    await page.evaluate(() => document.fonts.ready);
    assert.equal(await detail.locator(".junction-arc").count(), 3);
    assert.equal(await detail.locator(".read-row").count(), 306);
    assert.equal(await detail.locator(".read-row:visible").count(), 26);
    assert.equal(await detail.locator(".read-hatched").count(), 1);
    assert.equal(await detail.locator(".read-hatched:visible").count(), 1);
    const widths = await detail.locator(".junction-arc").evaluateAll(arcs =>
      arcs.map(arc => Number(arc.getAttribute("stroke-width"))));
    assert.ok(widths[1] > widths[0] + 1.5);
    assert.ok(widths[2] > widths[1] + 1);
    await page.mouse.move(0, 0);
    await detail.screenshot({path: process.argv[3], animations: "disabled"});
    assert.deepEqual(errors, []);
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
