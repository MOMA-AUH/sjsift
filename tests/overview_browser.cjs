/* Overview, navigation, and progressive enhancement at the rendered-report seam. */
const assert = require("node:assert/strict");
const {pathToFileURL} = require("node:url");

async function snapshotEvidence(page) {
  return page.evaluate(() => ({
    support: [...document.querySelectorAll(".variant-detail table")].map(table => table.textContent),
    overview: [...document.querySelectorAll("#overview-table tbody tr")]
      .map(row => [row.dataset.index, row.textContent]).sort(),
    reads: [...document.querySelectorAll(".read-row")].map(row => ({
      variant: row.closest(".variant-detail").id, record: row.dataset.recordId,
      membership: row.dataset.support, mapping: row.dataset.mapping, pair: row.dataset.pair,
      label: row.querySelector("summary").getAttribute("aria-label"),
    })).sort((a, b) => JSON.stringify(a).localeCompare(JSON.stringify(b))),
  }));
}

async function checkOverview(page, staticPage, file) {
  const url = pathToFileURL(file).href;
  async function waitForView(id) {
    await page.waitForURL(url + "#" + id);
    // Hash navigation completes before the hashchange handler renders the view.
    await page.locator("#" + id).waitFor({state: "visible"});
  }
  await page.goto(url);
  const evidence = await snapshotEvidence(page);
  assert.equal(evidence.reads.length, 10);
  const names = () => page.locator("#overview-table tbody tr:visible a").allTextContents();
  assert.deepEqual(await names(), ["Larger", "Large", "Beta", "Alpha", "Zero"]);
  await page.getByRole("combobox", {name: "Order", exact: true}).selectOption("catalog");
  assert.deepEqual(await names(), ["Alpha", "Zero", "Beta", "Large", "Larger"]);
  await page.getByRole("checkbox", {name: "Show only variants with support"}).check();
  assert.deepEqual(await names(), ["Alpha", "Beta", "Large", "Larger"]);
  await page.getByRole("searchbox").fill("LARGE");
  assert.deepEqual(await names(), ["Large", "Larger"]);
  await page.getByRole("searchbox").fill("no-such-variant");
  assert.deepEqual(await names(), []);
  assert.equal(await page.locator("#no-matches").isVisible(), true);
  await page.getByRole("searchbox").fill("");
  await page.getByRole("combobox", {name: "Order", exact: true}).selectOption("support");
  assert.deepEqual(await names(), ["Larger", "Large", "Beta", "Alpha"]);
  await page.getByRole("link", {name: "Larger", exact: true}).click();
  await waitForView("variant-4");
  assert.equal(await page.locator("#variant-4").isVisible(), true);
  assert.equal(await page.locator("#previous").isDisabled(), true);
  await page.locator("#next").click();
  await waitForView("variant-3");
  await page.goBack();
  await waitForView("variant-4");
  assert.equal(await page.locator("#variant-4").isVisible(), true);
  await page.goForward();
  await waitForView("variant-3");
  assert.equal(await page.locator("#variant-3").isVisible(), true);
  await page.locator("#variant-picker").selectOption("0");
  await waitForView("variant-0");
  assert.equal(await page.locator("#next").isDisabled(), true);
  await page.locator("#previous").click();
  await waitForView("variant-2");
  assert.equal(await page.locator("#variant-2").isVisible(), true);
  await page.getByRole("link", {name: "All variants"}).click();
  await waitForView("overview");
  assert.deepEqual(await names(), ["Larger", "Large", "Beta", "Alpha"]);
  await page.getByRole("checkbox", {name: "Show only variants with support"}).uncheck();
  await page.getByRole("combobox", {name: "Order", exact: true}).selectOption("catalog");
  assert.deepEqual(await names(), ["Alpha", "Zero", "Beta", "Large", "Larger"]);
  assert.deepEqual(await snapshotEvidence(page), evidence);
  await page.screenshot({path: file + ".png", fullPage: true});

  await staticPage.goto(url);
  assert.equal(await staticPage.locator("#controls").isVisible(), false);
  assert.equal(await staticPage.locator("#overview-table tbody tr:visible").count(), 5);
  assert.equal(await staticPage.locator(".variant-detail:visible").count(), 5);
  assert.deepEqual(await staticPage.locator("#variant-0 .counts").allTextContents(), ["8 unique + 2 multi"]);
  assert.deepEqual(await staticPage.locator("#variant-1 .counts").allTextContents(), ["0 unique + 0 multi"]);
  assert.match(await staticPage.locator("#variant-0").innerText(), /chr1:101–200 \(\+\)/);
  assert.match(await staticPage.locator(".provenance").innerText(), /synthetic.SJ.out.tab/);
  assert.match(await staticPage.locator(".provenance").innerText(), /overview.bam/);
  await staticPage.getByRole("link", {name: "Alpha", exact: true}).click();
  await staticPage.waitForURL(url + "#variant-0");
  const inspection = staticPage.locator("#variant-0 .read-inspection").first();
  await inspection.locator("summary").focus();
  await staticPage.keyboard.press("Enter");
  assert.equal(await inspection.getAttribute("open"), "");
  assert.equal(await inspection.locator(".read-label").innerText(), "overview-pair");
  await inspection.locator("summary").click();
  assert.equal(await inspection.getAttribute("open"), null);
  assert.deepEqual(await snapshotEvidence(staticPage), evidence);
  await staticPage.screenshot({path: file + ".no-js.png", fullPage: true});
}

module.exports = {checkOverview, snapshotEvidence};
