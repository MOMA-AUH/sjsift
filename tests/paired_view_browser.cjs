/* Optional real-browser validation; see paired_view_browser.py. */
const assert = require("node:assert/strict");
const {pathToFileURL} = require("node:url");
const {chromium} = require("playwright");

(async () => {
  const browser = await chromium.launch({headless: true});
  try {
    if (process.argv[4]) {
      // Static report text and computed rendering must work without JavaScript.
      const legendPage = await browser.newPage({javaScriptEnabled: false, viewport: {width: 1400, height: 1000}});
      await legendPage.goto(pathToFileURL(process.argv[4]).href + "#variant-0");
      const styles = await legendPage.locator(".read-legend i").evaluateAll(nodes => nodes.map(node => ({
        color: getComputedStyle(node).backgroundColor, pattern: getComputedStyle(node).backgroundImage,
      })));
      assert.equal(styles.length, 3);
      const channels = styles[0].color.match(/\d+/g).map(Number);
      assert.equal(channels[0], channels[1], "Mapping legend must be achromatic");
      assert.equal(channels[1], channels[2], "Mapping legend must be achromatic");
      assert.equal(styles[1].color, styles[0].color);
      assert.equal(styles[0].pattern, "none");
      assert.match(styles[1].pattern, /repeating-linear-gradient/);
      assert.equal(styles[2].color, "rgb(146, 152, 158)");
      assert.equal(styles[2].pattern, "none");
      assert.notEqual(styles[0].color, styles[2].color);
      const arcs = await legendPage.locator(".junction-arc").evaluateAll(nodes => nodes.map(node => getComputedStyle(node).stroke));
      const bars = await legendPage.locator(".variant-detail .unique-support").evaluateAll(nodes => nodes.map(node => getComputedStyle(node).backgroundColor));
      assert.equal(new Set(arcs).size, 5);
      assert.deepEqual(new Set(arcs), new Set(bars));
      assert.deepEqual(bars.slice(0, 3), ["rgb(49, 89, 173)", "rgb(22, 131, 117)", "rgb(153, 80, 165)"]);
      for (const color of bars) assert.ok(!styles.some(style => style.color === color));
      for (let bucket = 0; bucket < 5; bucket++) {
        const rows = legendPage.locator(`.read-row[data-bucket="${bucket}"]`);
        assert.equal(await rows.count(), 7);
        for (const [kind, pattern] of [["solid", "none"], ["hatched", styles[1].pattern], ["muted", "none"]]) {
          const blocks = await legendPage.locator(`.read-row[data-bucket="${bucket}"].read-${kind} .op-M`).evaluateAll(nodes => nodes.map(node => ({
            color: getComputedStyle(node).backgroundColor, pattern: getComputedStyle(node).backgroundImage,
          })));
          assert.ok(blocks.length > 0);
          assert.ok(blocks.every(block => block.color === (kind === "muted" ? styles[2].color : bars[bucket]) && block.pattern === pattern));
        }
      }
      await legendPage.locator(".read-legend").screenshot({path: process.argv[4] + ".png"});
      await legendPage.close();
    }
    const page = await browser.newPage({viewport: {width: 1400, height: 1000}});
    const errors = [];
    const requests = [];
    page.on("pageerror", error => errors.push(error.message));
    page.on("request", request => requests.push(request.url()));
    await page.goto(pathToFileURL(process.argv[2]).href + "#variant-0");
    const group = page.locator(".alignment-group");
    assert.equal(await group.locator(".read-controls select").count(), 3);
    assert.deepEqual(await group.getByRole("combobox", {name: "Mapping class"}).locator("option").allTextContents(), ["All", "Unique"]);
    assert.deepEqual(await group.getByRole("combobox", {name: "Read layout"}).locator("option").allTextContents(), ["Single", "Pairs"]);
    assert.deepEqual(await group.getByRole("combobox", {name: "Junction support"}).locator("option").allTextContents(), ["All", "Variant", "next"]);
    const support = await page.locator(".variant-detail table").innerText();
    await group.locator(".read-expand").click();
    assert.equal(await group.locator(".read-row:visible").count(), 51);
    const geometryBefore = await group.locator(".read-row").first().locator(".cigar-op").evaluateAll(nodes =>
      nodes.map(node => ({left: node.getBoundingClientRect().left, width: node.getBoundingClientRect().width})));
    await group.getByRole("combobox", {name: "Read layout"}).selectOption("pairs");
    const geometryAfter = await group.locator(".read-row").first().locator(".cigar-op").evaluateAll(nodes =>
      nodes.map(node => ({left: node.getBoundingClientRect().left, width: node.getBoundingClientRect().width})));
    assert.deepEqual(geometryAfter, geometryBefore);
    assert.equal(await group.locator(".read-pair:visible").count(), 50);
    assert.equal(await group.locator(".pair-label").count(), 0);
    assert.equal(await group.locator(".read-pair-help").isVisible(), true);
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

    const separated = group.locator(".read-pair:visible").filter({has: page.locator('.read-hit[title*="separated"]')});
    assert.equal(await separated.count(), 1);
    const positions = await separated.locator(".read-geometry").evaluateAll(nodes => nodes.map(node => node.getBoundingClientRect().top));
    assert.equal(positions[0], positions[1]);
    assert.equal(await separated.locator(".mate-connector").evaluate(node => getComputedStyle(node).borderTopStyle), "dotted");
    const overlapping = group.locator(".pair-overlap:visible").first();
    const overlapY = await overlapping.locator(".read-geometry").evaluateAll(nodes => nodes.map(node => node.getBoundingClientRect().top));
    assert.ok(overlapY[1] > overlapY[0]);
    assert.ok(await group.locator(".read-row").evaluateAll(rows => rows.every(row =>
      row.querySelectorAll(".read-end-forward, .read-end-reverse").length === 1)));
    // At a narrow viewport, a one-base block is smaller than its arrowhead.
    await page.setViewportSize({width: 800, height: 1000});
    await page.waitForFunction(() => {
      const pair = document.querySelector('.read-hit[title*="separated"]').closest('.read-pair');
      const rows = pair.querySelectorAll('.read-geometry');
      return rows[1].getBoundingClientRect().top > rows[0].getBoundingClientRect().top;
    }, null, {timeout: 2000});
    for (const direction of ["forward", "reverse"]) {
      const terminal = group.locator(`.read-hit[title*="pair-61"] .read-end-${direction}`);
      const arrow = await terminal.evaluate(node => {
        const head = getComputedStyle(node, "::after");
        return {blockWidth: node.getBoundingClientRect().width, width: parseFloat(head.width),
          height: parseFloat(head.height), left: parseFloat(head.left), shape: head.clipPath};
      });
      assert.ok(arrow.blockWidth < 5);
      assert.ok(arrow.width >= 6 && arrow.height >= 10, "Short blocks need a full-size arrowhead");
      assert.match(arrow.shape, /^polygon\(/);
      if (direction === "forward") assert.ok(arrow.left >= arrow.blockWidth - 0.1);
      else assert.ok(arrow.left <= -arrow.width);
      // The protruding arrow itself still opens the read's inspection panel.
      await terminal.scrollIntoViewIfNeeded();
      const block = await terminal.boundingBox();
      await page.mouse.click(direction === "forward" ? block.x + block.width + 2 : block.x - 2,
        block.y + block.height / 2);
      assert.equal(await group.locator('.read-inspection[open] .read-label').innerText(), "pair-61");
      await group.locator('.read-inspection[open] > summary').click();
    }
    await group.locator('.read-pair').filter({has: page.locator('.read-hit[title*="pair-61"]')})
      .screenshot({path: process.argv[2] + ".short-arrows.png"});
    await page.setViewportSize({width: 1400, height: 1000});
    await page.waitForFunction(() => {
      const pair = document.querySelector('.read-hit[title*="separated"]').closest('.read-pair');
      const rows = pair.querySelectorAll('.read-geometry');
      return rows[1].getBoundingClientRect().top === rows[0].getBoundingClientRect().top;
    }, null, {timeout: 2000});
    const band = await overlapping.evaluate(node => getComputedStyle(node).backgroundColor);
    assert.notEqual(band, "rgba(0, 0, 0, 0)");
    assert.ok(await overlapping.evaluate(node => parseFloat(getComputedStyle(node, "::before").borderTopWidth) > 0));
    await overlapping.locator(".op-M").first().hover();
    assert.notEqual(await overlapping.evaluate(node => getComputedStyle(node).backgroundColor), band);
    await group.locator(".read-pairs").focus();
    await group.locator(".read-pairs").hover();
    await overlapping.locator(".read-hit").first().focus();
    assert.notEqual(await overlapping.evaluate(node => getComputedStyle(node).backgroundColor), band);
    const multi = group.locator(".read-hatched .op-M").first();
    assert.match(await multi.getAttribute("title"), /multimapping/);
    assert.match(await multi.evaluate(node => getComputedStyle(node).backgroundImage), /repeating-linear-gradient/);
    const gray = await group.locator(".read-muted .cigar-op").evaluateAll(nodes => nodes.map(node => getComputedStyle(node).backgroundColor));
    assert.ok(gray.length > 0 && gray.every(color => color === "rgb(146, 152, 158)"));
    const secondary = group.locator('.read-hit[title*="Flags: secondary"]');
    assert.match(await secondary.getAttribute("title"), /multimapping.*NH: 2/);
    assert.match(await secondary.locator(".op-M").first().getAttribute("title"), /multimapping.*Flags: secondary.*NH: 2/);
    await secondary.focus();
    await page.keyboard.press("Enter");
    assert.equal(await group.locator(".read-inspection[open]").count(), 1);
    // A separated mate remains clickable despite sharing a display lane.
    await separated.locator(".read-row").first().locator(".op-M").first().click();
    assert.equal(await separated.locator(".read-inspection[open]").count(), 1);
    assert.ok(await separated.evaluate(pair => ["::before", "::after"].every(side =>
      getComputedStyle(pair, side).display === "none")));
    await overlapping.locator(".op-M").first().click();
    assert.equal(await overlapping.locator(".read-inspection[open]").count(), 1);
    assert.ok(await overlapping.evaluate(pair => ["::before", "::after"].every(side =>
      getComputedStyle(pair, side).display === "none")));
    await overlapping.screenshot({path: process.argv[2] + ".details.png"});

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
    await group.getByRole("combobox", {name: "Read layout"}).selectOption("single");
    assert.equal(await group.locator(".read-pair").count(), 0);
    assert.equal(await group.locator(".read-pair-help").isVisible(), false);
    assert.equal(await group.locator(".read-row:visible").count(), 11);
    assert.equal(await group.locator(".read-row").count(), 129);
    assert.equal(await page.locator(".variant-detail table").innerText(), support);
    await group.locator(".read-filter").selectOption("all");
    await group.getByRole("combobox", {name: "Read layout"}).selectOption("pairs");
    await group.locator(".read-expand").click();
    await group.locator(".read-next").click();
    const brackets = await group.locator(".read-pair:visible").evaluateAll(pairs => pairs.map(pair => {
      const style = getComputedStyle(pair, "::before");
      const left = pair.getBoundingClientRect().left + parseFloat(style.left);
      const readLeft = Math.min(...Array.from(pair.querySelectorAll(".cigar-op"), op => op.getBoundingClientRect().left));
      const rightStyle = getComputedStyle(pair, "::after");
      const right = pair.getBoundingClientRect().left + parseFloat(rightStyle.left);
      const readRight = Math.max(...Array.from(pair.querySelectorAll(".cigar-op"), op => op.getBoundingClientRect().right));
      return {left, gap: readLeft - left - parseFloat(style.width), rightGap: right - readRight,
        rightBorder: parseFloat(rightStyle.borderRightWidth)};
    }));
    // Leave clearance for six-pixel arrows protruding beyond either endpoint.
    assert.ok(brackets.every(bracket => bracket.gap >= 8 && bracket.gap <= 12));
    assert.ok(brackets.every(bracket => bracket.rightGap >= 8 && bracket.rightGap <= 12 && bracket.rightBorder > 0));
    assert.ok(new Set(brackets.map(bracket => Math.round(bracket.left))).size > 1);
    await page.screenshot({path: process.argv[2] + ".png", fullPage: true});
    // Shared physical records retain role-specific strand labels after either
    // support filter, in both catalog orders and both layout modes.
    if (process.argv[3]) {
      await page.goto(pathToFileURL(process.argv[3]).href + "#variant-0");
      for (const index of [0, 1]) {
        await page.goto(pathToFileURL(process.argv[3]).href + `#variant-${index}`);
        const shared = page.locator(`#pileup-v${index}`);
        const counts = await page.locator(`#variant-${index} table`).innerText();
        for (const layout of ["single", "pairs"]) {
          await shared.locator(".read-pairs").selectOption(layout);
          for (const label of ["Same donor", "Same acceptor"]) {
            await shared.locator(".support-filter").selectOption({label});
            assert.equal(await shared.locator(".read-row:visible").count(), 2);
            assert.equal(await shared.locator(".read-pair:visible").count(), layout === "pairs" ? 1 : 0);
            for (const read of await shared.locator(".read-row:visible .read-hit").all()) {
              for (const attribute of ["title", "aria-label"]) {
                const text = await read.getAttribute(attribute);
                assert.match(text, /Same donor: strand evidence agrees/);
                assert.match(text, /Same acceptor: strand unverified/);
                assert.match(text, /multimapping/);
              }
            }
          }
        }
        assert.equal(await page.locator(`#variant-${index} table`).innerText(), counts);
      }
    }
    assert.ok(requests.every(url => url.startsWith("file:")));
    assert.deepEqual(errors, []);
    console.log("Paired-view browser checks passed.");
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
