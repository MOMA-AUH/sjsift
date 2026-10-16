/* Progressive enhancement: the generated HTML already contains every count. */
(() => {
  "use strict";
  const byId = id => document.getElementById(id);
  const tbody = byId("overview-table").tBodies[0];
  const rows = Array.from(tbody.rows);
  const details = Array.from(document.querySelectorAll(".variant-detail"));
  const picker = byId("variant-picker");
  let visible = [];

  function showView(focus = false) {
    const match = /^#variant-(\d+)$/.exec(location.hash);
    const index = match ? match[1] : null;
    const selected = visible.find(row => row.dataset.index === index);
    const detail = selected ? byId(`variant-${index}`) : null;
    byId("overview").hidden = Boolean(detail);
    byId("details").hidden = !detail;
    for (const article of details) article.hidden = article !== detail;
    if (detail) {
      picker.value = index;
      const position = visible.indexOf(selected);
      byId("previous").disabled = position === 0;
      byId("next").disabled = position === visible.length - 1;
      byId("view-status").textContent = `Showing ${picker.selectedOptions[0].textContent}, variant ${position + 1} of ${visible.length}.`;
      if (focus) detail.querySelector("h2").focus();
    } else {
      byId("view-status").textContent = `${visible.length} variants match the filters.`;
      if (focus) byId("overview-heading").focus();
    }
  }

  function update() {
    const search = byId("search").value.toLowerCase();
    const observedOnly = byId("observed-only").checked;
    const supportOrder = byId("order").value === "support";
    const ordered = rows.slice().sort((a, b) => {
      if (supportOrder) {
        // Counts can exceed JavaScript's safe integer range; compare them exactly.
        for (const field of ["unique", "total"]) {
          const left = BigInt(a.dataset[field]);
          const right = BigInt(b.dataset[field]);
          if (left !== right) return left > right ? -1 : 1;
        }
      }
      return Number(a.dataset.index) - Number(b.dataset.index);
    });
    visible = [];
    picker.replaceChildren();
    for (const row of ordered) {
      const name = row.querySelector("a").textContent;
      row.hidden = !name.toLowerCase().includes(search) || (observedOnly && row.dataset.total === "0");
      tbody.append(row);
      if (!row.hidden) {
        visible.push(row);
        const option = document.createElement("option");
        option.value = row.dataset.index;
        option.textContent = name;
        picker.append(option);
      }
    }
    byId("visible-count").textContent = String(visible.length);
    byId("no-matches").hidden = visible.length !== 0;
    showView();
  }

  for (const id of ["controls", "visible-summary", "detail-navigation"]) byId(id).hidden = false;
  byId("search").addEventListener("input", update);
  byId("order").addEventListener("change", update);
  byId("observed-only").addEventListener("change", update);
  picker.addEventListener("change", () => { location.hash = `variant-${picker.value}`; });
  for (const [id, step] of [["previous", -1], ["next", 1]]) {
    byId(id).addEventListener("click", () => {
      const position = visible.findIndex(row => row.dataset.index === picker.value);
      const next = visible[position + step];
      if (next) location.hash = `variant-${next.dataset.index}`;
    });
  }
  window.addEventListener("hashchange", () => showView(true));
  for (const group of document.querySelectorAll(".alignment-group")) {
    const find = selector => group.querySelector(selector);
    const reads = Array.from(group.querySelectorAll(".read-row"));
    const tracks = Array.from(group.querySelectorAll(".read-rows"), list => ({
      list, rows: Array.from(list.children)
    }));
    // Keep inspection compact while retaining native, keyboard-accessible details.
    for (const row of reads) {
      const inspection = row.querySelector(".read-inspection");
      inspection.addEventListener("toggle", () => {
        if (!inspection.open) return;
        for (const other of reads) {
          if (other !== row) other.querySelector(".read-inspection").open = false;
        }
      });
    }
    const filter = find(".read-filter");
    const supportFilter = find(".support-filter");
    const pairs = find(".read-pairs");
    const expand = find(".read-expand");
    const previous = find(".read-previous");
    const next = find(".read-next");
    let expanded = false;
    let page = 0;
    function arrangeReads(matching) {
      // Move the existing records; never duplicate geometry or inspection controls.
      for (const track of tracks) {
        const fragment = document.createDocumentFragment();
        for (const row of track.rows) fragment.append(row);
        track.list.append(fragment);
      }
      for (const pair of group.querySelectorAll(".read-pair")) pair.remove();
      const mates = new Map();
      if (pairs.value === "pairs") {
        for (const row of matching) {
          if (!row.dataset.pair) continue;
          if (!mates.has(row.dataset.pair)) mates.set(row.dataset.pair, []);
          mates.get(row.dataset.pair).push(row);
        }
      }
      const used = new Set();
      const units = [];
      for (const row of matching) {
        if (used.has(row)) continue;
        const members = mates.get(row.dataset.pair);
        if (!members || members.length !== 2) {
          units.push({element: row, rows: [row], bucket: Number(row.dataset.bucket)});
          continue;
        }
        const wrapper = document.createElement("li");
        wrapper.className = "read-pair";
        wrapper.style.setProperty("--pair-left", `${Math.min(...members.map(member => Number(member.dataset.left)))}%`);
        wrapper.style.setProperty("--pair-right", `${Math.max(...members.map(member => Number(member.dataset.right)))}%`);
        const list = document.createElement("ol");
        list.className = "read-pair-members";
        list.setAttribute("aria-label", "Associated mates");
        const left = Math.min(...members.map(member => Number(member.dataset.right)));
        const right = Math.max(...members.map(member => Number(member.dataset.left)));
        if (left < right) {
          const connector = document.createElement("span");
          connector.className = "mate-connector";
          connector.setAttribute("aria-hidden", "true");
          connector.style.left = `${left}%`;
          connector.style.width = `${right - left}%`;
          wrapper.append(connector);
        } else {
          wrapper.classList.add("pair-overlap");
        }
        row.before(wrapper);
        list.append(...members);
        wrapper.append(list);
        members.forEach(member => used.add(member));
        units.push({element: wrapper, rows: members,
          bucket: Math.min(...members.map(member => Number(member.dataset.bucket)))});
      }
      return units;
    }
    function updatePairLanes() {
      for (const connector of group.querySelectorAll(".mate-connector")) {
        const pair = connector.parentElement;
        if (!pair.clientWidth) continue;
        const gap = pair.clientWidth * parseFloat(connector.style.width) / 100;
        const arrowSize = parseFloat(getComputedStyle(pair).getPropertyValue("--read-arrow-size"));
        // Even nonoverlapping alignments need separate lanes if their arrows touch.
        const stacked = gap <= 2 * arrowSize + 2;
        pair.classList.toggle("pair-overlap", stacked);
        connector.hidden = stacked;
      }
    }
    function updateReads() {
      find(".read-pair-help").hidden = pairs.value !== "pairs";
      const matching = reads.filter(row =>
        (filter.value === "all" || row.dataset.mapping === "unique") &&
        (supportFilter.value === "all" || row.dataset.support.split(" ").includes(supportFilter.value))
      );
      const units = arrangeReads(matching);
      const size = expanded ? 50 : 10;
      // Give each support group a preview so abundant reference reads cannot
      // push the variant evidence out of the initial figure (or vice versa).
      const buckets = new Map();
      for (const unit of units) {
        const key = unit.bucket;
        if (!buckets.has(key)) buckets.set(key, []);
        buckets.get(key).push(unit);
      }
      const pages = Math.max(1, ...Array.from(buckets.values(), rows => Math.ceil(rows.length / size)));
      page = Math.min(page, pages - 1);
      const shownUnits = Array.from(buckets.values()).flatMap(rows => rows.slice(page * size, (page + 1) * size));
      const shown = new Set(shownUnits.flatMap(unit => unit.rows));
      const visibleUnits = new Set(shownUnits);
      for (const unit of units) unit.element.hidden = !visibleUnits.has(unit);
      for (const row of reads) {
        row.hidden = !shown.has(row);
        if (row.hidden) row.querySelector(".read-inspection").open = false;
      }
      updatePairLanes();
      find(".read-filter-empty").hidden = matching.length !== 0;
      expand.textContent = expanded ? "Collapse reads" : "Expand reads";
      expand.setAttribute("aria-expanded", String(expanded));
      previous.hidden = next.hidden = !expanded;
      previous.disabled = page === 0;
      next.disabled = page === pages - 1;
      const pairedRows = shownUnits.filter(unit => unit.rows.length === 2).length;
      const layout = pairs.value === "pairs" ? ` · ${pairedRows} pairs and ${shownUnits.length - pairedRows} individual rows` : "";
      const unitLabel = pairs.value === "pairs" ? "rows" : "records";
      find(".read-page").textContent = !matching.length ? "No matching reads." : expanded
        ? `Page ${page + 1} of ${pages} · showing ${shown.size} of ${matching.length} matching records${layout} · up to 50 ${unitLabel} per support group${page === pages - 1 ? " · End of embedded reads" : ""}`
        : `Preview: ${shown.size} of ${matching.length} embedded records${layout} · up to 10 ${unitLabel} per support group. Expand for more.`;
    }
    filter.addEventListener("change", () => { page = 0; updateReads(); });
    supportFilter.addEventListener("change", () => { page = 0; updateReads(); });
    pairs.addEventListener("change", () => { page = 0; updateReads(); });
    expand.addEventListener("click", () => { expanded = !expanded; page = 0; updateReads(); });
    previous.addEventListener("click", () => { page -= 1; updateReads(); });
    next.addEventListener("click", () => { page += 1; updateReads(); });
    const resizeObserver = new ResizeObserver(updatePairLanes);
    for (const track of tracks) resizeObserver.observe(track.list);
    find(".read-controls").hidden = false;
    updateReads();
  }
  update();
})();
