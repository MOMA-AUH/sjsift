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
    const expand = find(".read-expand");
    const previous = find(".read-previous");
    const next = find(".read-next");
    let expanded = false;
    let page = 0;
    function updateReads() {
      const matching = reads.filter(row =>
        (filter.value === "all" || row.dataset.mapping === "unique") &&
        (supportFilter.value === "all" || row.dataset.support.split(" ").includes(supportFilter.value))
      );
      const size = expanded ? 50 : 10;
      // Give each support group a preview so abundant reference reads cannot
      // push the variant evidence out of the initial figure (or vice versa).
      const buckets = new Map();
      for (const row of matching) {
        const key = row.dataset.bucket;
        if (!buckets.has(key)) buckets.set(key, []);
        buckets.get(key).push(row);
      }
      const pages = Math.max(1, ...Array.from(buckets.values(), rows => Math.ceil(rows.length / size)));
      page = Math.min(page, pages - 1);
      const shown = new Set(Array.from(buckets.values()).flatMap(rows => rows.slice(page * size, (page + 1) * size)));
      for (const row of reads) {
        row.hidden = !shown.has(row);
        if (row.hidden) row.querySelector(".read-inspection").open = false;
      }
      find(".shown-count").textContent = `Currently shown: ${shown.size} (${matching.length} embedded match filter)`;
      find(".read-filter-empty").hidden = matching.length !== 0;
      expand.textContent = expanded ? "Collapse reads" : "Expand reads";
      expand.setAttribute("aria-expanded", String(expanded));
      previous.hidden = next.hidden = !expanded;
      previous.disabled = page === 0;
      next.disabled = page === pages - 1;
      find(".read-page").textContent = !matching.length ? "No matching reads." : expanded
        ? `Page ${page + 1} of ${pages} · showing ${shown.size} of ${matching.length} matching records · up to 50 per support group${page === pages - 1 ? " · End of embedded reads" : ""}`
        : `Preview: ${shown.size} of ${matching.length} embedded records · up to 10 per support group. Expand for more.`;
    }
    filter.addEventListener("change", () => { page = 0; updateReads(); });
    supportFilter.addEventListener("change", () => { page = 0; updateReads(); });
    expand.addEventListener("click", () => { expanded = !expanded; page = 0; updateReads(); });
    previous.addEventListener("click", () => { page -= 1; updateReads(); });
    next.addEventListener("click", () => { page += 1; updateReads(); });
    find(".read-controls").hidden = false;
    updateReads();
  }
  update();
})();
