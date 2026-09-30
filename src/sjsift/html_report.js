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
  update();
})();
