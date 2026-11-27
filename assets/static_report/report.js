const rows = Array.from(document.querySelectorAll(".sample-row"));
const filters = Array.from(document.querySelectorAll("[data-filter-key]"));
const panels = Array.from(document.querySelectorAll(".details-panel:not([data-organism-panel])"));
const organismPanels = Array.from(document.querySelectorAll("[data-organism-panel]"));
const organismPanel = document.getElementById("organism-detail-panel");
const sampleSortHeaders = Array.from(document.querySelectorAll(".sample-sort"));
const sampleTableBody = document.querySelector(".panel .table-wrap table tbody");
let sampleSortState = { key: "sample", direction: "asc" };
const detailCache = new Map();

function escapeHtml(value) {
  return String(value ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/\"/g, "&quot;")
    .replace(/'/g, "&#39;");
}

function normaliseHeaderKey(header) {
  return String(header || "")
    .trim()
    .toLowerCase();
}

function advanceAlignmentCoord(start, segment, step) {
  let next = Number(start);
  let last = null;
  for (const character of String(segment || "")) {
    if (character !== "-") {
      last = next;
      next += step;
    }
  }
  return {
    last: last ?? Number(start),
    next,
  };
}

// BLASTX query previews use the standard genetic code unless alternate
// translation tables are added to the report metadata later.
const STANDARD_GENETIC_CODE = {
  TTT: "F",
  TTC: "F",
  TTA: "L",
  TTG: "L",
  TCT: "S",
  TCC: "S",
  TCA: "S",
  TCG: "S",
  TAT: "Y",
  TAC: "Y",
  TAA: "*",
  TAG: "*",
  TGT: "C",
  TGC: "C",
  TGA: "*",
  TGG: "W",
  CTT: "L",
  CTC: "L",
  CTA: "L",
  CTG: "L",
  CCT: "P",
  CCC: "P",
  CCA: "P",
  CCG: "P",
  CAT: "H",
  CAC: "H",
  CAA: "Q",
  CAG: "Q",
  CGT: "R",
  CGC: "R",
  CGA: "R",
  CGG: "R",
  ATT: "I",
  ATC: "I",
  ATA: "I",
  ATG: "M",
  ACT: "T",
  ACC: "T",
  ACA: "T",
  ACG: "T",
  AAT: "N",
  AAC: "N",
  AAA: "K",
  AAG: "K",
  AGT: "S",
  AGC: "S",
  AGA: "R",
  AGG: "R",
  GTT: "V",
  GTC: "V",
  GTA: "V",
  GTG: "V",
  GCT: "A",
  GCC: "A",
  GCA: "A",
  GCG: "A",
  GAT: "D",
  GAC: "D",
  GAA: "E",
  GAG: "E",
  GGT: "G",
  GGC: "G",
  GGA: "G",
  GGG: "G",
};
function translateCodon(codon, geneticCode = STANDARD_GENETIC_CODE) {
  return geneticCode[codon] || "X";
}

function translateBlastxQuery(sequence) {
  const cleaned = String(sequence || "")
    .toUpperCase()
    .replace(/\s+/g, "");
  let translated = "";
  let codon = "";

  for (const character of cleaned) {
    if (character === "-") {
      if (codon.length === 0) {
        translated += "-";
      }
      continue;
    }
    codon += character;
    if (codon.length === 3) {
      translated += translateCodon(codon);
      codon = "";
    }
  }

  return translated;
}

function isBlastxLikeAlignment(qseq, sseq) {
  const query = String(qseq || "");
  const subject = String(sseq || "");
  if (!query || !subject) {
    return false;
  }
  const aaLike = /^[A-Z*\-]+$/i.test(subject);
  return aaLike && query.length > subject.length * 2;
}

function buildProteinAlignment(queryAa, subjectAa, qstart, qend, sstart, send, width = 60) {
  const querySequence = String(queryAa || "");
  const subjectSequence = String(subjectAa || "");
  const queryNtStart = Number(qstart);
  const queryNtEnd = Number(qend);
  const subjectStart = Number(sstart);
  const subjectEnd = Number(send);

  if (!querySequence || !subjectSequence) {
    return '<div class="blast-align-empty">No translated alignment available for this hit.</div>';
  }

  const subjectStep = subjectEnd >= subjectStart ? 1 : -1;
  const labelWidth = 5;
  const coordinateWidth = Math.max(
    String(subjectStart).length,
    String(subjectEnd).length,
    String(Math.ceil(queryNtStart / 3)).length,
    String(Math.ceil(queryNtEnd / 3)).length,
    3,
  );
  const padLabel = (value) => String(value).padEnd(labelWidth, " ");
  const padCoordinate = (value) => String(value).padStart(coordinateWidth, " ");
  const lines = [`Query nt range: ${escapeHtml(`${queryNtStart}-${queryNtEnd}`)}`, ""];
  let queryAaCursor = 1;
  let subjectCursor = subjectStart;
  const maxLength = Math.max(querySequence.length, subjectSequence.length);

  for (let offset = 0; offset < maxLength; offset += width) {
    const queryChunk = querySequence.slice(offset, offset + width);
    const subjectChunk = subjectSequence.slice(offset, offset + width);
    const matchChunk = Array.from({ length: Math.max(queryChunk.length, subjectChunk.length) }, (_, index) => {
      const queryBase = queryChunk[index] || " ";
      const subjectBase = subjectChunk[index] || " ";
      return queryBase !== "-" && subjectBase !== "-" && queryBase.toUpperCase() === subjectBase.toUpperCase()
        ? "|"
        : " ";
    }).join("");

    const queryRange = advanceAlignmentCoord(queryAaCursor, queryChunk, 1);
    const subjectRange = advanceAlignmentCoord(subjectCursor, subjectChunk, subjectStep);

    lines.push(
      `${padLabel("Query")} ${padCoordinate(queryAaCursor)}  ${escapeHtml(queryChunk)}  ${padCoordinate(queryRange.last)}`,
      `${padLabel("")} ${" ".repeat(coordinateWidth)}  ${escapeHtml(matchChunk)}`,
      `${padLabel("Sbjct")} ${padCoordinate(subjectCursor)}  ${escapeHtml(subjectChunk)}  ${padCoordinate(subjectRange.last)}`,
    );

    queryAaCursor = queryRange.next;
    subjectCursor = subjectRange.next;

    if (offset + width < maxLength) {
      lines.push("");
    }
  }

  return `<pre class="blast-align-pre">${lines.join("\n")}</pre>`;
}

function buildBlastAlignment(rowMap, width = 60) {
  const qseq = String(rowMap.qseq || "");
  const sseq = String(rowMap.sseq || "");
  const qstart = Number(rowMap.qstart);
  const qend = Number(rowMap.qend);
  const sstart = Number(rowMap.sstart);
  const send = Number(rowMap.send);

  if (
    !qseq ||
    !sseq ||
    !Number.isFinite(qstart) ||
    !Number.isFinite(qend) ||
    !Number.isFinite(sstart) ||
    !Number.isFinite(send)
  ) {
    return '<div class="blast-align-empty">No alignment sequence available for this hit.</div>';
  }

  if (isBlastxLikeAlignment(qseq, sseq)) {
    return buildProteinAlignment(translateBlastxQuery(qseq), sseq, qstart, qend, sstart, send, width);
  }

  const queryStep = qend >= qstart ? 1 : -1;
  const subjectStep = send >= sstart ? 1 : -1;
  const labelWidth = 5;
  const coordinateWidth = Math.max(
    String(qstart).length,
    String(qend).length,
    String(sstart).length,
    String(send).length,
    3,
  );
  const padLabel = (value) => String(value).padEnd(labelWidth, " ");
  const padCoordinate = (value) => String(value).padStart(coordinateWidth, " ");
  const lines = [];
  let queryCursor = qstart;
  let subjectCursor = sstart;
  const maxLength = Math.max(qseq.length, sseq.length);

  for (let offset = 0; offset < maxLength; offset += width) {
    const queryChunk = qseq.slice(offset, offset + width);
    const subjectChunk = sseq.slice(offset, offset + width);
    const matchChunk = Array.from({ length: Math.max(queryChunk.length, subjectChunk.length) }, (_, index) => {
      const queryBase = queryChunk[index] || " ";
      const subjectBase = subjectChunk[index] || " ";
      return queryBase !== "-" && subjectBase !== "-" && queryBase.toUpperCase() === subjectBase.toUpperCase()
        ? "|"
        : " ";
    }).join("");

    const queryRange = advanceAlignmentCoord(queryCursor, queryChunk, queryStep);
    const subjectRange = advanceAlignmentCoord(subjectCursor, subjectChunk, subjectStep);

    lines.push(
      `${padLabel("Query")} ${padCoordinate(queryCursor)}  ${escapeHtml(queryChunk)}  ${padCoordinate(queryRange.last)}`,
      `${padLabel("")} ${" ".repeat(coordinateWidth)}  ${escapeHtml(matchChunk)}`,
      `${padLabel("Sbjct")} ${padCoordinate(subjectCursor)}  ${escapeHtml(subjectChunk)}  ${padCoordinate(subjectRange.last)}`,
    );

    queryCursor = queryRange.next;
    subjectCursor = subjectRange.next;

    if (offset + width < maxLength) {
      lines.push("");
    }
  }

  return `<pre class="blast-align-pre">${lines.join("\n")}</pre>`;
}

function detailSortCellIndex(visibleIndex, alignmentInsertIndex) {
  if (alignmentInsertIndex > 0 && visibleIndex >= alignmentInsertIndex) {
    return visibleIndex + 1;
  }
  return visibleIndex;
}

function renderDetailTable(headers, rows, options = {}) {
  if (!headers || !headers.length) {
    return options.emptyMessage ? `<div class="muted">${escapeHtml(options.emptyMessage)}</div>` : "";
  }

  const headerClass = options.sortable ? ' class="detail-sort"' : "";
  const tableClass = options.coverage ? ' class="coverage-help"' : "";
  const headerKeys = headers.map(normaliseHeaderKey);
  const hiddenColumns = options.blastAlignment
    ? new Set(["qstart", "qend", "sstart", "send", "qseq", "sseq"])
    : new Set();
  const alignmentKey = options.blastAlignmentKey || "blast";
  const visibleIndices = headers
    .map((header, index) => ({ header, index, key: headerKeys[index] }))
    .filter((item) => !hiddenColumns.has(item.key));
  const scientificNameIndex = visibleIndices.findIndex((item) => item.key === "sscinames");
  const alignmentInsertIndex = options.blastAlignment ? scientificNameIndex + 1 : -1;

  return `
        <div class="detail-table-wrap">
          <table${options.sortable ? " data-sortable-blast" : ""}${tableClass}${options.blastAlignment ? ' data-blast-alignment="true"' : ""}>
            <thead>
              <tr>
                ${visibleIndices
                  .map(
                    (item, visibleIndex) => `
                  <th${headerClass}${options.sortable ? ` data-sort-index="${detailSortCellIndex(visibleIndex, alignmentInsertIndex)}"` : ""} title="${escapeHtml(options.headerHelp?.[item.key] || item.header)}">${escapeHtml(item.header)}</th>
                  ${options.blastAlignment && visibleIndex + 1 === alignmentInsertIndex ? '<th title="alignment preview">alignment</th>' : ""}
                `,
                  )
                  .join("")}
                ${options.blastAlignment && alignmentInsertIndex === 0 ? '<th title="alignment preview">alignment</th>' : ""}
              </tr>
            </thead>
            <tbody>
              ${rows
                .map((row, rowIndex) => {
                  if (!options.blastAlignment) {
                    return `<tr>${row.map((cell) => `<td>${escapeHtml(cell)}</td>`).join("")}</tr>`;
                  }

                  const rowMap = Object.fromEntries(
                    headers.map((header, headerIndex) => [headerKeys[headerIndex], row[headerIndex] ?? ""]),
                  );
                  return `
                  <tr class="blast-summary-row">
                    ${visibleIndices
                      .map(
                        (item, visibleIndex) => `
                      <td>${escapeHtml(row[item.index])}</td>
                      ${
                        visibleIndex + 1 === alignmentInsertIndex
                          ? `
                        <td>
                          <button
                            type="button"
                            class="blast-align-button"
                            data-blast-align-toggle
                            data-blast-align-target="${alignmentKey}-align-${rowIndex}"
                            aria-expanded="false"
                          >View alignment</button>
                        </td>
                      `
                          : ""
                      }
                    `,
                      )
                      .join("")}
                    ${
                      alignmentInsertIndex === 0
                        ? `
                      <td>
                        <button
                          type="button"
                          class="blast-align-button"
                          data-blast-align-toggle
                          data-blast-align-target="${alignmentKey}-align-${rowIndex}"
                          aria-expanded="false"
                        >View alignment</button>
                      </td>
                    `
                        : ""
                    }
                  </tr>
                  <tr class="blast-align-row" id="${alignmentKey}-align-${rowIndex}" hidden>
                    <td colspan="${visibleIndices.length + 1}">
                      <div class="blast-align-panel">${buildBlastAlignment(rowMap)}</div>
                    </td>
                  </tr>
                `;
                })
                .join("")}
            </tbody>
          </table>
        </div>
      `;
}

function renderReadsSections(readsSections) {
  if (!readsSections || !readsSections.length) {
    return '<div class="muted">No extracted reads found.</div>';
  }

  return readsSections
    .map(
      (section) => `
        <div class="subsection">
          <div class="reads-toolbar">
            <button type="button" class="reads-action" data-select-all-reads>Select all</button>
            <button type="button" class="reads-action" data-copy-reads>Copy selected reads</button>
            <button type="button" class="reads-action" data-blast-reads>NCBI BLAST</button>
          </div>
          <div class="reads-list">
            ${(section.entries || [])
              .map(
                (entry) => `
              <div class="read-entry">
                <label class="read-entry-label">
                  <input type="checkbox" data-read-checkbox />
                  <span>${escapeHtml(entry.header)}</span>
                </label>
                <pre data-read-content>${escapeHtml(entry.content)}</pre>
              </div>
            `,
              )
              .join("")}
          </div>
        </div>
      `,
    )
    .join("");
}

function decodeBase64ToBytes(base64Text) {
  if (typeof Uint8Array.fromBase64 === "function") {
    return Uint8Array.fromBase64(base64Text);
  }
  const binary = atob(base64Text);
  const bytes = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index += 1) {
    bytes[index] = binary.charCodeAt(index);
  }
  return bytes;
}

async function gunzipBytes(bytes) {
  if (!("DecompressionStream" in window)) {
    throw new Error("This browser does not support gzip detail payloads.");
  }
  const response = new Response(new Response(bytes).body.pipeThrough(new DecompressionStream("gzip")));
  return response.text();
}

async function gunzipBase64(text) {
  return gunzipBytes(decodeBase64ToBytes(text.trim()));
}

async function loadDetailData(panelId) {
  if (detailCache.has(panelId)) {
    return detailCache.get(panelId);
  }
  const node = document.getElementById(`detail-data-${panelId}`);
  if (!node) {
    return null;
  }
  const detail = JSON.parse(await gunzipBase64(node.textContent));
  const plots = Array.from(document.querySelectorAll(`script[data-plot-for="${CSS.escape(panelId)}"]`));
  (detail.coverage_plot_sections || []).forEach((section, index) => {
    section.src = plots[index] ? plots[index].textContent.trim() : "";
  });
  detailCache.set(panelId, detail);
  return detail;
}

function bindDetailPanelActions(panel) {
  panel.querySelectorAll("[data-sortable-blast]").forEach((table) => {
    initSortableBlastTable(table);
  });

  panel.querySelectorAll("[data-blast-align-toggle]").forEach((button) => {
    button.addEventListener("click", () => {
      const target = panel.querySelector(`#${CSS.escape(button.dataset.blastAlignTarget || "")}`);
      if (!target) {
        return;
      }
      const isExpanded = button.getAttribute("aria-expanded") === "true";
      target.hidden = isExpanded;
      button.setAttribute("aria-expanded", isExpanded ? "false" : "true");
      button.textContent = isExpanded ? "View alignment" : "Hide alignment";
    });
  });

  panel.querySelectorAll("[data-copy-reads]").forEach((button) => {
    button.addEventListener("click", async () => {
      const section = button.closest(".subsection");
      if (!section) {
        return;
      }
      await copyText(getSelectedReads(section));
      button.textContent = "Copied";
      window.setTimeout(() => {
        button.textContent = "Copy selected reads";
      }, 1200);
    });
  });

  panel.querySelectorAll("[data-blast-reads]").forEach((button) => {
    button.addEventListener("click", () => {
      const section = button.closest(".subsection");
      if (!section) {
        return;
      }
      const reads = getSelectedReads(section);
      if (reads) {
        openBlastWithReads(reads);
      }
    });
  });

  panel.querySelectorAll("[data-select-all-reads]").forEach((button) => {
    button.addEventListener("click", () => {
      const section = button.closest(".subsection");
      if (!section) {
        return;
      }
      section.querySelectorAll("[data-read-checkbox]").forEach((checkbox) => {
        checkbox.checked = true;
      });
    });
  });

  panel.querySelectorAll("[data-close-organism-panel]").forEach((button) => {
    button.addEventListener("click", () => {
      openPanel(button.dataset.sampleId, button.dataset.classifier, false, "auto");
    });
  });
}

function renderOrganismPanel(detail) {
  if (!organismPanel || !detail) {
    return;
  }

  organismPanel.dataset.sampleId = detail.sample_id;
  organismPanel.dataset.classifier = detail.classifier_key;
  organismPanel.dataset.panelId = detail.panel_id;

  const root = organismPanel.querySelector("[data-organism-detail-root]");
  if (!root) {
    return;
  }

  const coveragePlots =
    detail.coverage_plot_sections && detail.coverage_plot_sections.length
      ? `
            <div class="coverage-plots-scroll">
              ${detail.coverage_plot_sections
                .map(
                  (section) => `
              <div class="plot">
                <img src="${escapeHtml(section.src)}" alt="${escapeHtml(section.name || section.src)}" />
              </div>
              `,
                )
                .join("")}
            </div>
          `
      : '<div class="muted">No coverage plots found.</div>';

  const coverageTables =
    detail.coverage_tables && detail.coverage_tables.length
      ? detail.coverage_tables
          .map(
            (table) => `
            <div class="subsection">
              ${renderDetailTable(table.headers || [], table.rows || [], {
                coverage: true,
                headerHelp: detailHeaderHelp,
              })}
            </div>
          `,
          )
          .join("")
      : '<div class="muted">No coverage stats found.</div>';

  root.innerHTML = `
        <div class="organism-detail-header">
          <div class="organism-back-bar">
            <button
              type="button"
              class="organism-back"
              data-close-organism-panel="${escapeHtml(detail.panel_id)}"
              data-sample-id="${escapeHtml(detail.sample_id)}"
              data-classifier="${escapeHtml(detail.classifier_key)}"
            >&#8592; Back</button>
            <span class="detail-separator">/</span>
            <span class="detail-sample">${escapeHtml(detail.sample)}</span>
          </div>
          <div class="organism-eyebrow">Metaval details</div>
          <h2 class="organism-title">${escapeHtml(detail.organism)}</h2>
          <div class="organism-meta">
            <div>Sample: ${escapeHtml(detail.sample)}</div>
            <div>Taxid: ${escapeHtml(detail.taxid)}</div>
            <div>Assigned reads: ${escapeHtml(detail.assigned_reads || "NA")}</div>
            <div>Classifier: ${escapeHtml(detail.classifier)}</div>
          </div>
        </div>

        <div class="organism-sections">
          <details open>
            <summary>BLAST section</summary>
            <div class="section-body">
              <div class="subsection">
                <h3>blastn</h3>
                ${renderDetailTable((detail.blastn || {}).headers || [], (detail.blastn || {}).rows || [], {
                  sortable: true,
                  blastAlignment: true,
                  blastAlignmentKey: "blastn",
                  emptyMessage: "No blastn results found.",
                  headerHelp: detailHeaderHelp,
                })}
              </div>
              <div class="subsection">
                <h3>blastx</h3>
                ${renderDetailTable((detail.blastx || {}).headers || [], (detail.blastx || {}).rows || [], {
                  sortable: true,
                  blastAlignment: true,
                  blastAlignmentKey: "blastx",
                  emptyMessage: "No blastx results found.",
                  headerHelp: detailHeaderHelp,
                })}
              </div>
            </div>
          </details>

          <details open>
            <summary>Mapping section</summary>
            <div class="section-body">
              <div class="mapping-section-scroll">${coverageTables}</div>
            </div>
          </details>

          <details open>
            <summary>Coverage plots</summary>
            <div class="section-body">${coveragePlots}</div>
          </details>

          <details open>
            <summary>Reads section</summary>
            <div class="section-body">${renderReadsSections(detail.reads_sections || [])}</div>
          </details>
        </div>
      `;

  bindDetailPanelActions(organismPanel);
}

/**
 * Populate a sample-level filter select with unique values from the sample rows.
 */
function populateSampleFilterOptions(select) {
  const key = select.dataset.filterKey;
  const values = [...new Set(rows.map((row) => row.dataset[key]).filter(Boolean))].sort();

  for (const value of values) {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = key === "is_ntc" ? (value === "true" ? "NTC" : "Sample") : value;
    select.appendChild(option);
  }
}

function applyFilters() {
  rows.forEach((row) => {
    const visible = filters.every((select) => {
      return !select.value || row.dataset[select.dataset.filterKey] === select.value;
    });
    row.style.display = visible ? "" : "none";
  });

  const activeRow = rows.find((row) => row.classList.contains("is-active") && row.style.display !== "none");
  if (!activeRow) {
    rows.forEach((row) => row.classList.remove("is-active"));
    panels.forEach((panel) => {
      panel.hidden = true;
    });
    organismPanels.forEach((panel) => {
      panel.hidden = true;
    });
  }
}

function compareSampleRows(left, right, key, direction) {
  const leftValue = (left.dataset[key] || "").trim().toLowerCase();
  const rightValue = (right.dataset[key] || "").trim().toLowerCase();
  const order = direction === "asc" ? 1 : -1;
  return leftValue.localeCompare(rightValue, undefined, { numeric: true, sensitivity: "base" }) * order;
}

function updateSampleSortHeaders() {
  sampleSortHeaders.forEach((header) => {
    header.classList.toggle(
      "is-asc",
      header.dataset.sortKey === sampleSortState.key && sampleSortState.direction === "asc",
    );
    header.classList.toggle(
      "is-desc",
      header.dataset.sortKey === sampleSortState.key && sampleSortState.direction === "desc",
    );
  });
}

function sortSampleRows(key) {
  sampleSortState = {
    key,
    direction: sampleSortState.key === key && sampleSortState.direction === "asc" ? "desc" : "asc",
  };
  rows
    .slice()
    .sort((left, right) => compareSampleRows(left, right, sampleSortState.key, sampleSortState.direction))
    .forEach((row) => sampleTableBody.appendChild(row));
  updateSampleSortHeaders();
  applyFilters();
}

function setClassifier(panel, classifier) {
  panel.querySelectorAll("[data-classifier-tab]").forEach((button) => {
    button.classList.toggle("is-active", button.dataset.classifierTab === classifier);
  });
  panel.querySelectorAll("[data-classifier-view]").forEach((view) => {
    view.hidden = view.dataset.classifierView !== classifier;
  });
}

function closePanel(sampleId) {
  const panel = document.getElementById(sampleId);
  const row = rows.find((item) => item.dataset.sampleId === sampleId);
  if (panel) {
    panel.hidden = true;
  }
  organismPanels.forEach((item) => {
    item.hidden = true;
  });
  if (row) {
    row.classList.remove("is-active");
    row.scrollIntoView({ behavior: "smooth", block: "center" });
  }
}

// Taxpasta tables are embedded as JSON and rendered on demand, a page of rows at a time,
// so the browser never holds every table of every sample in the DOM at once.
const TAXPASTA_PAGE_SIZE = 200;
const CLASSIFIER_LABEL_TITLES = { K: "Kraken2", C: "Centrifuge", D: "Diamond" };
const taxpastaStates = new WeakMap();

function taxpastaCellText(cell) {
  return Array.isArray(cell) ? cell.map((pair) => pair.join(": ")).join(" ") : String(cell ?? "");
}

const taxpastaLoading = new WeakMap();

function getTaxpastaState(view) {
  return taxpastaStates.get(view) || { data: null, matches: [], rendered: 0, openLineages: new Set() };
}

// Tables are stored gzip-compressed and unpacked the first time they are shown.
function loadTaxpastaState(view) {
  if (taxpastaStates.has(view)) {
    return Promise.resolve(taxpastaStates.get(view));
  }
  if (!taxpastaLoading.has(view)) {
    const node = view.querySelector("script[data-taxpasta-gz]");
    const promise = (node ? gunzipBase64(node.textContent).then(JSON.parse) : Promise.resolve(null)).then((data) => {
      const state = { data, matches: [], rendered: 0, openLineages: new Set(), searchText: null };
      taxpastaStates.set(view, state);
      return state;
    });
    taxpastaLoading.set(view, promise);
  }
  return taxpastaLoading.get(view);
}

function taxpastaSearchText(state) {
  if (!state.searchText) {
    state.searchText = state.data.rows.map((row) => row.map(taxpastaCellText).join(" ").toLowerCase());
  }
  return state.searchText;
}

function renderTaxpastaLineageRow(data, index) {
  const lineage = data.lineages[index];
  return `<tr class="taxpasta-lineage-row" data-lineage-for="${index}"><td colspan="${data.headers.length + 1}"><div class="lineage-text"><strong>lineage:</strong> ${escapeHtml(lineage || "No lineage available.")}</div></td></tr>`;
}

function renderTaxpastaRow(state, index) {
  const { data } = state;
  const row = data.rows[index];
  const isOpen = state.openLineages.has(index);
  const lineageCell = data.lineages[index]
    ? `<button type="button" class="lineage-toggle${isOpen ? " is-open" : ""}" data-lineage-index="${index}">${isOpen ? "Hide" : "Show"}</button>`
    : "";
  const cells = row
    .map((cell, cellIndex) => {
      if (data.headers[cellIndex] === "classifiers" && Array.isArray(cell)) {
        const lines = cell
          .map(
            ([label, value]) =>
              `<div class="taxonomy-classifier-line"><span class="taxonomy-classifier-label" title="${CLASSIFIER_LABEL_TITLES[label] || ""}">${escapeHtml(label)}:</span><span class="taxonomy-classifier-value">${escapeHtml(value)}</span></div>`,
          )
          .join("");
        return `<td class="taxonomy-classifiers">${lines}</td>`;
      }
      const text = escapeHtml(cell);
      return `<td title="${text}">${text}</td>`;
    })
    .join("");
  const panelId = escapeHtml((data.metaval_panel_ids || [])[index] || "");
  const main = `<tr class="taxonomy-row" data-row-index="${index}" data-metaval-panel="${panelId}"><td>${lineageCell}</td>${cells}</tr>`;
  return isOpen ? main + renderTaxpastaLineageRow(data, index) : main;
}

function renderMoreTaxpastaRows(view) {
  const state = getTaxpastaState(view);
  const body = view.querySelector("[data-taxpasta-body]");
  if (!state.data || !body) {
    return;
  }
  const next = state.matches.slice(state.rendered, state.rendered + TAXPASTA_PAGE_SIZE);
  body.insertAdjacentHTML("beforeend", next.map((index) => renderTaxpastaRow(state, index)).join(""));
  state.rendered += next.length;

  const more = view.querySelector("[data-taxpasta-more]");
  if (more) {
    const remaining = state.matches.length - state.rendered;
    more.hidden = remaining <= 0;
    const label = more.querySelector("[data-taxpasta-more-label]");
    if (label) {
      label.textContent = `Showing ${state.rendered.toLocaleString()} of ${state.matches.length.toLocaleString()} rows`;
    }
    const button = more.querySelector("[data-taxpasta-more-button]");
    if (button) {
      button.textContent = `Show ${Math.min(remaining, TAXPASTA_PAGE_SIZE).toLocaleString()} more`;
    }
  }
}

async function applyTaxpastaFilter(view) {
  if (!view) {
    return;
  }
  const state = await loadTaxpastaState(view);
  const body = view.querySelector("[data-taxpasta-body]");
  if (!state.data || !body) {
    return;
  }

  const filter = view.querySelector("[data-superkingdom-filter]");
  const flagFilter = view.querySelector("[data-flag-filter]");
  const search = view.querySelector("[data-organism-search]");
  const metavalOnlyButton = view.querySelector("[data-metaval-only]");
  const counter = view.querySelector("[data-organisms-shown]");
  const superkingdom = filter ? filter.value : "all";
  const flag = flagFilter ? flagFilter.value : "all";
  const searchValue = search ? search.value.trim().toLowerCase() : "";
  const metavalOnly = metavalOnlyButton ? metavalOnlyButton.dataset.active === "true" : false;
  const { data } = state;
  const searchText = searchValue ? taxpastaSearchText(state) : null;

  const matches = [];
  for (let index = 0; index < data.rows.length; index += 1) {
    if (superkingdom !== "all" && data.superkingdoms[index] !== superkingdom) continue;
    if (flag !== "all" && data.flags[index] !== flag) continue;
    if (metavalOnly && !data.metaval_checked[index]) continue;
    if (searchText && !searchText[index].includes(searchValue)) continue;
    matches.push(index);
  }
  // As before, lineages of rows that are filtered out are closed.
  const matchSet = new Set(matches);
  state.openLineages.forEach((index) => {
    if (!matchSet.has(index)) state.openLineages.delete(index);
  });

  state.matches = matches;
  state.rendered = 0;
  body.textContent = "";
  renderMoreTaxpastaRows(view);
  const wrap = view.querySelector(".taxpasta-wrap");
  if (wrap) {
    wrap.scrollTop = 0;
  }
  if (counter) {
    counter.textContent = String(matches.length);
  }
}

function toggleTaxpastaLineage(button) {
  const view = button.closest(".taxpasta-view");
  const row = button.closest(".taxonomy-row");
  if (!view || !row) {
    return;
  }
  const state = getTaxpastaState(view);
  const index = Number(button.dataset.lineageIndex);
  const next = row.nextElementSibling;
  if (state.openLineages.has(index)) {
    state.openLineages.delete(index);
    if (next && next.classList.contains("taxpasta-lineage-row")) {
      next.remove();
    }
    button.textContent = "Show";
    button.classList.remove("is-open");
  } else {
    state.openLineages.add(index);
    row.insertAdjacentHTML("afterend", renderTaxpastaLineageRow(state.data, index));
    button.textContent = "Hide";
    button.classList.add("is-open");
  }
}

function ensureTaxpastaRendered(view) {
  if (!view) {
    return;
  }
  const state = taxpastaStates.get(view);
  if (!state || state.rendered === 0) {
    applyTaxpastaFilter(view);
  }
}

function debounce(fn, waitMs) {
  let timer = null;
  return (...args) => {
    window.clearTimeout(timer);
    timer = window.setTimeout(() => fn(...args), waitMs);
  };
}

function defaultTaxonomyWidth(headerName, measuredWidth) {
  let initialWidth = Math.max(measuredWidth + 16, 96);
  if (headerName === "name") {
    initialWidth = 190;
  } else if (headerName === "classifiers") {
    initialWidth = 118;
  } else if (headerName === "lineage") {
    initialWidth = 88;
  } else if (headerName === "taxid") {
    initialWidth = 90;
  } else if (headerName === "rank") {
    initialWidth = 116;
  } else if (headerName.includes("_vs_")) {
    initialWidth = 170;
  } else if (headerName.includes("_ntc")) {
    initialWidth = 140;
  } else if (headerName.startsWith("srr")) {
    initialWidth = 140;
  } else if (headerName.includes("-dna") || headerName.includes("-rna")) {
    initialWidth = 140;
  }
  return initialWidth;
}

function syncTaxonomyTableWidth(table) {
  const cols = Array.from(table.querySelectorAll("col"));
  const totalWidth = cols.reduce((sum, col) => {
    return sum + parseFloat(col.style.width || 0);
  }, 0);
  table.style.width = `${Math.max(totalWidth, table.parentElement.clientWidth)}px`;
}

function initResizableTable(table) {
  if (!table || table.dataset.resizableReady === "true") {
    return;
  }

  const cols = Array.from(table.querySelectorAll("col"));
  const headers = Array.from(table.querySelectorAll("thead th"));
  if (!cols.length || !headers.length) {
    return;
  }

  headers.forEach((header, index) => {
    const col = cols[index];
    if (!col) {
      return;
    }

    const headerName = header.textContent.trim().toLowerCase();
    const initialWidth = defaultTaxonomyWidth(headerName, header.getBoundingClientRect().width);
    col.style.width = `${initialWidth}px`;

    const handle = header.querySelector(".col-resizer");
    if (!handle) {
      return;
    }

    handle.addEventListener("mousedown", (event) => {
      event.preventDefault();
      const startX = event.clientX;
      const startWidth = parseFloat(col.style.width || initialWidth);

      document.body.classList.add("is-resizing");

      function onMouseMove(moveEvent) {
        const nextWidth = Math.max(80, startWidth + moveEvent.clientX - startX);
        col.style.width = `${nextWidth}px`;
        syncTaxonomyTableWidth(table);
      }

      function onMouseUp() {
        document.body.classList.remove("is-resizing");
        window.removeEventListener("mousemove", onMouseMove);
        window.removeEventListener("mouseup", onMouseUp);
      }

      window.addEventListener("mousemove", onMouseMove);
      window.addEventListener("mouseup", onMouseUp);
    });
  });

  syncTaxonomyTableWidth(table);
  table.dataset.resizableReady = "true";
}

function scrollPanelIntoView(panel, classifier, behavior = "smooth") {
  const classifierView = panel.querySelector(`[data-classifier-view="${classifier}"]`);
  const target = classifierView || panel;
  window.requestAnimationFrame(() => {
    target.scrollIntoView({ behavior, block: "start" });
  });
}

function openPanel(sampleId, classifier = "kraken2", allowToggle = true, scrollBehavior = "smooth") {
  const panel = document.getElementById(sampleId);
  const activeRow = rows.find((row) => row.classList.contains("is-active"));
  const isSameRow = activeRow && activeRow.dataset.sampleId === sampleId;

  if (allowToggle && isSameRow && panel && !panel.hidden) {
    closePanel(sampleId);
    return;
  }

  panels.forEach((panel) => {
    panel.hidden = panel.id !== sampleId;
  });
  organismPanels.forEach((panel) => {
    panel.hidden = true;
  });
  rows.forEach((row) => {
    row.classList.toggle("is-active", row.dataset.sampleId === sampleId);
  });

  if (!panel) {
    return;
  }

  setClassifier(panel, classifier);
  panel.querySelectorAll(".taxpasta-table").forEach((table) => initResizableTable(table));
  ensureTaxpastaRendered(panel.querySelector(`[data-classifier-view="${classifier}"]`));
  if (scrollBehavior) {
    scrollPanelIntoView(panel, classifier, scrollBehavior);
  }
}

async function openOrganismPanel(panelId, scrollBehavior = "smooth") {
  const detail = await loadDetailData(panelId);
  if (!organismPanel || !detail) {
    return;
  }

  renderOrganismPanel(detail);
  panels.forEach((item) => {
    item.hidden = true;
  });
  organismPanels.forEach((item) => {
    item.hidden = item !== organismPanel;
  });
  rows.forEach((row) => {
    row.classList.toggle("is-active", row.dataset.sampleId === detail.sample_id);
  });

  if (scrollBehavior) {
    window.requestAnimationFrame(() => {
      organismPanel.scrollIntoView({ behavior: scrollBehavior, block: "start" });
    });
  }
}

function sortDetailTable(table, index, direction) {
  const tbody = table.tBodies[0];
  if (!tbody) {
    return;
  }
  const order = direction === "asc" ? 1 : -1;

  if (table.dataset.blastAlignment === "true") {
    Array.from(tbody.querySelectorAll(".blast-summary-row"))
      .map((summaryRow) => ({
        summaryRow,
        detailRow: summaryRow.nextElementSibling?.classList.contains("blast-align-row")
          ? summaryRow.nextElementSibling
          : null,
      }))
      .sort((left, right) => {
        const leftValue = (left.summaryRow.cells[index]?.textContent || "").trim();
        const rightValue = (right.summaryRow.cells[index]?.textContent || "").trim();
        return leftValue.localeCompare(rightValue, undefined, { numeric: true, sensitivity: "base" }) * order;
      })
      .forEach(({ summaryRow, detailRow }) => {
        tbody.appendChild(summaryRow);
        if (detailRow) {
          tbody.appendChild(detailRow);
        }
      });
    return;
  }

  Array.from(tbody.rows)
    .sort((left, right) => {
      const leftValue = (left.cells[index]?.textContent || "").trim();
      const rightValue = (right.cells[index]?.textContent || "").trim();
      return leftValue.localeCompare(rightValue, undefined, { numeric: true, sensitivity: "base" }) * order;
    })
    .forEach((row) => tbody.appendChild(row));
}

function initSortableBlastTable(table) {
  if (!table || table.dataset.sortableReady === "true") {
    return;
  }
  table.querySelectorAll(".detail-sort").forEach((header) => {
    header.addEventListener("click", () => {
      const direction = header.dataset.direction === "asc" ? "desc" : "asc";
      table.querySelectorAll(".detail-sort").forEach((item) => {
        item.classList.remove("is-asc", "is-desc");
        delete item.dataset.direction;
      });
      header.dataset.direction = direction;
      header.classList.add(direction === "asc" ? "is-asc" : "is-desc");
      sortDetailTable(table, Number(header.dataset.sortIndex), direction);
    });
  });
  table.dataset.sortableReady = "true";
}

async function copyText(text) {
  if (navigator.clipboard && window.isSecureContext) {
    await navigator.clipboard.writeText(text);
    return;
  }

  const area = document.createElement("textarea");
  area.value = text;
  area.setAttribute("readonly", "");
  area.style.position = "absolute";
  area.style.left = "-9999px";
  document.body.appendChild(area);
  area.select();
  document.execCommand("copy");
  area.remove();
}

function openBlastWithReads(reads) {
  const form = document.createElement("form");
  form.method = "post";
  form.action = "https://blast.ncbi.nlm.nih.gov/Blast.cgi";
  form.target = "_blank";
  form.style.display = "none";

  const fields = {
    PROGRAM: "blastn",
    PAGE_TYPE: "BlastSearch",
    LINK_LOC: "blasthome",
    QUERY: reads,
  };

  Object.entries(fields).forEach(([name, value]) => {
    const input = document.createElement(name === "QUERY" ? "textarea" : "input");
    input.name = name;
    input.value = value;
    form.appendChild(input);
  });

  document.body.appendChild(form);
  form.submit();
  form.remove();
}

function getSelectedReads(section) {
  const selected = Array.from(section.querySelectorAll("[data-read-checkbox]:checked"))
    .map((checkbox) => checkbox.closest(".read-entry")?.querySelector("[data-read-content]")?.textContent || "")
    .filter(Boolean);
  if (selected.length) {
    return selected.join("\n");
  }
  return Array.from(section.querySelectorAll("[data-read-content]"))
    .map((block) => block.textContent || "")
    .filter(Boolean)
    .join("\n");
}

function restorePanelFromUrl() {
  const hashId = window.location.hash.startsWith("#") ? window.location.hash.slice(1) : "";
  if (hashId.startsWith("detail-")) {
    openOrganismPanel(hashId, "auto");
    return;
  }
  const params = new URLSearchParams(window.location.search);
  const hashSampleId = window.location.hash.startsWith("#sample-") ? window.location.hash.slice(1) : "";
  const sampleId = params.get("sample") || hashSampleId;
  const classifier = params.get("classifier") || "kraken2";
  if (!sampleId) {
    return;
  }
  openPanel(sampleId, classifier, false, "auto");
  const panel = document.getElementById(sampleId);
  if (panel) {
    window.setTimeout(() => {
      scrollPanelIntoView(panel, classifier, "auto");
    }, 0);
  }
}

filters.forEach((select) => {
  populateSampleFilterOptions(select);
  select.addEventListener("change", applyFilters);
});

sampleSortHeaders.forEach((header) => {
  header.addEventListener("click", () => sortSampleRows(header.dataset.sortKey));
});

rows.forEach((row) => {
  row.addEventListener("click", () => openPanel(row.dataset.sampleId));
});

panels.forEach((panel) => {
  panel.querySelectorAll("[data-classifier-tab]").forEach((button) => {
    button.addEventListener("click", () => {
      setClassifier(panel, button.dataset.classifierTab);
      panel.querySelectorAll(".taxpasta-table").forEach((table) => initResizableTable(table));
      ensureTaxpastaRendered(panel.querySelector(`[data-classifier-view="${button.dataset.classifierTab}"]`));
    });
  });
});

document.querySelectorAll("[data-close-panel]").forEach((button) => {
  button.addEventListener("click", () => closePanel(button.dataset.closePanel));
});

document.querySelectorAll("[data-open-organism]").forEach((button) => {
  button.addEventListener("click", () => {
    openOrganismPanel(button.dataset.openOrganism);
  });
});

document.querySelectorAll("[data-toggle-organisms]").forEach((button) => {
  button.addEventListener("click", () => {
    const box = button.closest(".extracted-reads");
    if (!box) {
      return;
    }
    const hiddenItems = Array.from(box.querySelectorAll("[data-extra-organism]"));
    const expand = hiddenItems.some((item) => item.hidden);
    hiddenItems.forEach((item) => {
      item.hidden = !expand;
    });
    button.textContent = expand ? button.dataset.expandedLabel : button.dataset.collapsedLabel;
  });
});

document.querySelectorAll(".taxpasta-view").forEach((view) => {
  const filter = view.querySelector("[data-superkingdom-filter]");
  const flagFilter = view.querySelector("[data-flag-filter]");
  const search = view.querySelector("[data-organism-search]");
  const metavalOnlyButton = view.querySelector("[data-metaval-only]");
  if (filter) {
    filter.addEventListener("change", () => applyTaxpastaFilter(view));
  }
  if (flagFilter) {
    flagFilter.addEventListener("change", () => applyTaxpastaFilter(view));
  }
  if (search) {
    search.addEventListener("input", debounce(() => applyTaxpastaFilter(view), 150));
  }
  if (metavalOnlyButton) {
    metavalOnlyButton.dataset.active = "false";
    metavalOnlyButton.addEventListener("click", () => {
      const nextState = metavalOnlyButton.dataset.active !== "true";
      metavalOnlyButton.dataset.active = nextState ? "true" : "false";
      metavalOnlyButton.classList.toggle("is-active", nextState);
      metavalOnlyButton.textContent = nextState ? "On" : "Off";
      applyTaxpastaFilter(view);
    });
  }
  // Rows are rendered when the table is first shown (openPanel / tab switch), not here.
  const wrap = view.querySelector(".taxpasta-wrap");
  if (wrap) {
    wrap.addEventListener("scroll", () => {
      const state = taxpastaStates.get(view);
      if (state && state.rendered < state.matches.length && wrap.scrollTop + wrap.clientHeight >= wrap.scrollHeight - 200) {
        renderMoreTaxpastaRows(view);
      }
    }, { passive: true });
  }
});

// One delegated listener replaces per-row listeners on every taxpasta row.
document.addEventListener("click", (event) => {
  const lineageButton = event.target.closest(".lineage-toggle");
  if (lineageButton) {
    toggleTaxpastaLineage(lineageButton);
    return;
  }
  const moreButton = event.target.closest("[data-taxpasta-more-button]");
  if (moreButton) {
    renderMoreTaxpastaRows(moreButton.closest(".taxpasta-view"));
    return;
  }
  const taxonomyRow = event.target.closest(".taxonomy-row");
  if (taxonomyRow && taxonomyRow.dataset.metavalPanel) {
    openOrganismPanel(taxonomyRow.dataset.metavalPanel);
  }
});

restorePanelFromUrl();
updateSampleSortHeaders();
