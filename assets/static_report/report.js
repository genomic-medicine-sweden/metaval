const rows = Array.from(document.querySelectorAll(".sample-row"));
const filters = Array.from(document.querySelectorAll("[data-filter-key]"));
const panels = Array.from(document.querySelectorAll(".details-panel:not([data-organism-panel])"));
const organismPanels = Array.from(document.querySelectorAll("[data-organism-panel]"));
const organismPanel = document.getElementById("organism-detail-panel");
const detailSectionsArchive = (() => {
  const node = document.getElementById("detail-sections-data");
  return node ? node.textContent.trim() : "";
})();
const sampleSortHeaders = Array.from(document.querySelectorAll(".sample-sort"));
const sampleTableBody = document.querySelector(".panel .table-wrap table tbody");
let sampleSortState = { key: "sample", direction: "asc" };
let detailZipStatePromise = null;
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
  const alignmentInsertIndex = options.blastAlignment
    ? Math.max(
        visibleIndices.findIndex((item) => item.key === "sscinames"),
        -1,
      ) + 1
    : -1;

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
  const binary = atob(base64Text);
  const bytes = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index += 1) {
    bytes[index] = binary.charCodeAt(index);
  }
  return bytes;
}

function parseZipStoredEntries(bytes) {
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let eocdOffset = -1;
  for (let offset = bytes.length - 22; offset >= 0; offset -= 1) {
    if (view.getUint32(offset, true) === 0x06054b50) {
      eocdOffset = offset;
      break;
    }
  }
  if (eocdOffset < 0) {
    throw new Error("ZIP end-of-central-directory not found.");
  }

  const totalEntries = view.getUint16(eocdOffset + 10, true);
  const centralDirectoryOffset = view.getUint32(eocdOffset + 16, true);
  const textDecoder = new TextDecoder();
  const entries = new Map();
  let cursor = centralDirectoryOffset;

  for (let index = 0; index < totalEntries; index += 1) {
    if (view.getUint32(cursor, true) !== 0x02014b50) {
      throw new Error("Invalid ZIP central directory entry.");
    }
    const compressionMethod = view.getUint16(cursor + 10, true);
    const compressedSize = view.getUint32(cursor + 20, true);
    const filenameLength = view.getUint16(cursor + 28, true);
    const extraLength = view.getUint16(cursor + 30, true);
    const commentLength = view.getUint16(cursor + 32, true);
    const localHeaderOffset = view.getUint32(cursor + 42, true);
    const filenameBytes = bytes.slice(cursor + 46, cursor + 46 + filenameLength);
    const filename = textDecoder.decode(filenameBytes);
    entries.set(filename, {
      compressionMethod,
      compressedSize,
      localHeaderOffset,
    });
    cursor += 46 + filenameLength + extraLength + commentLength;
  }

  return { bytes, entries };
}

async function loadDetailZipState() {
  if (!detailSectionsArchive) {
    return { bytes: new Uint8Array(), entries: new Map() };
  }
  if (!detailZipStatePromise) {
    detailZipStatePromise = Promise.resolve(parseZipStoredEntries(decodeBase64ToBytes(detailSectionsArchive)));
  }
  return detailZipStatePromise;
}

function readStoredZipEntry(archive, filename) {
  const entry = archive.entries.get(filename);
  if (!entry || entry.compressionMethod !== 0) {
    return null;
  }

  const view = new DataView(archive.bytes.buffer, archive.bytes.byteOffset, archive.bytes.byteLength);
  const offset = entry.localHeaderOffset;
  if (view.getUint32(offset, true) !== 0x04034b50) {
    throw new Error("Invalid ZIP local file header.");
  }
  const filenameLength = view.getUint16(offset + 26, true);
  const extraLength = view.getUint16(offset + 28, true);
  const dataStart = offset + 30 + filenameLength + extraLength;
  return archive.bytes.slice(dataStart, dataStart + entry.compressedSize);
}

async function gunzipBytes(bytes) {
  if (!("DecompressionStream" in window)) {
    throw new Error("This browser does not support gzip detail payloads.");
  }
  const response = new Response(new Response(bytes).body.pipeThrough(new DecompressionStream("gzip")));
  return response.text();
}

async function loadDetailData(panelId) {
  if (detailCache.has(panelId)) {
    return detailCache.get(panelId);
  }
  const archive = await loadDetailZipState();
  const entryBytes = readStoredZipEntry(archive, `details/${panelId}.json.gz`);
  if (!entryBytes) {
    return null;
  }
  const detail = JSON.parse(await gunzipBytes(entryBytes));
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

function populateFilterOptions(select) {
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

function applyTaxpastaFilter(view) {
  if (!view) {
    return;
  }

  const filter = view.querySelector("[data-superkingdom-filter]");
  const flagFilter = view.querySelector("[data-flag-filter]");
  const search = view.querySelector("[data-organism-search]");
  const metavalOnlyButton = view.querySelector("[data-metaval-only]");
  const counter = view.querySelector("[data-organisms-shown]");
  const taxonomyRows = Array.from(view.querySelectorAll(".taxonomy-row"));
  let shown = 0;
  const searchValue = search ? search.value.trim().toLowerCase() : "";
  const metavalOnly = metavalOnlyButton ? metavalOnlyButton.dataset.active === "true" : false;

  taxonomyRows.forEach((row) => {
    const matchesSuperkingdom = !filter || filter.value === "all" || row.dataset.superkingdom === filter.value;
    const matchesFlag = !flagFilter || flagFilter.value === "all" || row.dataset.flag === flagFilter.value;
    const matchesSearch = !searchValue || row.textContent.toLowerCase().includes(searchValue);
    const matchesMetaval = !metavalOnly || row.dataset.metavalChecked === "true";
    const matches = matchesSuperkingdom && matchesFlag && matchesSearch && matchesMetaval;
    row.hidden = !matches;

    const lineageRow = row.nextElementSibling;
    if (lineageRow && lineageRow.classList.contains("taxpasta-lineage-row") && !matches) {
      lineageRow.hidden = true;
      const toggle = row.querySelector(".lineage-toggle");
      if (toggle) {
        toggle.textContent = "Show";
      }
    }

    if (matches) {
      shown += 1;
    }
  });

  if (counter) {
    counter.textContent = String(shown);
  }
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
  applyTaxpastaFilter(panel.querySelector(`[data-classifier-view="${classifier}"]`));
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
  populateFilterOptions(select);
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
      applyTaxpastaFilter(panel.querySelector(`[data-classifier-view="${button.dataset.classifierTab}"]`));
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
    search.addEventListener("input", () => applyTaxpastaFilter(view));
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
  applyTaxpastaFilter(view);
});

document.querySelectorAll(".lineage-toggle").forEach((button) => {
  button.addEventListener("click", (event) => {
    event.stopPropagation();
    const target = document.getElementById(button.dataset.lineageTarget);
    if (!target) {
      return;
    }
    const isHidden = target.hidden;
    target.hidden = !isHidden;
    button.textContent = isHidden ? "Hide" : "Show";
    button.classList.toggle("is-open", isHidden);
  });
});

document.querySelectorAll(".taxonomy-row").forEach((row) => {
  if (!row.dataset.metavalPanel) {
    return;
  }
  row.addEventListener("click", () => {
    openOrganismPanel(row.dataset.metavalPanel);
  });
});

restorePanelFromUrl();
updateSampleSortHeaders();
