import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectRoot = process.cwd();
const presentationDir = path.join(projectRoot, "Presentation 8-17-2026");
const sourceNewDesign = path.join(
  presentationDir,
  "AIDE_phase_I_II_modelsadditive_newdesign_shared_N30_ncycle2_rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv",
);
const sourceBoin12 = path.join(
  presentationDir,
  "BOIN12_v1.4.2.0_Operating Characteristics_2026-08-09 170536.992978_fut0.85.csv",
);
const outputDir = path.join(projectRoot, "outputs", "alpha_comparison_nontite_newdesign_boin_20260816");
const outputPath = path.join(outputDir, "nontite_newdesign_alpha_and_boin_comparison.xlsx");

const scenarios = [1, 16, 20, 24, 27, 38];
const alphas = [0, 0.3, 0.6, 0.9];
const designOrder = ["one_stage", "two_stage_highest_utility", "two_stage_top2_randomized"];
const colors = {
  navy: "#173F5F",
  teal: "#0F766E",
  blue: "#DCEAF7",
  green: "#E2F0D9",
  slate: "#F1F5F9",
  gray: "#E5E7EB",
  ink: "#1F2937",
  note: "#FFF7ED",
};

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i += 1) {
    const c = text[i];
    if (c === '"') {
      if (quoted && text[i + 1] === '"') {
        field += '"';
        i += 1;
      } else {
        quoted = !quoted;
      }
    } else if (c === "," && !quoted) {
      row.push(field);
      field = "";
    } else if ((c === "\n" || c === "\r") && !quoted) {
      if (c === "\r" && text[i + 1] === "\n") i += 1;
      row.push(field);
      if (row.some((value) => value !== "")) rows.push(row);
      row = [];
      field = "";
    } else {
      field += c;
    }
  }
  if (field !== "" || row.length) {
    row.push(field);
    if (row.some((value) => value !== "")) rows.push(row);
  }
  return rows;
}

function csvObjects(text) {
  const [headers, ...data] = parseCsv(text);
  return data.map((values) => Object.fromEntries(headers.map((header, i) => [header, values[i] ?? ""])));
}

function num(value) {
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : null;
}

function designInfo(row) {
  if (row.Allocation === "one_stage") {
    return { key: "one_stage", label: "One-stage (highest utility)", stageRule: "Highest utility throughout" };
  }
  if (row.Stage2_Allocation === "highest_utility") {
    return { key: "two_stage_highest_utility", label: "Two-stage (highest utility)", stageRule: "Stage II: highest utility" };
  }
  if (row.Stage2_Allocation === "top2_randomized") {
    return { key: "two_stage_top2_randomized", label: "Two-stage (top-2 randomized)", stageRule: "Stage II: top-2 randomized" };
  }
  throw new Error(`Unknown design: ${row.Allocation} / ${row.Stage2_Allocation}`);
}

function alphaKey(value) {
  return Number(value).toFixed(1);
}

function colLetter(index) {
  let n = index;
  let output = "";
  while (n > 0) {
    const remainder = (n - 1) % 26;
    output = String.fromCharCode(65 + remainder) + output;
    n = Math.floor((n - 1) / 26);
  }
  return output;
}

function applyTitle(sheet, title, subtitle, note, lastColumn) {
  sheet.getRange(`A1:${lastColumn}1`).merge();
  sheet.getRange("A1").values = [[title]];
  sheet.getRange("A1").format = {
    fill: colors.navy,
    font: { bold: true, color: "#FFFFFF", size: 16 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
  };
  sheet.getRange(`A2:${lastColumn}2`).merge();
  sheet.getRange("A2").values = [[subtitle]];
  sheet.getRange("A2").format = {
    fill: "#EAF2F8",
    font: { bold: true, color: colors.navy, size: 10 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
  };
  sheet.getRange(`A4:${lastColumn}4`).merge();
  sheet.getRange("A4").values = [[note]];
  sheet.getRange("A4").format = {
    fill: colors.note,
    font: { italic: true, color: "#7C2D12", size: 9 },
    wrapText: true,
    verticalAlignment: "center",
  };
  sheet.getRange("A1").format.rowHeight = 26;
  sheet.getRange("A2").format.rowHeight = 20;
  sheet.getRange("A4").format.rowHeight = 30;
}

function styleTable(sheet, startRow, lastRow, lastColumn, tableName, widths = {}) {
  const header = sheet.getRange(`A${startRow}:${lastColumn}${startRow}`);
  header.format = {
    fill: colors.teal,
    font: { bold: true, color: "#FFFFFF", size: 9 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
    borders: { preset: "outside", style: "thin", color: "#0B5E57" },
  };
  header.format.rowHeight = 30;
  const body = sheet.getRange(`A${startRow + 1}:${lastColumn}${lastRow}`);
  body.format = {
    font: { color: colors.ink, size: 9 },
    verticalAlignment: "center",
    borders: { insideHorizontal: { style: "thin", color: "#E2E8F0" } },
  };
  const table = sheet.tables.add(`A${startRow}:${lastColumn}${lastRow}`, true, tableName);
  table.style = "TableStyleMedium2";
  table.showBandedColumns = false;
  for (const [column, width] of Object.entries(widths)) {
    sheet.getRange(`${column}:${column}`).format.columnWidth = width;
  }
}

function setAlphaHighlight(sheet, range) {
  sheet.getRange(range).conditionalFormats.add("cellIs", {
    operator: "equal",
    formula: 0,
    format: { fill: "#E2F0D9", font: { bold: true, color: "#166534" } },
  });
  sheet.getRange(range).conditionalFormats.add("cellIs", {
    operator: "equal",
    formula: 0.3,
    format: { fill: "#DBEAFE", font: { bold: true, color: "#1D4ED8" } },
  });
  sheet.getRange(range).conditionalFormats.add("cellIs", {
    operator: "equal",
    formula: 0.6,
    format: { fill: "#FEF3C7", font: { bold: true, color: "#92400E" } },
  });
  sheet.getRange(range).conditionalFormats.add("cellIs", {
    operator: "equal",
    formula: 0.9,
    format: { fill: "#F3E8FF", font: { bold: true, color: "#7E22CE" } },
  });
}

const newDesignRows = csvObjects(await fs.readFile(sourceNewDesign, "utf8"));
const alphaRows = newDesignRows.filter((row) => {
  const toxAlpha = num(row.Toxicity_IPDE_Alpha);
  const efficacyAlpha = num(row.Efficacy_IPDE_Alpha);
  return scenarios.includes(num(row.Scenario)) && alphas.includes(toxAlpha) && toxAlpha === efficacyAlpha;
});
if (alphaRows.length !== 360) throw new Error(`Expected 360 selected source rows; found ${alphaRows.length}.`);

const grouped = new Map();
for (const row of alphaRows) {
  const design = designInfo(row);
  const key = [row.Scenario, design.key, alphaKey(row.Toxicity_IPDE_Alpha)].join("|");
  if (!grouped.has(key)) grouped.set(key, []);
  grouped.get(key).push(row);
}
if (grouped.size !== 72 || [...grouped.values()].some((rows) => rows.length !== 5)) {
  throw new Error("The requested source file does not contain a complete 6-scenario × 3-design × 4-alpha grid.");
}

const settings = [];
for (const rows of grouped.values()) {
  const sorted = [...rows].sort((a, b) => num(a.Dose) - num(b.Dose));
  const row = sorted[0];
  const design = designInfo(row);
  settings.push({
    scenario: num(row.Scenario),
    designKey: design.key,
    design: design.label,
    stageRule: design.stageRule,
    alpha: num(row.Toxicity_IPDE_Alpha),
    trueObd: num(row.True_OBD),
    dlt: sorted.map((item) => num(item.True_DLT_rate)),
    efficacy: sorted.map((item) => num(item.True_Efficacy_rate)),
    obd: sorted.map((item) => num(item.OBD_Selection_pct)),
    treated: sorted.map((item) => num(item.Pts_Treated)),
    noObd: num(row.No_OBD_Selection_pct),
    earlyStop: num(row.Early_Stopping_pct),
    totalAdministrations: num(row.Total_Administrations),
    uniquePatients: num(row.Total_Unique_Patients),
    duration: num(row.Duration),
    ntrial: num(row.ntrial_from_files),
  });
}
settings.sort((a, b) => a.scenario - b.scenario || designOrder.indexOf(a.designKey) - designOrder.indexOf(b.designKey) || a.alpha - b.alpha);

const trueObdByScenario = new Map();
for (const setting of settings) {
  if (trueObdByScenario.has(setting.scenario) && trueObdByScenario.get(setting.scenario) !== setting.trueObd) {
    throw new Error(`True OBD is inconsistent within Scenario ${setting.scenario}.`);
  }
  trueObdByScenario.set(setting.scenario, setting.trueObd);
}

const uboin = {
  1: { treated: [9, 5.4, 4.7, 3.9, 3.7], selection: [15.6, 17.4, 20.6, 18.7, 13.7], outcome: 14, source: "U-BOIN Sce 1-10 fut 0.85.pdf (Scenario 1)" },
  16: { treated: [6.1, 6.2, 6.4, 6.4, 5.4], selection: [28.1, 26.8, 21.1, 16.3, 7.7], outcome: 0, source: "U-BOIN Sce 11-20 fut 0.85.pdf (local Scenario 6)" },
  20: { treated: [5.9, 10.8, 7, 4, 2.4], selection: [8.5, 47.5, 18.6, 12.4, 13], outcome: 0, source: "U-BOIN Sce 11-20 fut 0.85.pdf (local Scenario 10)" },
  24: { treated: [3.5, 4, 5.1, 9.8, 7.9], selection: [1.5, 3.2, 8.2, 60.1, 27], outcome: 0, source: "U-BOIN Sce 21-30 fut 0.85.pdf (local Scenario 4)" },
  27: { treated: [6.9, 8.1, 8.8, 4.3, 1.8], selection: [8.7, 20.9, 38.9, 16.5, 14.5], outcome: 0.5, source: "U-BOIN Sce 21-30 fut 0.85.pdf (local Scenario 7)" },
  38: { treated: [5.2, 5.5, 5.6, 6.3, 7.8], selection: [6.2, 8.9, 13, 13.3, 58.6], outcome: 0, source: "UBOIN Sce38.pdf" },
};

const boinRows = parseCsv(await fs.readFile(sourceBoin12, "utf8"));
const boin = {};
let currentScenario = null;
for (const row of boinRows) {
  const first = (row[0] ?? "").trim();
  if (/^Scenario\s+\d+$/.test(first)) {
    currentScenario = Number(first.replace("Scenario", "").trim());
    if (!boin[currentScenario]) boin[currentScenario] = {};
  } else if (currentScenario !== null && first === "No. Pts treated") {
    boin[currentScenario].treated = row.slice(1, 6).map(num);
  } else if (currentScenario !== null && first === "Select %") {
    boin[currentScenario].selection = row.slice(1, 6).map(num);
    boin[currentScenario].outcome = num(row[6]);
    boin[currentScenario].source = "BOIN12_v1.4.2.0_Operating Characteristics_2026-08-09 170536.992978_fut0.85.csv";
  }
}
boin[38] = {
  treated: [4.2, 4.7, 5.2, 6.1, 9.6],
  selection: [13.2, 13.8, 16.5, 13.6, 41.5],
  outcome: 1.4,
  source: "BOIN12 Sce38.pdf",
};
for (const scenario of scenarios) {
  if (!uboin[scenario] || !boin[scenario]?.selection || !boin[scenario]?.treated) {
    throw new Error(`Missing comparator data for Scenario ${scenario}.`);
  }
}

const workbook = Workbook.create();
const overview = workbook.worksheets.add("Read Me");
const selectionSheet = workbook.worksheets.add("Alpha Selection");
const allocationSheet = workbook.worksheets.add("Alpha Allocation");
const methodsSheet = workbook.worksheets.add("Methods alpha 0");
const rawNewSheet = workbook.worksheets.add("New-design source");
const rawComparatorSheet = workbook.worksheets.add("Comparator sources");
for (const sheet of [overview, selectionSheet, allocationSheet, methodsSheet, rawNewSheet, rawComparatorSheet]) sheet.showGridLines = false;

// Read Me
overview.getRange("A1:B1").merge();
overview.getRange("A1").values = [["Non-TITE New-Design Alpha and Comparator Tables"]];
overview.getRange("A1").format = { fill: colors.navy, font: { bold: true, color: "#FFFFFF", size: 16 }, horizontalAlignment: "center" };
overview.getRange("A1").format.rowHeight = 28;
overview.getRange("A3:B7").values = [
  ["Scope", "Scenarios 1, 16, 20, 24, 27, and 38; 1,000 simulations per new-design row."],
  ["Alpha comparison", "New design only. Toxicity and efficacy IPDE alpha are matched at 0, 0.3, 0.6, and 0.9."],
  ["New-design variants", "One-stage highest utility; two-stage highest utility; two-stage top-2 randomized."],
  ["Comparator set", "New-design alpha = 0 rows are shown beside U-BOIN and BOIN12."],
  ["Important definition", "The outcome column is source-specific: new design = No OBD selection; U-BOIN = Stop %; BOIN12 = No selection %. Do not treat these as identical estimands."],
];
overview.getRange("A3:A7").format = { fill: colors.teal, font: { bold: true, color: "#FFFFFF" }, verticalAlignment: "center" };
overview.getRange("B3:B7").format = { fill: "#F8FAFC", wrapText: true, verticalAlignment: "center" };
overview.getRange("A3:B7").format.borders = { preset: "outside", style: "thin", color: "#CBD5E1" };
overview.getRange("A3:B7").format.rowHeight = 36;
overview.getRange("A:A").format.columnWidth = 22;
overview.getRange("B:B").format.columnWidth = 105;
overview.getRange("A10:B14").values = [
  ["Sheet", "Purpose"],
  ["Alpha Selection", "OBD-selection and no-OBD comparison across alpha values."],
  ["Alpha Allocation", "Mean dose allocation, total administrations, participants, duration, and early stopping across alpha values."],
  ["Methods alpha 0", "New design (all three variants) vs U-BOIN vs BOIN12 at new-design alpha = 0."],
  ["Source sheets", "Filtered long-form new-design source and the extracted comparator values."],
];
overview.getRange("A10:B10").format = { fill: colors.teal, font: { bold: true, color: "#FFFFFF" } };
overview.getRange("A10:B14").format.borders = { insideHorizontal: { style: "thin", color: "#E2E8F0" }, outside: { style: "thin", color: "#CBD5E1" } };

// Alpha Selection
const selectionHeaders = ["Scenario", "New-design allocation", "IPDE alpha", "True OBD", "Correct OBD selection (%)", "No OBD selection (%)", "OBD D1 (%)", "OBD D2 (%)", "OBD D3 (%)", "OBD D4 (%)", "OBD D5 (%)", "Delta correct OBD vs alpha 0 (pp)", "Delta no OBD vs alpha 0 (pp)"];
applyTitle(selectionSheet, "New Design: Alpha Sensitivity - OBD Selection", "Non-TITE results from the specified new-design summary; toxicity IPDE alpha = efficacy IPDE alpha.", "Rows are filterable. Delta columns are formula-driven and compare each scenario/design with its alpha = 0 row.", colLetter(selectionHeaders.length));
selectionSheet.getRange("A7:M7").values = [selectionHeaders];
const selectionValues = settings.map((setting) => [setting.scenario, setting.design, setting.alpha, setting.trueObd, null, setting.noObd, ...setting.obd, null, null]);
selectionSheet.getRange(`A8:M${7 + selectionValues.length}`).values = selectionValues;
for (let row = 8; row < 8 + settings.length; row += 1) {
  selectionSheet.getRange(`E${row}`).formulas = [[`=INDEX($G${row}:$K${row},1,$D${row})`]];
  selectionSheet.getRange(`L${row}`).formulas = [[`=IF($C${row}=0,0,$E${row}-SUMIFS($E$8:$E$79,$A$8:$A$79,$A${row},$B$8:$B$79,$B${row},$C$8:$C$79,0))`]];
  selectionSheet.getRange(`M${row}`).formulas = [[`=IF($C${row}=0,0,$F${row}-SUMIFS($F$8:$F$79,$A$8:$A$79,$A${row},$B$8:$B$79,$B${row},$C$8:$C$79,0))`]];
}
styleTable(selectionSheet, 7, 79, "M", "AlphaSelectionTable", { A: 10, B: 31, C: 12, D: 10, E: 19, F: 18, L: 19, M: 19 });
selectionSheet.getRange("C8:C79").format.numberFormat = "0.0";
selectionSheet.getRange("D8:D79").format.numberFormat = "0";
selectionSheet.getRange("E8:M79").format.numberFormat = "0.0";
setAlphaHighlight(selectionSheet, "C8:C79");
selectionSheet.freezePanes.freezeRows(7);
selectionSheet.freezePanes.freezeColumns(2);

// Alpha Allocation
const allocationHeaders = ["Scenario", "New-design allocation", "IPDE alpha", "Mean patients D1", "Mean patients D2", "Mean patients D3", "Mean patients D4", "Mean patients D5", "Total administrations", "Unique patients", "Duration (days)", "Early stopping (%)", "Delta total administrations vs alpha 0", "Delta duration vs alpha 0 (days)"];
applyTitle(allocationSheet, "New Design: Alpha Sensitivity - Allocation and Trial Burden", "Non-TITE results from the specified new-design summary; toxicity IPDE alpha = efficacy IPDE alpha.", "Mean patients treated by dose and run-level quantities. Delta columns are formula-driven versus alpha = 0 within each scenario/design.", colLetter(allocationHeaders.length));
allocationSheet.getRange("A7:N7").values = [allocationHeaders];
const allocationValues = settings.map((setting) => [setting.scenario, setting.design, setting.alpha, ...setting.treated, setting.totalAdministrations, setting.uniquePatients, setting.duration, setting.earlyStop, null, null]);
allocationSheet.getRange(`A8:N${7 + allocationValues.length}`).values = allocationValues;
for (let row = 8; row < 8 + settings.length; row += 1) {
  allocationSheet.getRange(`M${row}`).formulas = [[`=IF($C${row}=0,0,$I${row}-SUMIFS($I$8:$I$79,$A$8:$A$79,$A${row},$B$8:$B$79,$B${row},$C$8:$C$79,0))`]];
  allocationSheet.getRange(`N${row}`).formulas = [[`=IF($C${row}=0,0,$K${row}-SUMIFS($K$8:$K$79,$A$8:$A$79,$A${row},$B$8:$B$79,$B${row},$C$8:$C$79,0))`]];
}
styleTable(allocationSheet, 7, 79, "N", "AlphaAllocationTable", { A: 10, B: 31, C: 12, I: 18, J: 15, K: 15, L: 15, M: 23, N: 22 });
allocationSheet.getRange("C8:C79").format.numberFormat = "0.0";
allocationSheet.getRange("D8:N79").format.numberFormat = "0.0";
setAlphaHighlight(allocationSheet, "C8:C79");
allocationSheet.freezePanes.freezeRows(7);
allocationSheet.freezePanes.freezeColumns(2);

// Comparator table
const methodHeaders = ["Scenario", "Method", "New-design stage II rule", "IPDE alpha", "True OBD", "Correct OBD selection (%)", "No OBD / stop (%)", "Outcome definition", "Select D1 (%)", "Select D2 (%)", "Select D3 (%)", "Select D4 (%)", "Select D5 (%)", "Mean patients D1", "Mean patients D2", "Mean patients D3", "Mean patients D4", "Mean patients D5", "Source"];
applyTitle(methodsSheet, "Alpha = 0: New Design vs U-BOIN vs BOIN12", "All three non-TITE new-design variants at alpha = 0 are shown alongside U-BOIN and BOIN12 for the same six scenarios.", "Comparator outcomes are reported exactly as their sources label them. U-BOIN and BOIN12 do not use/report the new design's IPDE alpha parameter.", colLetter(methodHeaders.length));
methodsSheet.getRange("A7:S7").values = [methodHeaders];
const alphaZeroSettings = settings.filter((setting) => setting.alpha === 0);
const methodRows = [];
for (const scenario of scenarios) {
  for (const setting of alphaZeroSettings.filter((item) => item.scenario === scenario)) {
    methodRows.push({
      scenario,
      method: `New Design - ${setting.design}`,
      stage: setting.stageRule,
      alpha: 0,
      trueObd: setting.trueObd,
      outcome: setting.noObd,
      outcomeDefinition: "No OBD selection %",
      selection: setting.obd,
      treated: setting.treated,
      source: "New-design summary",
    });
  }
  methodRows.push({
    scenario,
    method: "U-BOIN",
    stage: "Not applicable",
    alpha: null,
    trueObd: trueObdByScenario.get(scenario),
    outcome: uboin[scenario].outcome,
    outcomeDefinition: "Stop %",
    selection: uboin[scenario].selection,
    treated: uboin[scenario].treated,
    source: "U-BOIN - see Comparator sources",
  });
  methodRows.push({
    scenario,
    method: "BOIN12",
    stage: "Not applicable",
    alpha: null,
    trueObd: trueObdByScenario.get(scenario),
    outcome: boin[scenario].outcome,
    outcomeDefinition: "No selection %",
    selection: boin[scenario].selection,
    treated: boin[scenario].treated,
    source: "BOIN12 - see Comparator sources",
  });
}
const methodValues = methodRows.map((item) => [item.scenario, item.method, item.stage, item.alpha, item.trueObd, null, item.outcome, item.outcomeDefinition, ...item.selection, ...item.treated, item.source]);
methodsSheet.getRange(`A8:S${7 + methodValues.length}`).values = methodValues;
for (let row = 8; row < 8 + methodRows.length; row += 1) {
  methodsSheet.getRange(`F${row}`).formulas = [[`=INDEX($I${row}:$M${row},1,$E${row})`]];
}
styleTable(methodsSheet, 7, 37, "S", "MethodsAlphaZeroTable", { A: 10, B: 35, C: 28, D: 12, E: 10, F: 19, G: 16, H: 22, N: 16, O: 16, P: 16, Q: 16, R: 16, S: 30 });
methodsSheet.getRange("D8:D37").format.numberFormat = "0.0";
methodsSheet.getRange("E8:E37").format.numberFormat = "0";
methodsSheet.getRange("F8:R37").format.numberFormat = "0.0";
methodsSheet.getRange("D8:D25").format = { fill: "#E2F0D9", font: { bold: true, color: "#166534" } };
methodsSheet.freezePanes.freezeRows(7);
methodsSheet.freezePanes.freezeColumns(2);

// Filtered raw new-design source
const rawHeaders = ["Task ID", "Scenario", "New-design allocation", "Stage II allocation", "IPDE alpha", "True OBD", "Dose", "True DLT rate", "True efficacy rate", "OBD selection (%)", "Mean patients treated", "No OBD selection (%)", "Total administrations", "Unique patients", "Duration (days)", "Early stopping (%)", "ntrial"];
applyTitle(rawNewSheet, "Filtered New-Design Source Data", "Long-form records retained from the user-specified non-TITE summary file.", "Six scenarios x three new-design allocation variants x four matched alpha values x five dose rows = 360 source records.", colLetter(rawHeaders.length));
rawNewSheet.getRange("A7:Q7").values = [rawHeaders];
const rawRows = alphaRows
  .map((row) => {
    const design = designInfo(row);
    return [num(row.Task_ID), num(row.Scenario), design.label, row.Stage2_Allocation, num(row.Toxicity_IPDE_Alpha), num(row.True_OBD), num(row.Dose), num(row.True_DLT_rate), num(row.True_Efficacy_rate), num(row.OBD_Selection_pct), num(row.Pts_Treated), num(row.No_OBD_Selection_pct), num(row.Total_Administrations), num(row.Total_Unique_Patients), num(row.Duration), num(row.Early_Stopping_pct), num(row.ntrial_from_files)];
  })
  .sort((a, b) => a[1] - b[1] || designOrder.indexOf(designInfo({ Allocation: a[2].startsWith("One") ? "one_stage" : "two_stage", Stage2_Allocation: a[3] }).key) - designOrder.indexOf(designInfo({ Allocation: b[2].startsWith("One") ? "one_stage" : "two_stage", Stage2_Allocation: b[3] }).key) || a[4] - b[4] || a[6] - b[6]);
rawNewSheet.getRange(`A8:Q${7 + rawRows.length}`).values = rawRows;
styleTable(rawNewSheet, 7, 367, "Q", "NewDesignSourceTable", { A: 10, B: 10, C: 31, D: 20, E: 12, F: 10, G: 8, J: 18, K: 20, L: 19, M: 18, N: 15, O: 15, P: 16, Q: 10 });
rawNewSheet.getRange("E8:E367").format.numberFormat = "0.0";
rawNewSheet.getRange("F8:G367").format.numberFormat = "0";
rawNewSheet.getRange("H8:P367").format.numberFormat = "0.0";
rawNewSheet.getRange("Q8:Q367").format.numberFormat = "#,##0";
setAlphaHighlight(rawNewSheet, "E8:E367");
rawNewSheet.freezePanes.freezeRows(7);
rawNewSheet.freezePanes.freezeColumns(2);

// Comparator source transcription with auditable source labels.
const comparatorHeaders = ["Scenario", "Method", "Source", "Outcome definition", "Outcome (%)", "Select D1 (%)", "Select D2 (%)", "Select D3 (%)", "Select D4 (%)", "Select D5 (%)", "Mean patients D1", "Mean patients D2", "Mean patients D3", "Mean patients D4", "Mean patients D5"];
applyTitle(rawComparatorSheet, "Comparator Source Values", "U-BOIN and BOIN12 values transcribed from the Presentation 8-17-2026 source PDFs/CSV.", "U-BOIN ranges 11-20 and 21-30 number scenarios locally; the source references identify the local scenario where applicable.", colLetter(comparatorHeaders.length));
rawComparatorSheet.getRange("A7:O7").values = [comparatorHeaders];
const comparatorValues = [];
for (const scenario of scenarios) {
  comparatorValues.push([scenario, "U-BOIN", uboin[scenario].source, "Stop %", uboin[scenario].outcome, ...uboin[scenario].selection, ...uboin[scenario].treated]);
  comparatorValues.push([scenario, "BOIN12", boin[scenario].source, "No selection %", boin[scenario].outcome, ...boin[scenario].selection, ...boin[scenario].treated]);
}
rawComparatorSheet.getRange(`A8:O${7 + comparatorValues.length}`).values = comparatorValues;
styleTable(rawComparatorSheet, 7, 19, "O", "ComparatorSourceTable", { A: 10, B: 14, C: 58, D: 22, E: 14, K: 16, L: 16, M: 16, N: 16, O: 16 });
rawComparatorSheet.getRange("E8:O19").format.numberFormat = "0.0";
rawComparatorSheet.freezePanes.freezeRows(7);
rawComparatorSheet.freezePanes.freezeColumns(2);

for (const sheet of [selectionSheet, allocationSheet, methodsSheet, rawNewSheet, rawComparatorSheet]) {
  sheet.getUsedRange().format.autofitRows();
}

const selectionPreview = await workbook.render({ sheetName: "Alpha Selection", range: "A1:M34", scale: 1.2, format: "png" });
await fs.mkdir(outputDir, { recursive: true });
await fs.writeFile(path.join(outputDir, "alpha_selection_preview.png"), new Uint8Array(await selectionPreview.arrayBuffer()));
const methodsPreview = await workbook.render({ sheetName: "Methods alpha 0", range: "A1:S37", scale: 0.8, format: "png" });
await fs.writeFile(path.join(outputDir, "methods_alpha0_preview.png"), new Uint8Array(await methodsPreview.arrayBuffer()));
const sourcePreview = await workbook.render({ sheetName: "New-design source", range: "A1:Q24", scale: 0.8, format: "png" });
await fs.writeFile(path.join(outputDir, "source_preview.png"), new Uint8Array(await sourcePreview.arrayBuffer()));

const inspect = await workbook.inspect({
  kind: "table",
  range: "Alpha Selection!A7:M16",
  include: "values,formulas",
  tableMaxRows: 10,
  tableMaxCols: 13,
});
console.log(inspect.ndjson);
const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 200 },
  summary: "final formula error scan",
});
console.log(errors.ndjson);

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);
console.log(`OUTPUT=${outputPath}`);
