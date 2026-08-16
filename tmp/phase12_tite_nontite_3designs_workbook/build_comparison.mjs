import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const ROOT = process.cwd();
const OUTPUT_DIR = path.join(ROOT, "Presentation 8-17-2026", "Table and Plots");
const STEM = "phase12_tite_vs_nontite_3designs_N30_fut0p85_scenarios_1_16_20_24_27_38";
const CSV_PATH = path.join(OUTPUT_DIR, `${STEM}_comparison.csv`);
const XLSX_PATH = path.join(OUTPUT_DIR, `${STEM}_comparison.xlsx`);
const PREVIEW_PATH = path.join(ROOT, "tmp", "phase12_tite_nontite_3designs_workbook", "comparison_preview.png");
const SCENARIOS = [1, 16, 20, 24, 27, 38];
const DESIGNS = [
  { key: "one_stage", short: "one-stage" },
  { key: "two_stage_highest_utility", short: "2-stage highest utility" },
  { key: "two_stage_top2_randomized", short: "2-stage top-2 randomized" },
];
const NAVY = "#173F5F";
const TEAL = "#0F766E";
const GREY = "#E7E6E6";
const PALE_BLUE = "#DDEBF7";
const PALE_GREEN = "#E2F0D9";
const PALE_RED = "#FCE4D6";
const TEXT = "#1F2937";
const GRID = "#CBD5E1";

function parseCsvLine(line) {
  const values = [];
  let value = "";
  let quoted = false;
  for (let i = 0; i < line.length; i += 1) {
    const character = line[i];
    if (character === '"') {
      if (quoted && line[i + 1] === '"') {
        value += '"';
        i += 1;
      } else {
        quoted = !quoted;
      }
    } else if (character === "," && !quoted) {
      values.push(value);
      value = "";
    } else {
      value += character;
    }
  }
  values.push(value);
  return values;
}

function colName(number) {
  let n = number;
  let name = "";
  while (n > 0) {
    const remainder = (n - 1) % 26;
    name = String.fromCharCode(65 + remainder) + name;
    n = Math.floor((n - 1) / 26);
  }
  return name;
}

function style(sheet, address, format) {
  sheet.getRange(address).format = format;
}

async function main() {
  const csvText = await fs.readFile(CSV_PATH, "utf8");
  const lines = csvText.trim().split(/\r?\n/);
  const headers = parseCsvLine(lines[0]);
  const records = lines.slice(1).map((line, index) => {
    const values = parseCsvLine(line);
    return { excelRow: index + 2, record: Object.fromEntries(headers.map((header, i) => [header, values[i] ?? ""])) };
  });
  const rowMap = new Map(records.map(({ excelRow, record }) => [`${record.Version}|${record.Design_Key}|${record.Scenario}`, excelRow]));
  const sourceRow = (version, designKey, scenario) => {
    const row = rowMap.get(`${version}|${designKey}|${scenario}`);
    if (!row) throw new Error(`Missing source row: ${version}, ${designKey}, scenario ${scenario}.`);
    return row;
  };
  const sourceRecord = (version, designKey, scenario) => records.find(({ record }) => record.Version === version && record.Design_Key === designKey && Number(record.Scenario) === scenario).record;
  const headerColumn = (header) => {
    const index = headers.indexOf(header);
    if (index < 0) throw new Error(`Missing column: ${header}`);
    return colName(index + 1);
  };
  const sourceFormula = (row, header) => `='Source Data'!$${headerColumn(header)}$${row}`;

  const workbook = await Workbook.fromCSV(csvText, { sheetName: "Source Data" });
  const source = workbook.worksheets.getItem("Source Data");
  const sourceLastColumn = colName(headers.length);
  const sourceRange = `A1:${sourceLastColumn}${records.length + 1}`;
  source.showGridLines = false;
  style(source, sourceRange, { font: { color: TEXT, size: 10 }, verticalAlignment: "center", borders: { preset: "all", style: "thin", color: GRID } });
  style(source, `A1:${sourceLastColumn}1`, { fill: NAVY, font: { color: "#FFFFFF", bold: true, size: 10 }, horizontalAlignment: "center", verticalAlignment: "center", wrapText: true, borders: { preset: "all", style: "thin", color: NAVY } });
  source.getRange(`A1:${sourceLastColumn}${records.length + 1}`).format.autofitColumns();
  source.getRange("A:A").format.columnWidth = 12;
  source.getRange("B:B").format.columnWidth = 14;
  source.getRange("C:C").format.columnWidth = 30;
  source.getRange("D:D").format.columnWidth = 30;
  source.getRange("E:E").format.columnWidth = 36;
  source.getRange("F:F").format.columnWidth = 11;
  source.getRange("G:G").format.columnWidth = 55;
  source.getRange(`H:${sourceLastColumn}`).format.columnWidth = 15;
  source.getRange(`A1:${sourceLastColumn}1`).format.rowHeight = 32;
  source.freezePanes.freezeRows(1);
  source.tables.add(sourceRange, true, "ThreeDesignSource");

  const comparison = workbook.worksheets.add("Comparison");
  comparison.showGridLines = false;
  comparison.getRange("A1:G1").merge();
  comparison.getRange("A1").values = [["Phase I/II Operating Characteristics"]];
  style(comparison, "A1:G1", { fill: NAVY, font: { color: "#FFFFFF", bold: true, size: 16 }, horizontalAlignment: "center", verticalAlignment: "center" });
  comparison.getRange("A1:G1").format.rowHeight = 30;
  comparison.getRange("A2:G2").merge();
  comparison.getRange("A2").values = [["TITE vs new non-TITE | one-stage, two-stage highest utility, and two-stage top-2 randomized"]];
  style(comparison, "A2:G2", { fill: "#EDF2F7", font: { color: NAVY, bold: true, size: 10 }, horizontalAlignment: "center", verticalAlignment: "center" });
  comparison.getRange("A4:G4").merge();
  comparison.getRange("A4").values = [["Comparison scope: matched N = 30 configurations with 1,000 simulations per populated row; three allocation designs are shown for both versions."]];
  comparison.getRange("A5:G5").merge();
  comparison.getRange("A5").values = [["Sources: Presentation 8-17 TITE IDX 1001-2000 and new non-TITE IDX 0001-1000; all three designs and all six scenarios are available."]];
  style(comparison, "A4:G5", { fill: "#FFF7ED", font: { color: "#7C2D12", italic: true, size: 9 }, wrapText: true, verticalAlignment: "center" });
  comparison.getRange("A4:G5").format.rowHeight = 22;
  comparison.getRange("A7:G7").values = [["Metric", "Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "No OBD %"]];
  style(comparison, "A7:G7", { fill: TEAL, font: { color: "#FFFFFF", bold: true, size: 10 }, horizontalAlignment: "center", verticalAlignment: "center", borders: { preset: "all", style: "thin", color: TEAL } });
  comparison.getRange("A7:G7").format.rowHeight = 22;

  function section(row, label, labels = null) {
    comparison.getRange(`A${row}`).values = [[label]];
    if (labels) comparison.getRange(`B${row}:D${row}`).values = [[...labels]];
    style(comparison, `A${row}:G${row}`, { fill: GREY, font: { color: TEXT, bold: true, size: 10 }, borders: { preset: "all", style: "thin", color: GRID }, verticalAlignment: "center" });
    style(comparison, `B${row}:G${row}`, { horizontalAlignment: "center" });
    comparison.getRange(`A${row}:G${row}`).format.rowHeight = 18;
  }

  function dataRow(row, label, formulas, fill, bold = false, unavailable = false, noObd = true) {
    comparison.getRange(`A${row}`).values = [[label]];
    if (unavailable) {
      comparison.getRange(`B${row}:G${row}`).values = [["N/A", "N/A", "N/A", "N/A", "N/A", noObd ? "N/A" : ""]];
    } else {
      comparison.getRange(`B${row}:G${row}`).formulas = [formulas];
    }
    style(comparison, `A${row}:G${row}`, { fill, font: { color: TEXT, bold, size: 10 }, borders: { preset: "all", style: "thin", color: GRID }, verticalAlignment: "center" });
    style(comparison, `B${row}:G${row}`, { horizontalAlignment: "center" });
    comparison.getRange(`A${row}:G${row}`).format.rowHeight = 18;
  }

  for (let scenarioIndex = 0; scenarioIndex < SCENARIOS.length; scenarioIndex += 1) {
    const scenario = SCENARIOS[scenarioIndex];
    const start = 8 + (scenarioIndex * 26);
    const truthRow = sourceRow("New non-TITE", "two_stage_top2_randomized", scenario);
    comparison.getRange(`A${start}:G${start}`).merge();
    comparison.getRange(`A${start}`).values = [[`Scenario ${scenario}`]];
    style(comparison, `A${start}:G${start}`, { fill: "#D9E2F3", font: { color: NAVY, bold: true, size: 11 }, borders: { preset: "all", style: "thin", color: GRID }, verticalAlignment: "center" });
    comparison.getRange(`A${start}:G${start}`).format.rowHeight = 20;
    comparison.getRange(`A${start + 1}`).values = [["DLT rate"]];
    comparison.getRange(`B${start + 1}:F${start + 1}`).formulas = [[1, 2, 3, 4, 5].map((dose) => sourceFormula(truthRow, `True_DLT_Dose_${dose}`))];
    comparison.getRange(`A${start + 2}`).values = [["Efficacy rate"]];
    comparison.getRange(`B${start + 2}:F${start + 2}`).formulas = [[1, 2, 3, 4, 5].map((dose) => sourceFormula(truthRow, `True_Efficacy_Dose_${dose}`))];
    style(comparison, `A${start + 1}:G${start + 2}`, { font: { color: TEXT, size: 10 }, borders: { preset: "all", style: "thin", color: GRID }, verticalAlignment: "center" });
    style(comparison, `B${start + 1}:G${start + 2}`, { horizontalAlignment: "center" });
    comparison.getRange(`B${start + 1}:F${start + 2}`).format.numberFormat = "0.00";
    comparison.getRange(`A${start + 1}:G${start + 2}`).format.rowHeight = 18;

    section(start + 3, "OBD selection (%)");
    let row = start + 4;
    for (const design of DESIGNS) {
      for (const version of ["TITE", "New non-TITE"]) {
        const sourceDataRow = sourceRow(version, design.key, scenario);
        const record = sourceRecord(version, design.key, scenario);
        const isAvailable = record.Available === "True";
        const isTite = version === "TITE";
        const formulas = [...[1, 2, 3, 4, 5].map((dose) => sourceFormula(sourceDataRow, `OBD_Selection_Dose_${dose}_pct`)), sourceFormula(sourceDataRow, "No_OBD_Selection_pct")];
        dataRow(row, `${version} | ${design.short}${isAvailable ? "" : " (rerun)"}`, formulas, isAvailable ? (isTite ? PALE_BLUE : PALE_GREEN) : PALE_RED, !isTite, !isAvailable);
        comparison.getRange(`B${row}:G${row}`).format.numberFormat = "0.0";
        row += 1;
      }
    }
    section(row, "Mean administrations by dose");
    row += 1;
    for (const design of DESIGNS) {
      for (const version of ["TITE", "New non-TITE"]) {
        const sourceDataRow = sourceRow(version, design.key, scenario);
        const record = sourceRecord(version, design.key, scenario);
        const isAvailable = record.Available === "True";
        const isTite = version === "TITE";
        const formulas = [...[1, 2, 3, 4, 5].map((dose) => sourceFormula(sourceDataRow, `Mean_Administrations_Dose_${dose}`)), ""];
        dataRow(row, `${version} | ${design.short}${isAvailable ? "" : " (rerun)"}`, formulas, isAvailable ? (isTite ? PALE_BLUE : PALE_GREEN) : PALE_RED, !isTite, !isAvailable, false);
        comparison.getRange(`B${row}:F${row}`).format.numberFormat = "0.0";
        row += 1;
      }
    }
    section(row, "Run-level means", ["Total admin.", "Unique patients", "Duration days"]);
    row += 1;
    for (const design of DESIGNS) {
      for (const version of ["TITE", "New non-TITE"]) {
        const sourceDataRow = sourceRow(version, design.key, scenario);
        const record = sourceRecord(version, design.key, scenario);
        const isAvailable = record.Available === "True";
        const isTite = version === "TITE";
        const formulas = [sourceFormula(sourceDataRow, "Mean_Total_Administrations"), sourceFormula(sourceDataRow, "Mean_Total_Unique_Patients"), sourceFormula(sourceDataRow, "Mean_Duration_days"), "", "", ""];
        dataRow(row, `${version} | ${design.short}${isAvailable ? "" : " (rerun)"}`, formulas, isAvailable ? (isTite ? PALE_BLUE : PALE_GREEN) : PALE_RED, !isTite, !isAvailable, false);
        comparison.getRange(`B${row}:D${row}`).format.numberFormat = "0.0";
        row += 1;
      }
    }
  }

  comparison.getRange("A:A").format.columnWidth = 39;
  comparison.getRange("B:F").format.columnWidth = 12;
  comparison.getRange("G:G").format.columnWidth = 13;
  comparison.freezePanes.freezeRows(7);
  const used = comparison.getUsedRange();
  const formulaErrors = used.values.flat().filter((value) => typeof value === "string" && /^#(REF!|DIV\/0!|VALUE!|NAME\?|N\/A)$/.test(value));
  if (formulaErrors.length) throw new Error(`Workbook contains formula errors: ${formulaErrors.join(", ")}`);
  const inspection = await workbook.inspect({ kind: "region,formula", sheetId: "Comparison", range: "A1:G164", maxChars: 8000, tableMaxRows: 10, tableMaxCols: 7 });
  console.log(inspection.ndjson ?? inspection);
  const preview = await workbook.render({ sheetName: "Comparison", autoCrop: "all", scale: 1, format: "png" });
  await fs.writeFile(PREVIEW_PATH, new Uint8Array(await preview.arrayBuffer()));
  const output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(XLSX_PATH);
  console.log(XLSX_PATH);
  console.log(PREVIEW_PATH);
}

await main();
