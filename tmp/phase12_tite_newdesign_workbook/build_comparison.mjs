import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const ROOT = process.cwd();
const OUTPUT_DIR = path.join(ROOT, "Presentation 8-17-2026", "Table and Plots");
const STEM = "phase12_tite_vs_newdesign_N30_fut0p85_scenarios_1_16_20_24_27_38";
const CSV_PATH = path.join(OUTPUT_DIR, `${STEM}_comparison.csv`);
const XLSX_PATH = path.join(OUTPUT_DIR, `${STEM}_comparison.xlsx`);
const PREVIEW_PATH = path.join(ROOT, "tmp", "phase12_tite_newdesign_workbook", "comparison_preview.png");

const SCENARIOS = [1, 16, 20, 24, 27, 38];
const TITE = "TITE design";
const NEW = "New non-TITE design";
const NAVY = "#173F5F";
const TEAL = "#0F766E";
const GREY = "#E7E6E6";
const PALE_BLUE = "#DDEBF7";
const PALE_GREEN = "#E2F0D9";
const PALE_RED = "#FCE4D6";
const TEXT = "#1F2937";
const GRID = "#CBD5E1";

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

function sourceFormula(rowNumber, columnName) {
  return `='Source Data'!$${columnName}$${rowNumber}`;
}

function styleRange(sheet, address, format) {
  sheet.getRange(address).format = format;
}

async function main() {
  const csvText = await fs.readFile(CSV_PATH, "utf8");
  const lines = csvText.trim().split(/\r?\n/);
  const headers = lines[0].split(",");
  const records = lines.slice(1).map((line, index) => {
    const values = line.split(",");
    return { excelRow: index + 2, values, record: Object.fromEntries(headers.map((header, i) => [header, values[i]])) };
  });

  const sourceRows = new Map(records.map(({ excelRow, record }) => [`${record.Method}|${record.Scenario}`, excelRow]));
  const sourceRow = (method, scenario) => {
    const row = sourceRows.get(`${method}|${scenario}`);
    if (!row) throw new Error(`Missing source row for ${method}, scenario ${scenario}.`);
    return row;
  };
  const headerCol = (header) => {
    const index = headers.indexOf(header);
    if (index < 0) throw new Error(`Missing expected source column: ${header}`);
    return colName(index + 1);
  };
  const field = (row, header) => sourceFormula(row, headerCol(header));

  const workbook = await Workbook.fromCSV(csvText, { sheetName: "Source Data" });
  const source = workbook.worksheets.getItem("Source Data");
  const sourceLastColumn = colName(headers.length);
  const sourceRange = `A1:${sourceLastColumn}${records.length + 1}`;
  source.showGridLines = false;
  styleRange(source, sourceRange, {
    font: { color: TEXT, size: 10 },
    verticalAlignment: "center",
    borders: { preset: "all", style: "thin", color: GRID },
  });
  styleRange(source, `A1:${sourceLastColumn}1`, {
    fill: NAVY,
    font: { color: "#FFFFFF", bold: true, size: 10 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
    borders: { preset: "all", style: "thin", color: NAVY },
  });
  source.getRange(`A1:${sourceLastColumn}${records.length + 1}`).format.autofitColumns();
  source.getRange("A:A").format.columnWidth = 12;
  source.getRange("B:B").format.columnWidth = 30;
  source.getRange("C:C").format.columnWidth = 11;
  source.getRange("D:D").format.columnWidth = 57;
  source.getRange(`E:${sourceLastColumn}`).format.columnWidth = 15;
  source.getRange("A1:" + sourceLastColumn + "1").format.rowHeight = 32;
  source.freezePanes.freezeRows(1);
  source.tables.add(sourceRange, true, "TiteNewdesignSource");

  const comparison = workbook.worksheets.add("Comparison");
  comparison.showGridLines = false;
  comparison.getRange("A1:G1").merge();
  comparison.getRange("A1").values = [["Phase I/II Operating Characteristics"]];
  styleRange(comparison, "A1:G1", {
    fill: NAVY,
    font: { color: "#FFFFFF", bold: true, size: 16 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
  });
  comparison.getRange("A1:G1").format.rowHeight = 30;

  comparison.getRange("A2:G2").merge();
  comparison.getRange("A2").values = [["TITE versus new non-TITE design | N = 30 | Cycle_Max = 2 | efficacy futility cutoff = 0.85"]];
  styleRange(comparison, "A2:G2", {
    fill: "#EDF2F7",
    font: { color: NAVY, bold: true, size: 10 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
  });

  comparison.getRange("A4:G4").merge();
  comparison.getRange("A4").values = [["Comparison scope: matching two-stage/top-2-randomized configurations, with 1,000 simulations per available scenario."]];
  comparison.getRange("A5:G5").merge();
  comparison.getRange("A5").values = [["Sources: Presentation 8-17 TITE results (IDX 1001-2000) and new non-TITE results (IDX 0001-1000); all six scenarios are available."]];
  styleRange(comparison, "A4:G5", {
    fill: "#FFF7ED",
    font: { color: "#7C2D12", italic: true, size: 9 },
    wrapText: true,
    verticalAlignment: "center",
  });
  comparison.getRange("A4:G5").format.rowHeight = 22;

  comparison.getRange("A7:G7").values = [["Metric", "Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "No OBD %"]];
  styleRange(comparison, "A7:G7", {
    fill: TEAL,
    font: { color: "#FFFFFF", bold: true, size: 10 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    borders: { preset: "all", style: "thin", color: TEAL },
  });
  comparison.getRange("A7:G7").format.rowHeight = 22;

  function setFormulaRow(row, label, formulas, fill, bold = false, unavailable = false) {
    comparison.getRange(`A${row}`).values = [[label]];
    if (unavailable) {
      comparison.getRange(`B${row}:G${row}`).values = [["N/A", "N/A", "N/A", "N/A", "N/A", "N/A"]];
    } else {
      comparison.getRange(`B${row}:G${row}`).formulas = [formulas];
    }
    styleRange(comparison, `A${row}:G${row}`, {
      fill,
      font: { color: TEXT, bold, size: 10 },
      borders: { preset: "all", style: "thin", color: GRID },
      verticalAlignment: "center",
    });
    styleRange(comparison, `B${row}:G${row}`, { horizontalAlignment: "center" });
    comparison.getRange(`A${row}:G${row}`).format.rowHeight = 18;
  }

  function setSection(row, label, headers = null) {
    comparison.getRange(`A${row}`).values = [[label]];
    if (headers) comparison.getRange(`B${row}:D${row}`).values = [[...headers]];
    styleRange(comparison, `A${row}:G${row}`, {
      fill: GREY,
      font: { color: TEXT, bold: true, size: 10 },
      borders: { preset: "all", style: "thin", color: GRID },
      verticalAlignment: "center",
    });
    styleRange(comparison, `B${row}:G${row}`, { horizontalAlignment: "center" });
    comparison.getRange(`A${row}:G${row}`).format.rowHeight = 18;
  }

  for (let index = 0; index < SCENARIOS.length; index += 1) {
    const scenario = SCENARIOS[index];
    const start = 8 + (index * 13);
    const titeRow = sourceRow(TITE, scenario);
    const newRow = sourceRow(NEW, scenario);
    const titeAvailable = records.find(({ excelRow }) => excelRow === titeRow).record.Available === "True";

    comparison.getRange(`A${start}:G${start}`).merge();
    comparison.getRange(`A${start}`).values = [[`Scenario ${scenario}`]];
    styleRange(comparison, `A${start}:G${start}`, {
      fill: "#D9E2F3",
      font: { color: NAVY, bold: true, size: 11 },
      borders: { preset: "all", style: "thin", color: GRID },
      verticalAlignment: "center",
    });
    comparison.getRange(`A${start}:G${start}`).format.rowHeight = 20;

    comparison.getRange(`A${start + 1}`).values = [["DLT rate"]];
    comparison.getRange(`B${start + 1}:F${start + 1}`).formulas = [[1, 2, 3, 4, 5].map((dose) => field(newRow, `True_DLT_Dose_${dose}`))];
    comparison.getRange(`A${start + 2}`).values = [["Efficacy rate"]];
    comparison.getRange(`B${start + 2}:F${start + 2}`).formulas = [[1, 2, 3, 4, 5].map((dose) => field(newRow, `True_Efficacy_Dose_${dose}`))];
    styleRange(comparison, `A${start + 1}:G${start + 2}`, {
      font: { color: TEXT, size: 10 },
      borders: { preset: "all", style: "thin", color: GRID },
      verticalAlignment: "center",
    });
    styleRange(comparison, `B${start + 1}:G${start + 2}`, { horizontalAlignment: "center" });
    comparison.getRange(`A${start + 1}:G${start + 2}`).format.rowHeight = 18;
    comparison.getRange(`B${start + 1}:F${start + 2}`).format.numberFormat = "0.00";

    setSection(start + 3, "OBD selection (%)");
    setFormulaRow(
      start + 4,
      titeAvailable ? TITE : `${TITE} - unavailable`,
      [...[1, 2, 3, 4, 5].map((dose) => field(titeRow, `OBD_Selection_Dose_${dose}_pct`)), field(titeRow, "No_OBD_Selection_pct")],
      titeAvailable ? PALE_BLUE : PALE_RED,
      false,
      !titeAvailable,
    );
    setFormulaRow(
      start + 5,
      NEW,
      [...[1, 2, 3, 4, 5].map((dose) => field(newRow, `OBD_Selection_Dose_${dose}_pct`)), field(newRow, "No_OBD_Selection_pct")],
      PALE_GREEN,
      true,
    );
    comparison.getRange(`B${start + 4}:G${start + 5}`).format.numberFormat = "0.0";

    setSection(start + 6, "Mean administrations by dose");
    setFormulaRow(
      start + 7,
      titeAvailable ? TITE : `${TITE} - unavailable`,
      [...[1, 2, 3, 4, 5].map((dose) => field(titeRow, `Mean_Administrations_Dose_${dose}`)), ""],
      titeAvailable ? PALE_BLUE : PALE_RED,
      false,
      !titeAvailable,
    );
    setFormulaRow(
      start + 8,
      NEW,
      [...[1, 2, 3, 4, 5].map((dose) => field(newRow, `Mean_Administrations_Dose_${dose}`)), ""],
      PALE_GREEN,
      true,
    );
    comparison.getRange(`B${start + 7}:F${start + 8}`).format.numberFormat = "0.0";

    setSection(start + 9, "Run-level means", ["Total admin.", "Unique patients", "Duration days"]);
    const runFormulas = (row) => [
      field(row, "Mean_Total_Administrations"),
      field(row, "Mean_Total_Unique_Patients"),
      field(row, "Mean_Duration_days"),
      "",
      "",
      "",
    ];
    setFormulaRow(start + 10, titeAvailable ? TITE : `${TITE} - unavailable`, runFormulas(titeRow), titeAvailable ? PALE_BLUE : PALE_RED, false, !titeAvailable);
    setFormulaRow(start + 11, NEW, runFormulas(newRow), PALE_GREEN, true);
    comparison.getRange(`B${start + 10}:D${start + 11}`).format.numberFormat = "0.0";
  }

  comparison.getRange("A:A").format.columnWidth = 33;
  comparison.getRange("B:F").format.columnWidth = 12;
  comparison.getRange("G:G").format.columnWidth = 13;
  comparison.freezePanes.freezeRows(7);

  const usedComparison = comparison.getUsedRange();
  const errorCells = usedComparison.values.flat().filter((value) => typeof value === "string" && /^#(REF!|DIV\/0!|VALUE!|NAME\?|N\/A)$/.test(value));
  if (errorCells.length > 0) throw new Error(`Workbook contains formula errors: ${errorCells.join(", ")}`);

  const inspection = await workbook.inspect({
    kind: "region,formula",
    sheetId: "Comparison",
    range: "A1:G85",
    maxChars: 8000,
    tableMaxRows: 12,
    tableMaxCols: 7,
  });
  console.log(inspection.ndjson ?? inspection);

  const preview = await workbook.render({ sheetName: "Comparison", autoCrop: "all", scale: 1, format: "png" });
  await fs.writeFile(PREVIEW_PATH, new Uint8Array(await preview.arrayBuffer()));

  const output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(XLSX_PATH);
  console.log(XLSX_PATH);
  console.log(PREVIEW_PATH);
}

await main();
