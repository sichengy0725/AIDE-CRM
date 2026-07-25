import fs from "node:fs/promises";
import { Workbook } from "@oai/artifact-tool";

const templatePath = "BOIN12_scenario_template.csv";
const sourcePath = "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv";
const outputDir = "outputs/boin12_scenario_template_filled";
const outputPath = `${outputDir}/BOIN12_scenario_template_filled.csv`;

const templateCsv = await fs.readFile(templatePath, "utf8");
const sourceCsv = await fs.readFile(sourcePath, "utf8");
const templateWorkbook = await Workbook.fromCSV(templateCsv, { sheetName: "BOIN12 scenarios" });
const sourceWorkbook = await Workbook.fromCSV(sourceCsv, { sheetName: "Source summary" });

const templateSheet = templateWorkbook.worksheets.getItem("BOIN12 scenarios");
const sourceSheet = sourceWorkbook.worksheets.getItem("Source summary");
const sourceValues = sourceSheet.getUsedRange().values;
const headers = sourceValues[0];
const headerIndex = new Map(headers.map((header, index) => [String(header), index]));
const probabilityColumns = [
  "Tox_Dose1", "Tox_Dose2", "Tox_Dose3", "Tox_Dose4", "Tox_Dose5",
  "Eff_Dose1", "Eff_Dose2", "Eff_Dose3", "Eff_Dose4", "Eff_Dose5",
];

for (const column of ["Scenario", ...probabilityColumns]) {
  if (!headerIndex.has(column)) {
    throw new Error(`Source summary is missing required column: ${column}`);
  }
}

const scenarioRows = sourceValues.slice(1).map((row) => {
  const scenario = row[headerIndex.get("Scenario")];
  const toxicity = probabilityColumns.slice(0, 5).map((column) => row[headerIndex.get(column)]);
  const efficacy = probabilityColumns.slice(5).map((column) => row[headerIndex.get(column)]);
  if (toxicity.some((value) => value === null || value === "") || efficacy.some((value) => value === null || value === "")) {
    throw new Error(`Scenario ${scenario} has missing toxicity or efficacy probabilities.`);
  }
  return { scenario, toxicity, efficacy };
});

const outputValues = [[" ", "d1", "d2", "d3", "d4", "d5"]];
for (const { scenario, toxicity, efficacy } of scenarioRows) {
  outputValues.push([`Scenario ${scenario}`, " ", " ", " ", " ", " "]);
  outputValues.push(["DLT probability", ...toxicity]);
  outputValues.push(["Efficacy probability", ...efficacy]);
}

templateSheet.getRange("A1:F200").clear({ applyTo: "contents" });
templateSheet.getRange(`A1:F${outputValues.length}`).values = outputValues;
templateSheet.getRange(`A1:F${outputValues.length}`).format.wrapText = false;
templateSheet.getRange("A1:F1").format = {
  fill: "#D9EAF7",
  font: { bold: true },
  horizontalAlignment: "center",
};
templateSheet.getRange("A1").format.font = { bold: true, color: "#D9EAF7" };
for (let row = 2; row <= outputValues.length; row += 3) {
  templateSheet.getRange(`A${row}:F${row}`).format = {
    fill: "#F2F2F2",
    font: { bold: true },
  };
  templateSheet.getRange(`B${row}:F${row}`).format.font = { bold: true, color: "#F2F2F2" };
  templateSheet.getRange(`B${row + 1}:F${row + 2}`).format.numberFormat = "0.00";
}
templateSheet.getRange(`A1:A${outputValues.length}`).format.columnWidth = 22;
templateSheet.getRange(`B1:F${outputValues.length}`).format.columnWidth = 11;
templateSheet.showGridLines = false;

const check = await templateWorkbook.inspect({
  kind: "table",
  range: `BOIN12 scenarios!A1:F${outputValues.length}`,
  include: "values",
  tableMaxRows: 8,
  tableMaxCols: 6,
});
console.log(check.ndjson);

const preview = await templateWorkbook.render({
  sheetName: "BOIN12 scenarios",
  range: "A1:F10",
  scale: 2,
  format: "png",
});

await fs.mkdir(outputDir, { recursive: true });
await fs.writeFile(`${outputDir}/BOIN12_scenario_template_preview.png`, new Uint8Array(await preview.arrayBuffer()));

const csvEscape = (value) => {
  const text = value === null || value === undefined ? "" : String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
};
const csv = `${outputValues.map((row) => row.map(csvEscape).join(",")).join("\r\n")}\r\n`;
await fs.writeFile(outputPath, csv, "utf8");

console.log(JSON.stringify({ outputPath, scenarioCount: scenarioRows.length, rowCount: outputValues.length }));
