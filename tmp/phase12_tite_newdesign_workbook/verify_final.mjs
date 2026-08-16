import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.cwd();
const xlsxPath = path.join(
  root,
  "Presentation 8-17-2026",
  "Table and Plots",
  "phase12_tite_vs_newdesign_N30_fut0p85_scenarios_1_16_20_24_27_38_comparison.xlsx",
);
const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(xlsxPath));
const sheetCheck = await workbook.inspect({ kind: "sheet", include: "id,name" });
const errorCheck = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 300 },
  summary: "final formula error scan",
});
console.log(sheetCheck.ndjson ?? sheetCheck);
console.log(errorCheck.ndjson ?? errorCheck);
const preview = await workbook.render({ sheetName: "Source Data", range: "A1:AD13", scale: 0.55, format: "png" });
const previewPath = path.join(root, "tmp", "phase12_tite_newdesign_workbook", "source_data_preview.png");
await fs.writeFile(previewPath, new Uint8Array(await preview.arrayBuffer()));
console.log(previewPath);
