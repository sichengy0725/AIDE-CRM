import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.cwd();
const xlsxPath = path.join(root, "Presentation 8-17-2026", "Table and Plots", "phase12_tite_vs_nontite_3designs_N30_fut0p85_scenarios_1_16_20_24_27_38_comparison.xlsx");
const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(xlsxPath));
const sheets = await workbook.inspect({ kind: "sheet", include: "id,name" });
const errors = await workbook.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 300 }, summary: "final formula error scan" });
console.log(sheets.ndjson ?? sheets);
console.log(errors.ndjson ?? errors);
const sourcePreview = await workbook.render({ sheetName: "Source Data", range: "A1:AI37", scale: 0.42, format: "png" });
const sourcePreviewPath = path.join(root, "tmp", "phase12_tite_nontite_3designs_workbook", "source_data_preview.png");
await fs.writeFile(sourcePreviewPath, new Uint8Array(await sourcePreview.arrayBuffer()));
console.log(sourcePreviewPath);
