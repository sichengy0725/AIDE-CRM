import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const projectRoot = process.cwd();
const outputDir = path.join(projectRoot, "outputs", "alpha_comparison_nontite_newdesign_boin_20260816");
const outputPath = path.join(outputDir, "nontite_newdesign_alpha_and_boin_comparison.xlsx");
const renderDir = path.join(projectRoot, "tmp", "nontite_alpha_comparison_workbook", "rendered_final");
await fs.mkdir(renderDir, { recursive: true });

const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(outputPath));
const sheets = await workbook.inspect({ kind: "sheet", include: "id,name" });
console.log(sheets.ndjson);
const selection = await workbook.inspect({
  kind: "table",
  range: "Alpha Selection!A7:M16",
  include: "values,formulas",
  tableMaxRows: 10,
  tableMaxCols: 13,
});
console.log(selection.ndjson);
const methods = await workbook.inspect({
  kind: "table",
  range: "Methods alpha 0!A7:S17",
  include: "values,formulas",
  tableMaxRows: 11,
  tableMaxCols: 19,
});
console.log(methods.ndjson);
const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 200 },
  summary: "post-export formula error scan",
});
console.log(errors.ndjson);

const renderSpecs = [
  ["Read Me", "A1:B14", "readme.png", 1.3],
  ["Alpha Selection", "A1:M34", "alpha_selection.png", 1.1],
  ["Alpha Allocation", "A1:N34", "alpha_allocation.png", 1.0],
  ["Methods alpha 0", "A1:S37", "methods_alpha0.png", 0.8],
  ["New-design source", "A1:Q24", "newdesign_source.png", 0.8],
  ["Comparator sources", "A1:O19", "comparator_sources.png", 1.0],
];
for (const [sheetName, range, filename, scale] of renderSpecs) {
  const image = await workbook.render({ sheetName, range, scale, format: "png" });
  await fs.writeFile(path.join(renderDir, filename), new Uint8Array(await image.arrayBuffer()));
}
console.log(`RENDER_DIR=${renderDir}`);
