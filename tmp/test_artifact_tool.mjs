import { Workbook } from "@oai/artifact-tool";

const workbook = Workbook.create();
const sheet = workbook.worksheets.add("Test");
sheet.getRange("A1:B2").values = [["a", "b"], [1, 2]];
console.log("artifact-tool basic workbook succeeded");
