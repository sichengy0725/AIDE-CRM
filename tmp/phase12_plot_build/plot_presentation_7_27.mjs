// Build comparison figures for the Presentation 7-27-2026 result files.
// The script is intentionally self-contained: it uses only the bundled Node
// runtime modules (pdfjs-dist and sharp) plus the result files in this repo.

import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import sharp from 'sharp';
import * as pdfjsLib from 'pdfjs-dist/legacy/build/pdf.mjs';
import { PDFDocument } from 'pdf-lib';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..', '..');
const inputDir = path.join(root, 'Presentation 7-27-2026');
const rawDir = path.join(inputDir, 'Raw Result');
const outDir = path.join(inputDir, 'Table and Plots', 'Prior');

function csvRows(text) {
  const rows = [];
  let row = [];
  let value = '';
  let quote = false;
  for (let i = 0; i < text.length; i += 1) {
    const char = text[i];
    if (char === '"') {
      if (quote && text[i + 1] === '"') { value += '"'; i += 1; }
      else quote = !quote;
    } else if (char === ',' && !quote) {
      row.push(value); value = '';
    } else if ((char === '\n' || char === '\r') && !quote) {
      if (char === '\r' && text[i + 1] === '\n') i += 1;
      row.push(value); value = '';
      if (row.length > 1 || row[0] !== '') rows.push(row);
      row = [];
    } else value += char;
  }
  if (value.length || row.length) { row.push(value); rows.push(row); }
  const [headers, ...body] = rows;
  return body.map((cells) => Object.fromEntries(headers.map((header, index) => [header, cells[index] ?? ''])));
}

function num(value) {
  if (value === undefined || value === null || value === '' || value === 'NA' || value === 'NaN') return Number.NaN;
  return Number(value);
}

function esc(text) {
  return String(text)
    .replaceAll('&', '&amp;')
    .replaceAll('<', '&lt;')
    .replaceAll('>', '&gt;')
    .replaceAll('"', '&quot;');
}

function tag(value) {
  return String(value).replace('.', 'p');
}

function fmt(value, digits = 1) {
  return Number.isFinite(value) ? Number(value).toFixed(digits).replace(/\.0$/, '') : 'NA';
}

async function readCsv(file) {
  return csvRows(await fs.readFile(file, 'utf8'));
}

async function writeText(file, text) {
  await fs.mkdir(path.dirname(file), { recursive: true });
  await fs.writeFile(file, text, 'utf8');
}

async function writeSvgOutputs(svg, stem) {
  const svgBuffer = Buffer.from(svg);
  await fs.mkdir(outDir, { recursive: true });
  await fs.writeFile(path.join(outDir, `${stem}.svg`), svgBuffer);
  const pngBuffer = await sharp(svgBuffer).png().toBuffer();
  await fs.writeFile(path.join(outDir, `${stem}.png`), pngBuffer);
  const document = await PDFDocument.create();
  const image = await document.embedPng(pngBuffer);
  const page = document.addPage([1296, 1008]);
  page.drawImage(image, { x: 0, y: 0, width: 1296, height: 1008 });
  await fs.writeFile(path.join(outDir, `${stem}.pdf`), await document.save());
}

async function mergePdfs(stems, outputName) {
  const merged = await PDFDocument.create();
  for (const stem of stems) {
    const source = await PDFDocument.load(await fs.readFile(path.join(outDir, `${stem}.pdf`)));
    const pages = await merged.copyPages(source, source.getPageIndices());
    pages.forEach((page) => merged.addPage(page));
  }
  await fs.writeFile(path.join(outDir, outputName), await merged.save());
}

function svgText(x, y, text, options = {}) {
  const { size = 18, anchor = 'start', weight = 400, fill = '#1a1a1a', rotate = null, family = 'Arial, Helvetica, sans-serif' } = options;
  const transform = rotate === null ? '' : ` transform="rotate(${rotate} ${x} ${y})"`;
  return `<text x="${x}" y="${y}" text-anchor="${anchor}" font-family="${family}" font-size="${size}" font-weight="${weight}" fill="${fill}"${transform}>${esc(text)}</text>`;
}

function parseNumberList(line, count) {
  const tokens = String(line).match(/-?(?:\d+\.?\d*|\.\d+)/g) ?? [];
  return tokens.slice(-count).map(Number);
}

async function pdfLines(pdfPath) {
  const data = new Uint8Array(await fs.readFile(pdfPath));
  const pdf = await pdfjsLib.getDocument({ data, disableWorker: true }).promise;
  const pages = [];
  for (let index = 1; index <= pdf.numPages; index += 1) {
    const page = await pdf.getPage(index);
    const content = await page.getTextContent();
    const items = content.items
      .filter((item) => item.str && item.str.trim())
      .map((item) => ({ text: item.str.trim(), x: item.transform[4], y: item.transform[5] }));
    const sorted = items.sort((a, b) => b.y - a.y || a.x - b.x);
    const groups = [];
    for (const item of sorted) {
      const existing = groups.find((group) => Math.abs(group.y - item.y) < 2.2);
      if (existing) existing.items.push(item);
      else groups.push({ y: item.y, items: [item] });
    }
    pages.push(groups
      .sort((a, b) => b.y - a.y)
      .map((group) => group.items.sort((a, b) => a.x - b.x).map((item) => item.text).join(' ')));
  }
  return pages;
}

function parsePdfOperatingCharacteristics(pageLines, source) {
  const scenarios = [];
  let current = null;
  for (const line of pageLines.flat()) {
    const scenarioMatch = line.match(/^Scenario\s+(\d+)$/i);
    if (scenarioMatch) {
      current = { scenario: Number(scenarioMatch[1]), source };
      scenarios.push(current);
    } else if (current && /^DLT\s+rate/i.test(line)) {
      current.dlt = parseNumberList(line, 5);
    } else if (current && /^Efficacy\s+rate/i.test(line)) {
      current.efficacy = parseNumberList(line, 5);
    } else if (current && /^Utility/i.test(line)) {
      current.utility = parseNumberList(line, 5);
    } else if (current && /^(?:No\.|#)\s*Pts\s+treated/i.test(line)) {
      current.treated = parseNumberList(line, 5);
    } else if (current && /^Select\s*%/i.test(line)) {
      const values = parseNumberList(line, 6);
      current.selected = values.slice(0, 5);
      current.stop = values[5] ?? 0;
    }
  }
  return scenarios;
}

function parseEffTox(html, source) {
  const rows = [];
  let current = null;
  const scenarioRe = /<tr><td\s+colspan\s*=\s*"8"><b>\s*(\d+)\s*<\/b><\/tr>/gi;
  const starts = [...html.matchAll(scenarioRe)];
  for (let index = 0; index < starts.length; index += 1) {
    const start = starts[index].index;
    const end = index + 1 < starts.length ? starts[index + 1].index : html.length;
    const block = html.slice(start, end);
    current = { scenario: Number(starts[index][1]), source };
    const cellsFor = (label) => {
      const found = block.search(label);
      if (found < 0) return [];
      const closing = block.indexOf('</tr>', found);
      const fragment = block.slice(Math.max(0, found - 120), closing >= 0 ? closing + 5 : block.length);
      return [...fragment.matchAll(/<td[^>]*>([\s\S]*?)<\/td>/gi)]
        .map((cell) => cell[1].replace(/<[^>]+>/g, ' ').replace(/&nbsp;/g, ' ').replace(/\s+/g, ' ').trim());
    };
    const trueCells = cellsFor(/True pT, pE/i);
    const trueValues = trueCells.slice(-6, -1).map((cell) => cell.split(',').map((entry) => Number(entry.trim())));
    current.dlt = trueValues.map((value) => value[0]);
    current.efficacy = trueValues.map((value) => value[1]);
    const utilityCells = cellsFor(/Trade-off Value/i);
    current.utility = utilityCells.slice(-6, -1).map(Number);
    const selectedCells = cellsFor(/% selected/i);
    const selected = selectedCells.slice(-6).map(Number);
    current.selected = selected.slice(0, 5);
    current.stop = selected[5] ?? 0;
    const treatedCells = cellsFor(/# Patients Treated/i);
    current.treated = treatedCells.slice(-6, -1).map(Number);
    rows.push(current);
  }
  return rows;
}

function simpleTable(rows, columns) {
  const header = columns.map((column) => `"${column.label}"`).join(',');
  const body = rows.map((row) => columns.map((column) => {
    const raw = typeof column.value === 'function' ? column.value(row) : row[column.value];
    return `"${String(raw ?? '').replaceAll('"', '""')}"`;
  }).join(','));
  return [header, ...body].join('\n');
}

function yDomain(values, min = 0, floorMax = 1) {
  const finite = values.filter(Number.isFinite);
  const max = Math.max(floorMax, ...finite);
  const roundedMax = max > 20 ? Math.ceil(max / 10) * 10 : Math.ceil(max / 2) * 2;
  return [min, roundedMax || floorMax];
}

function rangeForMetric(metric, values) {
  if (metric.key === 'mtdSelection') return [0, 100];
  if (metric.key === 'mtdAllocation') return [0, 100];
  if (metric.key === 'overdoseSelection') return [0, 100];
  if (metric.key === 'overdoseAllocation') return [0, 100];
  if (metric.key === 'sampleSize') return yDomain(values, 0, 25);
  if (metric.key === 'duration') return yDomain(values, 0, 10);
  if (metric.key === 'absoluteDifference') return yDomain(values, 0, 1);
  return yDomain(values, 0, 1);
}

const palette = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#56B4E9', '#E69F00', '#7F3C8D', '#4E79A7', '#59A14F'];

function drawBarPanel({ x, y, width, height, title, yLabel, categories, series, domain }) {
  const left = x + 76;
  const right = x + width - 22;
  const top = y + 34;
  const bottom = y + height - 45;
  const [min, max] = domain;
  const span = max - min || 1;
  const plotX = (idx) => left + ((idx + 0.5) * (right - left)) / categories.length;
  const plotY = (value) => bottom - ((value - min) / span) * (bottom - top);
  const ticks = 5;
  let svg = '';
  svg += `<rect x="${x}" y="${y}" width="${width}" height="${height}" fill="#EBEBEB"/>`;
  for (let tick = 0; tick <= ticks; tick += 1) {
    const value = min + (span * tick) / ticks;
    const yy = plotY(value);
    svg += `<line x1="${left}" y1="${yy}" x2="${right}" y2="${yy}" stroke="#ffffff" stroke-width="1.2"/>`;
    svg += svgText(left - 11, yy + 5, fmt(value, value < 1 ? 2 : 0), { size: 15, anchor: 'end', fill: '#4D4D4D' });
  }
  if (min < 0 && max > 0) svg += `<line x1="${left}" y1="${plotY(0)}" x2="${right}" y2="${plotY(0)}" stroke="#4D4D4D"/>`;
  const clusterWidth = (right - left) / categories.length;
  const usable = clusterWidth * 0.8;
  const barWidth = usable / series.length;
  series.forEach((line, seriesIndex) => {
    line.values.forEach((value, categoryIndex) => {
      if (!Number.isFinite(value)) return;
      const bx = plotX(categoryIndex) - usable / 2 + seriesIndex * barWidth;
      const by = plotY(Math.max(value, 0));
      const zero = plotY(0);
      svg += `<rect x="${bx}" y="${Math.min(by, zero)}" width="${Math.max(1, barWidth - 0.8)}" height="${Math.abs(zero - by)}" fill="${line.color}"/>`;
    });
  });
  categories.forEach((category, index) => {
    svg += svgText(plotX(index), bottom + 23, category, { size: 15, anchor: 'middle', fill: '#4D4D4D' });
  });
  svg += svgText((left + right) / 2, y + 21, title, { size: 18, anchor: 'middle', weight: 400 });
  svg += svgText(x + 20, (top + bottom) / 2, yLabel, { size: 15, anchor: 'middle', fill: '#333333', rotate: -90 });
  return svg;
}

function drawLinePanel({ x, y, width, height, title, yLabel, categories, series, domain }) {
  const left = x + 86;
  const right = x + width - 26;
  const top = y + 34;
  const bottom = y + height - 48;
  const [min, max] = domain;
  const span = max - min || 1;
  const plotX = (idx) => left + ((idx + 0.5) * (right - left)) / categories.length;
  const plotY = (value) => bottom - ((value - min) / span) * (bottom - top);
  let svg = `<rect x="${x}" y="${y}" width="${width}" height="${height}" fill="#EBEBEB"/>`;
  for (let tick = 0; tick <= 5; tick += 1) {
    const value = min + (span * tick) / 5;
    const yy = plotY(value);
    svg += `<line x1="${left}" y1="${yy}" x2="${right}" y2="${yy}" stroke="#ffffff" stroke-width="1.2"/>`;
    svg += svgText(left - 11, yy + 5, fmt(value, value < 1 ? 2 : 0), { size: 13, anchor: 'end', fill: '#4D4D4D' });
  }
  categories.forEach((category, index) => {
    const xx = plotX(index);
    svg += `<line x1="${xx}" y1="${top}" x2="${xx}" y2="${bottom}" stroke="#ffffff" stroke-width="0.8"/>`;
    svg += svgText(xx, bottom + 21, category, { size: 13, anchor: 'middle', fill: '#4D4D4D' });
  });
  for (const line of series) {
    const points = line.values.map((value, index) => Number.isFinite(value) ? `${plotX(index)},${plotY(value)}` : null);
    let segment = [];
    for (const point of [...points, null]) {
      if (point) segment.push(point);
      else if (segment.length) {
        svg += `<polyline fill="none" stroke="${line.color}" stroke-width="2.2" points="${segment.join(' ')}"/>`;
        segment = [];
      }
    }
    line.values.forEach((value, index) => {
      if (Number.isFinite(value)) svg += `<circle cx="${plotX(index)}" cy="${plotY(value)}" r="3.1" fill="${line.color}"/>`;
    });
  }
  svg += svgText((left + right) / 2, y + 21, title, { size: 17, anchor: 'middle', weight: 400 });
  svg += svgText(x + 22, (top + bottom) / 2, yLabel, { size: 14, anchor: 'middle', fill: '#333333', rotate: -90 });
  return svg;
}

function drawLegend({ labels, colors, x, y, width, fontSize = 17 }) {
  const totalChars = labels.reduce((sum, label) => sum + label.length + 3, 0);
  let cx = x;
  let cy = y;
  const unit = Math.max(85, Math.floor(width / Math.max(1, totalChars / 8)));
  let svg = '';
  labels.forEach((label, index) => {
    const itemWidth = Math.max(150, label.length * fontSize * 0.56 + 36);
    if (cx + itemWidth > x + width && cx > x) { cx = x; cy += fontSize + 12; }
    svg += `<rect x="${cx}" y="${cy - 13}" width="17" height="17" fill="${colors[index]}"/>`;
    svg += svgText(cx + 24, cy + 1, label, { size: fontSize, fill: '#222222' });
    cx += Math.max(unit, itemWidth + 18);
  });
  return svg;
}

function figureSvg({ title, subtitle, panels, legend, width = 2592, height = 2016, type = 'bar' }) {
  const sevenPanels = panels.length > 6;
  const margin = { left: 85, right: 85, top: 125, bottom: sevenPanels ? 205 : 160 };
  const cols = sevenPanels ? 4 : 3;
  const rows = Math.ceil(panels.length / cols);
  const gapX = 44;
  const gapY = 45;
  const panelWidth = (width - margin.left - margin.right - gapX * (cols - 1)) / cols;
  const panelHeight = (height - margin.top - margin.bottom - gapY * (rows - 1)) / rows;
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}"><rect width="100%" height="100%" fill="white"/>`;
  svg += svgText(width / 2, 56, title, { size: 34, anchor: 'middle', weight: 700 });
  if (subtitle) svg += svgText(width / 2, 88, subtitle, { size: 20, anchor: 'middle', fill: '#4D4D4D' });
  panels.forEach((panel, index) => {
    const col = index % cols;
    const row = Math.floor(index / cols);
    const x = margin.left + col * (panelWidth + gapX);
    const y = margin.top + row * (panelHeight + gapY);
    svg += (type === 'line' ? drawLinePanel : drawBarPanel)({ x, y, width: panelWidth, height: panelHeight, ...panel });
  });
  svg += drawLegend({ labels: legend.labels, colors: legend.colors, x: margin.left + 50, y: height - (sevenPanels ? 125 : 75), width: width - margin.left - margin.right - 100, fontSize: legend.fontSize ?? 17 });
  return `${svg}</svg>`;
}

function metricRowsForRandom(data, methodLabel, { oracle = false } = {}) {
  const byMetric = Object.fromEntries(data.map((row) => [row.Metric, row]));
  const totalDoseAdministrations = [1, 2, 3, 4, 5]
    .map((dose) => num(byMetric['Average dose allocation']?.[`D${dose}`]))
    .reduce((sum, value) => sum + value, 0);
  const sampleSize = oracle
    ? totalDoseAdministrations
    : num(byMetric['Average total unique patients']?.Total);
  return {
    methodLabel,
    mtdSelection: num(byMetric['Average MTD selection %']?.Total),
    mtdAllocation: num(byMetric['Average MTD allocation']?.Total) * 100 / num(data[0]?.Nmax_eff),
    overdoseSelection: num(byMetric['Average overdose selection %']?.Total),
    overdoseAllocation: num(byMetric['Average overdose allocation']?.Total) * 100 / num(data[0]?.Nmax_eff),
    sampleSize,
    duration: oracle ? (56 * sampleSize) / 7 : num(byMetric['Average trial duration']?.Duration) / 7,
    absoluteDifference: Math.abs(num(
      byMetric['Mean paired absolute MTD selection difference vs oracle %']?.Total ??
      byMetric['Mean paired MTD selection difference vs oracle %']?.Total
    )),
  };
}

async function buildRandomFigures() {
  const inputFiles = {
    '0.05': 'randomsce_targetgap0p05_batch10_j1to1000_n100_table_summary.csv',
    '0.10': 'randomsce_targetgap0p10_batch10_j1to1000_n100_table_summary.csv',
    '0.15': 'randomsce_targetgap0p15_batch10_j1to1000_n100_table_summary.csv',
  };
  const gapRows = {};
  for (const [gap, file] of Object.entries(inputFiles)) gapRows[gap] = await readCsv(path.join(rawDir, file));
  const priors = [
    { key: '0.15|0.85', label: 'Beta(0.15, 0.85)', color: palette[0] },
    { key: '0.3|0.7', label: 'Beta(0.3, 0.7)', color: palette[1] },
    { key: '0.5|0.5', label: 'Beta(0.5, 0.5)', color: palette[2] },
    { key: '1|1', label: 'Beta(1, 1)', color: palette[3] },
  ];
  const alphas = [0, 0.3, 0.6, 0.9];
  const metricSpecs = [
    { key: 'mtdSelection', panel: 'A', title: 'MTD selection', yLabel: 'Percentage (%)' },
    { key: 'mtdAllocation', panel: 'B', title: 'MTD allocation', yLabel: 'Percentage (%)' },
    { key: 'overdoseSelection', panel: 'C', title: 'Overdose selection', yLabel: 'Percentage (%)' },
    { key: 'overdoseAllocation', panel: 'D', title: 'Overdose allocation', yLabel: 'Percentage (%)' },
    { key: 'sampleSize', panel: 'E', title: 'Average sample size', yLabel: 'Patients' },
    { key: 'duration', panel: 'F', title: 'Average trial duration', yLabel: 'Weeks' },
    { key: 'absoluteDifference', panel: 'G', title: 'Absolute MTD-selection change vs oracle', yLabel: 'Percentage points' },
  ];
  const oracle = {
    key: 'oracle', label: 'Oracle: r-fixed, alpha = 0', color: '#222222', oracle: true,
    matches: (row) => row.Model === 'CRM' && row.CRM_r_model === 'r_fixed' && Math.abs(num(row.Alpha_true)) < 1e-12,
  };
  const priorMethods = priors.map((prior) => ({
    key: `random-${prior.key}`,
    label: `Discount CRM ${prior.label.replace('Beta', '')}`,
    color: prior.color,
    matches: (row, alpha) => row.Model === 'CRM' && row.CRM_r_model === 'random' &&
      Math.abs(num(row.Alpha_true) - alpha) < 1e-12 && `${row.CRM_r_Prior_a}|${row.CRM_r_Prior_b}` === prior.key,
  }));
  const allMethods = [
    { key: 'r-fixed', label: 'r-fixed CRM', color: palette[4], matches: (row, alpha) => row.Model === 'CRM' && row.CRM_r_model === 'r_fixed' && Math.abs(num(row.Alpha_true) - alpha) < 1e-12 },
    { key: 'alpha-crm', label: 'Alpha-CRM', color: palette[5], matches: (row, alpha) => row.Model === 'CRM' && row.CRM_r_model === 'alpha_crm' && Math.abs(num(row.Alpha_true) - alpha) < 1e-12 },
    { key: 'cumu-crm', label: 'Cumulative CRM', color: palette[6], matches: (row, alpha) => row.Model === 'CRM' && row.CRM_r_model === 'cumu_crm' && Math.abs(num(row.Alpha_true) - alpha) < 1e-12 },
    ...priorMethods,
  ];

  const selectRows = (rows, method, alpha, gap) => {
    const selected = rows.filter((row) => method.matches(row, alpha));
    if (selected.length !== 12) {
      throw new Error(`Expected 12 metric rows for alpha=${alpha}, gap=${gap}, ${method.label}; found ${selected.length}.`);
    }
    return selected;
  };

  const buildFigureSet = async ({ figure, methods, stemPrefix, title, subtitle, legendFontSize }) => {
    const extracted = [];
    const stems = [];
    for (const alpha of alphas) {
      const values = {};
      for (const method of methods) {
        values[method.key] = {};
        for (const [gap, rows] of Object.entries(gapRows)) {
          const selectedRows = selectRows(rows, method, alpha, gap);
          const summary = metricRowsForRandom(selectedRows, method.label, { oracle: Boolean(method.oracle) });
          values[method.key][gap] = summary;
          for (const metric of metricSpecs) {
            extracted.push({ Figure: figure, Alpha: alpha, Gap: gap, Method: method.label, Metric: metric.title, Value: summary[metric.key], Scenarios_found: selectedRows[0].n_scenarios_found, Derived_oracle_quantity: Boolean(method.oracle && (metric.key === 'sampleSize' || metric.key === 'duration')) });
          }
        }
      }
      const panels = metricSpecs.map((metric) => {
        const series = methods.map((method) => ({ label: method.label, color: method.color, values: Object.keys(inputFiles).map((gap) => values[method.key][gap][metric.key]) }));
        return { title: `(${metric.panel}) ${metric.title}`, yLabel: metric.yLabel, categories: Object.keys(inputFiles), series, domain: rangeForMetric(metric, series.flatMap((line) => line.values)) };
      });
      const stem = `${stemPrefix}_alpha${tag(alpha)}`;
      stems.push(stem);
      await writeSvgOutputs(figureSvg({
        title: `${title} (alpha = ${alpha})`,
        subtitle,
        panels,
        legend: { labels: methods.map((method) => method.label), colors: methods.map((method) => method.color), fontSize: legendFontSize },
      }), stem);
    }
    await mergePdfs(stems, `${stemPrefix}_all_alphas.pdf`);
    await writeText(path.join(outDir, `${stemPrefix}_extracted_data.csv`), simpleTable(extracted, [
      { label: 'Figure', value: 'Figure' }, { label: 'Alpha', value: 'Alpha' }, { label: 'Gap', value: 'Gap' }, { label: 'Method', value: 'Method' }, { label: 'Metric', value: 'Metric' }, { label: 'Value', value: 'Value' }, { label: 'Scenarios found', value: 'Scenarios_found' }, { label: 'Derived oracle quantity', value: 'Derived_oracle_quantity' },
    ]));
  };

  await buildFigureSet({
    figure: 'Random prior sensitivity',
    methods: [...priorMethods, oracle],
    stemPrefix: 'random_prior_sensitivity',
    title: 'Random scenarios - Discount CRM prior sensitivity',
    subtitle: 'Target gaps: 0.05, 0.10, and 0.15. Oracle = r-fixed CRM at alpha = 0; its duration is 56 days per administration.',
    legendFontSize: 16,
  });
  await buildFigureSet({
    figure: 'All methods and random priors',
    methods: [...allMethods, oracle],
    stemPrefix: 'random_all_methods_and_priors',
    title: 'Random scenarios - all methods and Discount CRM priors',
    subtitle: 'CFO excluded. Oracle = r-fixed CRM at alpha = 0; its duration is 56 days per administration.',
    legendFontSize: 13,
  });
}

function weightedMetric(values, weights) {
  const total = weights.reduce((sum, value) => sum + value, 0);
  return total > 0 ? values.reduce((sum, value, index) => sum + value * weights[index], 0) / total : Number.NaN;
}

function buildPhaseMetrics(record, truth, target = 0.3) {
  const treated = record.treated;
  const selected = record.selected;
  const obdIndex = truth.obd - 1;
  const totalTreated = treated.reduce((sum, value) => sum + value, 0);
  return {
    correctSelection: selected[obdIndex] ?? Number.NaN,
    overdoseSelection: selected.reduce((sum, value, index) => sum + (truth.dlt[index] > target + 1e-12 ? value : 0), 0),
    obdAllocation: totalTreated > 0 ? 100 * treated[obdIndex] / totalTreated : Number.NaN,
    dltRate: weightedMetric(truth.dlt, treated),
    efficacyRate: weightedMetric(truth.efficacy, treated),
    utility: weightedMetric(truth.utility, treated),
  };
}

function findAideRecord(rows, scenario, nmax, allocation, truth) {
  const records = rows.filter((row) => Number(row.Scenario) === scenario && Number(row.Nmax) === nmax && row.Allocation === allocation && Number(row.Utility_Type) === 3 && Math.abs(num(row.Lambda_T) - 1) < 1e-12);
  if (records.length !== 5) throw new Error(`Expected five AIDE rows for scenario=${scenario}, N=${nmax}, allocation=${allocation}; found ${records.length}.`);
  const ordered = [...records].sort((a, b) => Number(a.Dose) - Number(b.Dose));
  const selected = ordered.map((row) => num(row.OBD_Selection_pct));
  const treated = ordered.map((row) => num(row.Pts_Treated));
  return { scenario, source: allocation === 'one_stage' ? 'AIDE phase I' : 'AIDE two-stage', selected, treated, stop: num(ordered[0].No_OBD_Selection_pct), ...buildPhaseMetrics({ selected, treated }, truth) };
}

async function buildPhaseFigures() {
  const truthRows = await readCsv(path.join(inputDir, 'Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv'));
  const truth = Object.fromEntries(truthRows.map((row) => [Number(row.Scenario), {
    dlt: [1, 2, 3, 4, 5].map((dose) => num(row[`Tox_Dose${dose}`])),
    efficacy: [1, 2, 3, 4, 5].map((dose) => num(row[`Eff_Dose${dose}`])),
    utility: [1, 2, 3, 4, 5].map((dose) => num(row[`Utility3_Dose${dose}`])),
    obd: Number(row.OBD_Level_Utility3),
  }]));
  if (Object.keys(truth).length !== 37) throw new Error('The truth file must contain scenarios 1-37 exactly.');
  const aideRows = await readCsv(path.join(inputDir, 'AIDE_phase_I_II_IDX_1001_to_2000_dose_summary.csv'));
  const scenarios = Array.from({ length: 37 }, (_, index) => index + 1);
  const efftox30 = parseEffTox(await fs.readFile(path.join(inputDir, 'EffTox N = 30.html'), 'utf8'), 'EffTox');
  const efftox60 = parseEffTox(await fs.readFile(path.join(inputDir, 'EffTox N = 60.html'), 'utf8'), 'EffTox');
  const boin30 = parsePdfOperatingCharacteristics(await pdfLines(path.join(inputDir, 'BOIN12 Result N = 30.pdf')), 'BOIN12');
  const boin60 = parsePdfOperatingCharacteristics(await pdfLines(path.join(inputDir, 'BOIN12 Result N = 60.pdf')), 'BOIN12');
  const uboin30 = [];
  const uboin60 = [];
  for (const [range, count, offset] of [['Sce1-10', 10, 0], ['Sce11-20', 10, 10], ['Sce21-30', 10, 20], ['Sce31-37', 7, 30]]) {
    const n30 = path.join(inputDir, `U-BOIN ${range} N = 30 S1 = 9.pdf`);
    const n60 = path.join(inputDir, `U-BOIN ${range} N = 60 S1 = 12.pdf`);
    const parsed30 = parsePdfOperatingCharacteristics(await pdfLines(n30), 'U-BOIN');
    const parsed60 = parsePdfOperatingCharacteristics(await pdfLines(n60), 'U-BOIN');
    if (parsed30.length !== count || parsed60.length !== count) throw new Error(`Unexpected U-BOIN scenario count in ${range}.`);
    parsed30.forEach((record) => { record.scenario += offset; });
    parsed60.forEach((record) => { record.scenario += offset; });
    uboin30.push(...parsed30); uboin60.push(...parsed60);
  }
  await writeText(path.join(root, 'tmp', 'phase12_plot_build', 'pdf_parse_debug.json'), JSON.stringify({ boin30, boin60, uboin30, uboin60 }, null, 2));
  const sourceSets = [
    { n: 30, efftox: efftox30, boin: boin30, uboin: uboin30 },
    { n: 60, efftox: efftox60, boin: boin60, uboin: uboin60 },
  ];
  const metricSpecs = [
    { key: 'correctSelection', panel: 'A', title: 'True OBD selection', yLabel: 'Percentage (%)', domain: [0, 100] },
    { key: 'overdoseSelection', panel: 'B', title: 'Overdose selection', yLabel: 'Percentage (%)', domain: [0, 100] },
    { key: 'obdAllocation', panel: 'C', title: 'Patients allocated to true OBD', yLabel: 'Percentage (%)', domain: [0, 100] },
    { key: 'dltRate', panel: 'D', title: 'Treated-patient DLT rate', yLabel: 'Probability', domain: [0, 0.5] },
    { key: 'efficacyRate', panel: 'E', title: 'Treated-patient efficacy rate', yLabel: 'Probability', domain: [0, 0.65] },
    { key: 'utility', panel: 'F', title: 'Treated-patient utility', yLabel: 'Utility', domain: [0, 70] },
  ];
  const allExtracted = [];
  for (const set of sourceSets) {
    const phaseStems = [];
    const sourceLookup = {
      BOIN12: Object.fromEntries(set.boin.map((record) => [record.scenario, record])),
      'U-BOIN': Object.fromEntries(set.uboin.map((record) => [record.scenario, record])),
      EffTox: Object.fromEntries(set.efftox.map((record) => [record.scenario, record])),
    };
    for (const [name, lookup] of Object.entries(sourceLookup)) {
      if (Object.keys(lookup).length !== 37) throw new Error(`${name} N=${set.n} did not yield all 37 scenarios: ${Object.keys(lookup).join(', ')}`);
    }
    const methods = [
      { label: 'BOIN12', color: palette[0], values: {} },
      { label: 'U-BOIN', color: palette[1], values: {} },
      { label: 'EffTox', color: palette[2], values: {} },
      { label: 'AIDE phase I', color: palette[3], values: {} },
      { label: 'AIDE two-stage', color: palette[4], values: {} },
    ];
    for (const scenario of scenarios) {
      methods[0].values[scenario] = buildPhaseMetrics(sourceLookup.BOIN12[scenario], truth[scenario]);
      methods[1].values[scenario] = buildPhaseMetrics(sourceLookup['U-BOIN'][scenario], truth[scenario]);
      methods[2].values[scenario] = buildPhaseMetrics(sourceLookup.EffTox[scenario], truth[scenario]);
      methods[3].values[scenario] = findAideRecord(aideRows, scenario, set.n, 'one_stage', truth[scenario]);
      methods[4].values[scenario] = findAideRecord(aideRows, scenario, set.n, 'two_stage', truth[scenario]);
      for (const method of methods) for (const metric of metricSpecs) allExtracted.push({ Nmax: set.n, Scenario: scenario, Method: method.label, Metric: metric.title, Value: method.values[scenario][metric.key] });
    }
    for (const group of [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24], [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]]) {
      const panels = metricSpecs.map((metric) => ({
        title: `(${metric.panel}) ${metric.title}`,
        yLabel: metric.yLabel,
        categories: group.map(String),
        series: methods.map((method) => ({ label: method.label, color: method.color, values: group.map((scenario) => method.values[scenario][metric.key]) })),
        domain: metric.domain,
      }));
      const stem = `phase12_all_methods_N${set.n}_scenarios_${group[0]}_to_${group[group.length - 1]}`;
      phaseStems.push(stem);
      await writeSvgOutputs(figureSvg({
        title: `Phase I/II operating characteristics by scenario (N = ${set.n})`,
        subtitle: `Scenarios ${group[0]}-${group[group.length - 1]}; true OBD and dose-level utility from the lambda 1 truth summary.`,
        panels,
        legend: { labels: methods.map((method) => method.label), colors: methods.map((method) => method.color), fontSize: 17 },
        type: 'line',
      }), stem);
    }
    await mergePdfs(phaseStems, `phase12_all_methods_N${set.n}_scenarios_1_to_37.pdf`);
  }
  await writeText(path.join(outDir, 'phase12_all_methods_extracted_data.csv'), simpleTable(allExtracted, [
    { label: 'Nmax', value: 'Nmax' }, { label: 'Scenario', value: 'Scenario' }, { label: 'Method', value: 'Method' }, { label: 'Metric', value: 'Metric' }, { label: 'Value', value: 'Value' },
  ]));
}

async function main() {
  await fs.mkdir(outDir, { recursive: true });
  await buildRandomFigures();
  await writeText(path.join(outDir, 'README.txt'), [
    'Generated 2026-07-27 from the Presentation 7-27-2026 Raw Result files.',
    '',
    'Random-prior figures: four seven-panel figures (alpha 0, 0.3, 0.6, 0.9) comparing the four Beta priors for Discount CRM and the r-fixed alpha=0 oracle.',
    'All-method figures: four seven-panel figures (same alphas) comparing r-fixed CRM, Alpha-CRM, cumulative CRM, and all four Discount CRM prior settings, plus the oracle; CFO is excluded.',
    'Panel G is the mean paired absolute MTD-selection difference from the oracle.',
    'For the oracle only, sample size is the sum of the five average dose allocations (regular plus IPDE administrations) and duration is 56 days times that sample size, displayed in weeks.',
    '',
  ].join('\n'));
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
