const pptxgen = require("pptxgenjs");
const path = require("path");
const { applyTheme } = require("/root/.claude/skills/synced/865e8771-0238-452b-8e45-d6a2b77f2cfe_d7a95306-8d5f-463b-b5b9-558348628557/pptx/scripts/apply_theme.js");

const D = process.env.POSTER_ASSETS || path.join(__dirname, "assets");
const img = (f) => path.join(D, f);
const SLATE = "23343B", TEAL = "00857C", TERRA = "D4643F", BG = "EEF3F3", CARD = "FFFFFF", MUTED = "5B6B73";
const FONT = "Calibri";

const THEME = {
  name: "Inversion poster", headFontFace: FONT, bodyFontFace: FONT,
  colors: { dk1: "000000", lt1: "FFFFFF", dk2: SLATE, lt2: BG, accent1: TEAL, accent2: TERRA, accent3: "8A9BA3",
            accent4: "C3CED3", accent5: "2F3E46", accent6: "E3E9EB", hlink: TEAL, folHlink: "00635C" },
};

const pres = new pptxgen();
pres.defineLayout({ name: "POSTER", width: 48, height: 27 });
pres.layout = "POSTER";
pres.theme = { headFontFace: FONT, bodyFontFace: FONT };
pres.title = "Do bigger bacterial genomes carry bigger inversions?";
pres.author = "Leen Al Radi, Mathias Azevedo, Abishek Kumar, Wout Simoens";

pres.addSection({ title: "Poster" });
const s = pres.addSlide({ sectionTitle: "Poster" });
s.background = { color: BG };

const shadow = () => ({ type: "outer", blur: 8, offset: 2, angle: 90, color: "000000", opacity: 0.10 });
const card = (x, y, w, h, name) =>
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.35, fill: { color: CARD },
    line: { color: "DCE4E6", width: 1 }, shadow: shadow(), objectName: name });
const header = (x, y, num, title, sub) => {
  s.addShape(pres.shapes.OVAL, { x, y, w: 0.95, h: 0.95, fill: { color: TEAL }, line: { color: TEAL },
    objectName: `badge ${num}` });
  s.addText(String(num), { x, y, w: 0.95, h: 0.95, align: "center", valign: "middle", fontSize: 30, bold: true,
    color: "FFFFFF", margin: 0, isTextBox: true, objectName: `badge ${num} text` });
  s.addText(title, { x: x + 1.2, y: y - 0.05, w: 14, h: 1.05, fontSize: 38, bold: true, color: SLATE,
    valign: "middle", margin: 0, isTextBox: true, objectName: `${title} heading` });
  if (sub) s.addText(sub, { x: x + 1.2 + 9.0, y: y, w: 21.6, h: 0.95, fontSize: 22, color: MUTED, align: "right",
    valign: "middle", margin: 0, isTextBox: true, objectName: `${title} subtitle` });
};
const headline = (x, y, w, letter, text, name) =>
  s.addText([{ text: letter + "  ", options: { color: TEAL, bold: true } },
             { text, options: { color: SLATE, bold: true } }],
    { x, y, w, h: 0.62, fontSize: 25, margin: 0, valign: "top", isTextBox: true, objectName: name });
const caption = (x, y, w, h, runs, name) =>
  s.addText(runs, { x, y, w, h, fontSize: 18, color: MUTED, margin: 0, valign: "top", isTextBox: true,
    objectName: name });

// ---------------------------------------------------------------- header
s.addShape(pres.shapes.RECTANGLE, { x: 0, y: 0, w: 48, h: 4.1, fill: { color: SLATE }, line: { color: SLATE },
  objectName: "title band" });
s.addImage({ path: img("logo.png"), x: 0.7, y: 0.95, w: 6.0, h: 6.0 * 157 / 438, objectName: "KU Leuven logo" });
s.addText("Do bigger bacterial genomes carry bigger inversions?", { x: 7.4, y: 0.35, w: 40, h: 1.7, fontSize: 84,
  bold: true, color: "FFFFFF", margin: 0, valign: "middle", isTextBox: true, objectName: "title" });
s.addText("Leen Al Radi  ·  Mathias Azevedo  ·  Abishek Kumar  ·  Wout Simoens", { x: 7.4, y: 2.15, w: 40, h: 0.75,
  fontSize: 32, color: "FFFFFF", margin: 0, valign: "middle", isTextBox: true, objectName: "authors" });
s.addText("Master of Bioinformatics, Faculty of Bioscience Engineering, KU Leuven   |   Coach: Juliette Luiselli   |   Integrated Bioinformatics Project 2026",
  { x: 7.4, y: 2.95, w: 40, h: 0.6, fontSize: 22, color: "C9D6D9", margin: 0, valign: "middle", isTextBox: true,
    objectName: "affiliation" });

// ---------------------------------------------------------------- column A
const AX = 0.6, AW = 13.2;
// 1 Background
card(AX, 4.6, AW, 8.4, "card background");
header(AX + 0.45, 4.95, 1, "Why inversion size?");
s.addText([
  { text: "Inversions flip a segment of the chromosome. They change gene order and orientation and can affect fitness [1]. ", options: {} },
  { text: "They are not placed at random: many sit symmetrically around the origin of replication [2, 3]. ", options: {} },
  { text: "Bigger genomes offer more room, but do they also carry bigger inversions?", options: { bold: true, color: SLATE } },
], { x: AX + 0.5, y: 6.2, w: AW - 1.0, h: 2.3, fontSize: 23, color: "33444B", margin: 0, valign: "top",
  paraSpaceAfter: 6, isTextBox: true, objectName: "background text" });
s.addImage({ path: img("schematic.png"), x: AX + 0.35, y: 7.95, w: 7.3, h: 7.3 * 474 / 704, objectName: "question schematic" });
s.addText([
  { text: "Aim", options: { bold: true, color: TEAL, fontSize: 25, breakLine: true } },
  { text: "Test across 79 bacterial species whether inversion size scales with genome size, and whether it grows more than the extra room alone would predict.", options: {} },
], { x: AX + 7.95, y: 8.3, w: AW - 8.4, h: 4.2, fontSize: 21, color: "33444B", margin: 0, valign: "top",
  paraSpaceAfter: 4, isTextBox: true, objectName: "aim" });


// 2 Approach
const AY2 = 13.5, AH2 = 10.1;
card(AX, AY2, AW, AH2, "card approach");
header(AX + 0.45, AY2 + 0.35, 2, "Approach");
const pw = 10.3;
s.addImage({ path: img("pipeline_snake.png"), x: AX + (AW - pw) / 2, y: AY2 + 1.5, w: pw, h: pw * 931 / 1725,
  objectName: "pipeline" });
s.addImage({ path: img("fig_classes.png"), x: AX + 0.3, y: AY2 + 7.05, w: 3.6, h: 3.6 * 4.6 / 5.6,
  objectName: "spanning vs local diagram" });
s.addText([
  { text: "Spanning: ", options: { bold: true, color: TERRA } },
  { text: "breakpoints on opposite replication arms (crosses ori or ter).", options: { breakLine: true } },
  { text: "Local: ", options: { bold: true, color: TEAL } },
  { text: "both breakpoints on one arm.", options: { breakLine: true } },
  { text: "Room alone: ", options: { bold: true, color: SLATE } },
  { text: "real sizes shuffled across species, same rules, 200×.", options: {} },
], { x: AX + 4.05, y: AY2 + 7.35, w: AW - 4.5, h: 2.6, fontSize: 20, color: "33444B", margin: 0, valign: "top",
  paraSpaceAfter: 5, isTextBox: true, objectName: "definitions" });

// ---------------------------------------------------------------- results
const RX = 14.4, RW = 33.0;
card(RX, 4.6, RW, 9.8, "card results");
header(RX + 0.45, 4.95, 3, "Results", "79 species  ·  3,808 complete genomes  ·  7 phyla  ·  828 inversion families");
const colW = 15.7, c1 = RX + 0.55, c2 = RX + 0.55 + colW + 0.6;
headline(c1, 6.2, colW, "A", "Pooled, inversions grow only as much as the extra room allows", "headline A");
headline(c2, 6.2, colW, "B", "Split by position: spanning inversions grow, local ones do not", "headline B");
const fw = 15.7, fh = 6.5;
s.addImage({ path: img("fig_question.png"), x: c1, y: 6.95, w: fw, h: fh, objectName: "figure A" });
s.addImage({ path: img("fig_split_scatter.png"), x: c2, y: 6.95, w: fw, h: fh, objectName: "figure B" });
caption(c1, 13.55, colW, 0.7, [{ text: "Slope 0.46: ×1.4 per genome doubling, the same as real sizes shuffled across species (dashed, 0.52)." }], "caption A");
caption(c2, 13.55, colW, 0.7, [{ text: "Slopes differ by 0.93 [0.45–1.42], p = 0.0003; both lie outside every simulation of room alone." }], "caption B");

// C boxplots
const CY = 14.9, CH = 8.7, CW = 21.2;
card(RX, CY, CW, CH, "card boxplots");
headline(RX + 0.55, CY + 0.45, CW - 1.1, "C", "Spanning sizes shift up as a whole; local inversions stay small", "headline C");
const bw = CW - 1.0, bh = bw * 6.9 / 20.8;
s.addImage({ path: img("fig_split_boxes.png"), x: RX + 0.5, y: CY + 1.25, w: bw, h: bh, objectName: "figure C" });
caption(RX + 0.55, CY + 1.3 + bh + 0.05, CW - 1.1, 0.5,
  [{ text: "Spanning median and 90th percentile both grow (slopes 0.61, 0.82); removing each species’ largest 10% leaves 0.47 (p = 0.0002)." }],
  "caption C");

// 4 Take-home
const TX = RX + CW + 0.6, TW = RW - CW - 0.6;
card(TX, CY, TW, CH, "card take-home");
header(TX + 0.45, CY + 0.35, 4, "Take-home");
s.addText("Only inversions that cross the origin or terminus grow with the genome.", { x: TX + 0.5, y: CY + 1.5,
  w: TW - 1.0, h: 1.7, fontSize: 32, bold: true, color: TEAL, margin: 0, valign: "top", isTextBox: true,
  objectName: "take-home statement" });
const stat = (x, big, small, color, name) => {
  s.addText(big, { x, y: CY + 3.35, w: 4.9, h: 1.35, fontSize: 74, bold: true, color, margin: 0, valign: "middle",
    isTextBox: true, objectName: `${name} number` });
  s.addText(small, { x, y: CY + 4.7, w: 4.9, h: 1.0, fontSize: 20, color: MUTED, margin: 0, valign: "top",
    isTextBox: true, objectName: `${name} label` });
};
stat(TX + 0.5, "×1.4", "spanning size per genome doubling", TERRA, "stat growth");
stat(TX + 5.75, "3.6×", "symmetric spanning inversions, vs chance", TERRA, "stat symmetry");
s.addText([
  { text: "Local inversions show no increase at any genome size.", options: { bullet: true, breakLine: true } },
  { text: "Robust to terminus placement, relatedness and outliers.", options: { bullet: true, breakLine: true } },
  { text: "Symmetric inversions keep both replication arms balanced [2–5].", options: { bullet: true } },
], { x: TX + 0.5, y: CY + 6.0, w: TW - 0.9, h: 2.5, fontSize: 23, color: "33444B", margin: 0, valign: "top",
  paraSpaceAfter: 10, isTextBox: true, objectName: "take-home bullets" });

// ---------------------------------------------------------------- footer
const FY = 24.05, FH = 2.45;
card(0.6, FY, 46.8, FH, "card footer");
s.addText([
  { text: "References", options: { bold: true, color: SLATE, breakLine: true } },
  { text: "[1] Waters EV et al. (2025) Trends Microbiol 33:875–886", options: { breakLine: true } },
  { text: "[2] Repar J, Warnecke T (2017) Mol Biol Evol 34:1902–1911", options: { breakLine: true } },
  { text: "[3] D’Iorio M, Dewar K (2022) Life Sci Alliance 6:e202201434", options: {} },
], { x: 1.1, y: FY + 0.3, w: 15.5, h: 1.9, fontSize: 19, color: MUTED, margin: 0, valign: "top", paraSpaceAfter: 2,
  isTextBox: true, objectName: "references" });
s.addText([
  { text: " ", options: { breakLine: true } },
  { text: "[4] Bochkareva OO et al. (2018) BMC Genomics 19:965", options: { breakLine: true } },
  { text: "[5] Atre M et al. (2024) Nucleic Acids Res 52:3493", options: { breakLine: true } },
  { text: "Schematic in 1 created with BioRender.com", options: { italic: true } },
], { x: 16.9, y: FY + 0.3, w: 15.5, h: 1.9, fontSize: 19, color: MUTED, margin: 0, valign: "top", paraSpaceAfter: 2,
  isTextBox: true, objectName: "references 2" });
s.addText([
  { text: "Acknowledgements", options: { bold: true, color: SLATE, breakLine: true } },
  { text: "We thank our coach Juliette Luiselli for her guidance. Computing: VSC (Flemish Supercomputer Center), funded by FWO and the Flemish Government.", options: {} },
], { x: 33.0, y: FY + 0.3, w: 10.6, h: 1.9, fontSize: 18, color: MUTED, margin: 0, valign: "top", paraSpaceAfter: 3,
  isTextBox: true, objectName: "acknowledgements" });
s.addImage({ path: img("qr_repo.png"), x: 44.15, y: FY + 0.22, w: 2.0, h: 2.0, objectName: "QR code repository" });
s.addText("Code &\nmethods", { x: 46.2, y: FY + 0.5, w: 1.15, h: 1.4, fontSize: 16, color: MUTED, margin: 0,
  valign: "middle", isTextBox: true, objectName: "QR label" });

(async () => {
  const out = path.join(D, "Poster_IBP2026_inversions.pptx");
  await pres.writeFile({ fileName: out });
  await applyTheme(out, THEME);
  console.log("written", out);
})();
