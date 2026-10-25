const pptxgen = require("pptxgenjs");
const path = require("path");
const { applyTheme } = require("/root/.claude/skills/synced/865e8771-0238-452b-8e45-d6a2b77f2cfe_d7a95306-8d5f-463b-b5b9-558348628557/pptx/scripts/apply_theme.js");

const D = process.env.POSTER_ASSETS || path.join(__dirname, "assets");
const img = (f) => path.join(D, f);
// Blue theme; plot colours are the Okabe-Ito colour-blind safe pair (vermillion D55E00 / blue 0072B2)
const HEAD = "2F5F8F", INK = "1F2D3D", MUTED = "566577", BG = "EEF3F8", CARD = "FFFFFF", LINE = "D6E0EA";
const VERM = "D55E00", BLUE = "0072B2";
const FONT = "Calibri";

const THEME = {
  name: "Inversion poster", headFontFace: FONT, bodyFontFace: FONT,
  colors: { dk1: "000000", lt1: "FFFFFF", dk2: INK, lt2: BG, accent1: HEAD, accent2: VERM, accent3: BLUE,
            accent4: "7D8C99", accent5: "C5CFD8", accent6: LINE, hlink: BLUE, folHlink: HEAD },
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

const shadow = () => ({ type: "outer", blur: 8, offset: 2, angle: 90, color: "1F2D3D", opacity: 0.10 });
const card = (x, y, w, h, name) =>
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.35, fill: { color: CARD },
    line: { color: LINE, width: 1 }, shadow: shadow(), objectName: name });
const header = (x, y, num, title, sub, subRight) => {
  s.addShape(pres.shapes.OVAL, { x, y, w: 0.95, h: 0.95, fill: { color: HEAD }, line: { color: HEAD },
    objectName: `badge ${num}` });
  s.addText(String(num), { x, y, w: 0.95, h: 0.95, align: "center", valign: "middle", fontSize: 30, bold: true,
    color: "FFFFFF", margin: 0, isTextBox: true, objectName: `badge ${num} text` });
  s.addText(title, { x: x + 1.2, y: y - 0.05, w: 12, h: 1.05, fontSize: 38, bold: true, color: INK,
    valign: "middle", margin: 0, isTextBox: true, objectName: `${title} heading` });
  if (sub) s.addText(sub, { x: subRight - 22, y, w: 22, h: 0.95, fontSize: 22, color: MUTED, align: "right",
    valign: "middle", margin: 0, isTextBox: true, objectName: `${title} subtitle` });
};
const headline = (x, y, w, letter, text, name) =>
  s.addText([{ text: letter + "  ", options: { color: HEAD, bold: true } },
             { text, options: { color: INK, bold: true } }],
    { x, y, w, h: 0.95, fontSize: 24, margin: 0, valign: "top", isTextBox: true, objectName: name });
const caption = (x, y, w, text, name) =>
  s.addText(text, { x, y, w, h: 0.75, fontSize: 17, color: MUTED, margin: 0, valign: "top", isTextBox: true,
    objectName: name });

// ---------------------------------------------------------------- header
s.addShape(pres.shapes.RECTANGLE, { x: 0, y: 0, w: 48, h: 4.1, fill: { color: HEAD }, line: { color: HEAD },
  objectName: "title band" });
s.addImage({ path: img("logo.png"), x: 0.7, y: 0.95, w: 6.0, h: 6.0 * 157 / 438, objectName: "KU Leuven logo" });
s.addText("Do bigger bacterial genomes carry bigger inversions?", { x: 7.4, y: 0.35, w: 40, h: 1.7, fontSize: 84,
  bold: true, color: "FFFFFF", margin: 0, valign: "middle", isTextBox: true, objectName: "title" });
s.addText("Leen Al Radi  ·  Mathias Azevedo  ·  Abishek Kumar  ·  Wout Simoens", { x: 7.4, y: 2.15, w: 40, h: 0.75,
  fontSize: 32, color: "FFFFFF", margin: 0, valign: "middle", isTextBox: true, objectName: "authors" });
s.addText("Master of Bioinformatics, Faculty of Bioscience Engineering, KU Leuven   |   Coach: Juliette Luiselli   |   Integrated Bioinformatics Project 2026",
  { x: 7.4, y: 2.95, w: 40, h: 0.6, fontSize: 22, color: "DCE7F2", margin: 0, valign: "middle", isTextBox: true,
    objectName: "affiliation" });

// ---------------------------------------------------------------- row 1: background + approach
const R1 = 4.5, R1H = 8.0;
const AX = 0.6, AW = 13.2;
card(AX, R1, AW, R1H, "card background");
header(AX + 0.45, R1 + 0.35, 1, "Why inversion size?");
s.addText([
  { text: "Inversions flip a segment of the chromosome. They change gene order and orientation and can affect fitness [1]. ", options: {} },
  { text: "They are not placed at random: many sit symmetrically around the origin of replication [2, 3]. ", options: {} },
  { text: "Bigger genomes offer more room, but do they also carry bigger inversions?", options: { bold: true, color: INK } },
], { x: AX + 0.5, y: R1 + 1.55, w: AW - 1.0, h: 2.1, fontSize: 22, color: "33424F", margin: 0, valign: "top",
  isTextBox: true, objectName: "background text" });
s.addImage({ path: img("schematic.png"), x: AX + 0.35, y: R1 + 3.75, w: 6.0, h: 6.0 * 474 / 704,
  objectName: "question schematic" });
s.addText([
  { text: "Aim", options: { bold: true, color: HEAD, fontSize: 25, breakLine: true } },
  { text: "Test across 79 bacterial species whether inversion size scales with genome size, and whether it grows more than the extra room alone would predict.", options: {} },
], { x: AX + 6.75, y: R1 + 4.0, w: AW - 7.2, h: 3.7, fontSize: 21, color: "33424F", margin: 0, valign: "top",
  paraSpaceAfter: 4, isTextBox: true, objectName: "aim" });

const PX = AX + AW + 0.6, PW = 48 - 0.6 - PX;
card(PX, R1, PW, R1H, "card approach");
header(PX + 0.45, R1 + 0.35, 2, "Approach");
const pw = 25.8;
s.addImage({ path: img("pipeline_row.png"), x: PX + 0.45, y: R1 + 1.65, w: pw, h: pw * 784 / 3549,
  objectName: "pipeline" });
const KX = PX + 0.45 + pw + 0.5, KW = PX + PW - 0.4 - KX;
s.addImage({ path: img("fig_classes.png"), x: KX + (KW - 4.4) / 2, y: R1 + 0.3, w: 4.4, h: 4.4 * 4.6 / 5.6,
  objectName: "spanning vs local diagram" });
s.addText([
  { text: "Spanning: ", options: { bold: true, color: VERM } },
  { text: "breakpoints on opposite replication arms (crosses ori or ter).", options: { breakLine: true } },
  { text: "Local: ", options: { bold: true, color: BLUE } },
  { text: "both breakpoints on one arm.", options: { breakLine: true } },
  { text: "Room alone: ", options: { bold: true, color: INK } },
  { text: "real sizes shuffled across species, same rules, 200×.", options: {} },
], { x: KX, y: R1 + 4.0, w: KW, h: 3.8, fontSize: 19, color: "33424F", margin: 0, valign: "top",
  paraSpaceAfter: 6, isTextBox: true, objectName: "definitions" });

// ---------------------------------------------------------------- row 2: results + take-home
const R2 = 12.9, R2H = 10.8;
const RX = 0.6, RW = 34.8;
card(RX, R2, RW, R2H, "card results");
header(RX + 0.45, R2 + 0.35, 3, "Results", "79 species  ·  3,808 complete genomes  ·  7 phyla  ·  828 inversion families",
  RX + RW - 0.5);
const fh = 7.35;
const A = { x: RX + 0.5, w: fh * 11.8 / 7.6 }, B = { w: fh * 10.4 / 7.6 }, C = { w: fh * 10.4 / 7.6 };
const gap = (RW - 1.0 - A.w - B.w - C.w) / 2;
B.x = A.x + A.w + gap; C.x = B.x + B.w + gap;
const HY = R2 + 1.45, FY0 = R2 + 2.45, CYp = FY0 + fh + 0.1;
headline(A.x, HY, A.w, "A", "Pooled, inversions grow only as much as the extra room allows", "headline A");
headline(B.x, HY, B.w, "B", "Split by position: spanning grows, local does not", "headline B");
headline(C.x, HY, C.w, "C", "Spanning sizes shift as a whole; local stay small", "headline C");
s.addImage({ path: img("fig_question.png"), x: A.x, y: FY0, w: A.w, h: fh, objectName: "figure A" });
s.addImage({ path: img("fig_split_scatter.png"), x: B.x, y: FY0, w: B.w, h: fh, objectName: "figure B" });
s.addImage({ path: img("fig_split_boxes.png"), x: C.x, y: FY0, w: C.w, h: fh, objectName: "figure C" });
caption(A.x, CYp, A.w, "Slope 0.46 (×1.4 per genome doubling) matches real sizes shuffled across species (0.52).", "caption A");
caption(B.x, CYp, B.w, "Slopes differ by 0.93 [0.45–1.42], p = 0.0003; both lie outside all 200 room-alone simulations.", "caption B");
caption(C.x, CYp, C.w, "Spanning median and 90th percentile both grow (0.61, 0.82); without each species’ top 10%: 0.47.", "caption C");

const TX = RX + RW + 0.6, TW = 48 - 0.6 - TX;
card(TX, R2, TW, R2H, "card take-home");
header(TX + 0.45, R2 + 0.35, 4, "Take-home");
s.addText("Only inversions that cross the origin or terminus grow with the genome.", { x: TX + 0.5, y: R2 + 1.5,
  w: TW - 1.0, h: 1.4, fontSize: 30, bold: true, color: HEAD, margin: 0, valign: "top", isTextBox: true,
  objectName: "take-home statement" });
const stat = (x, big, small, name) => {
  s.addText(big, { x, y: R2 + 2.95, w: 5.0, h: 1.15, fontSize: 66, bold: true, color: VERM, margin: 0,
    valign: "middle", isTextBox: true, objectName: `${name} number` });
  s.addText(small, { x, y: R2 + 4.1, w: 5.0, h: 0.8, fontSize: 20, color: MUTED, margin: 0, valign: "top",
    isTextBox: true, objectName: `${name} label` });
};
stat(TX + 0.5, "×1.4", "spanning size per genome doubling", "stat growth");
stat(TX + 5.8, "3.6×", "symmetric spanning inversions, vs chance", "stat symmetry");
s.addText([
  { text: "Pooled, inversion size rises only as fast as the extra room in bigger genomes (0.46 vs 0.52).", options: { bullet: true, breakLine: true } },
  { text: "Spanning inversions grow ×1.4 per genome doubling, and their whole size distribution shifts up.", options: { bullet: true, breakLine: true } },
  { text: "Local inversions stay around 60 kb at any genome size.", options: { bullet: true, breakLine: true } },
  { text: "The growing inversions sit symmetrically around the origin or terminus, keeping both replication arms balanced [2–5].", options: { bullet: true, breakLine: true } },
  { text: "Holds when the terminus is moved, related species are removed or outliers are dropped.", options: { bullet: true } },
], { x: TX + 0.5, y: R2 + 5.15, w: TW - 0.9, h: 5.4, fontSize: 23, color: "33424F", margin: 0, valign: "top",
  paraSpaceAfter: 12, isTextBox: true, objectName: "take-home bullets" });

// ---------------------------------------------------------------- footer
const FY = 24.1, FH = 2.4;
card(0.6, FY, 46.8, FH, "card footer");
s.addText([
  { text: "References", options: { bold: true, color: INK, breakLine: true } },
  { text: "[1] Waters EV et al. (2025) Trends Microbiol 33:875–886", options: { breakLine: true } },
  { text: "[2] Repar J, Warnecke T (2017) Mol Biol Evol 34:1902–1911", options: { breakLine: true } },
  { text: "[3] D’Iorio M, Dewar K (2022) Life Sci Alliance 6:e202201434", options: {} },
], { x: 1.1, y: FY + 0.28, w: 15.5, h: 1.9, fontSize: 19, color: MUTED, margin: 0, valign: "top", paraSpaceAfter: 2,
  isTextBox: true, objectName: "references" });
s.addText([
  { text: " ", options: { breakLine: true } },
  { text: "[4] Bochkareva OO et al. (2018) BMC Genomics 19:965", options: { breakLine: true } },
  { text: "[5] Atre M et al. (2024) Nucleic Acids Res 52:3493", options: { breakLine: true } },
  { text: "Schematic in 1 created with BioRender.com", options: { italic: true } },
], { x: 16.9, y: FY + 0.28, w: 15.5, h: 1.9, fontSize: 19, color: MUTED, margin: 0, valign: "top", paraSpaceAfter: 2,
  isTextBox: true, objectName: "references 2" });
s.addText([
  { text: "Acknowledgements", options: { bold: true, color: INK, breakLine: true } },
  { text: "We thank our coach Juliette Luiselli for her guidance. Computing: VSC (Flemish Supercomputer Center), funded by FWO and the Flemish Government.", options: {} },
], { x: 33.0, y: FY + 0.28, w: 10.6, h: 1.9, fontSize: 18, color: MUTED, margin: 0, valign: "top", paraSpaceAfter: 3,
  isTextBox: true, objectName: "acknowledgements" });
s.addImage({ path: img("qr_repo.png"), x: 44.15, y: FY + 0.2, w: 2.0, h: 2.0, objectName: "QR code repository" });
s.addText("Code &\nmethods", { x: 46.2, y: FY + 0.5, w: 1.15, h: 1.4, fontSize: 16, color: MUTED, margin: 0,
  valign: "middle", isTextBox: true, objectName: "QR label" });

(async () => {
  const out = path.join(process.env.POSTER_OUT || D, "Poster_IBP2026_inversions.pptx");
  await pres.writeFile({ fileName: out });
  await applyTheme(out, THEME);
  console.log("written", out);
})();
