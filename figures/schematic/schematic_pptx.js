// Editable PowerPoint version of figure fig:schematic (arc_schematic.py drawn as native shapes).
//     node schematic_pptx.js            -> schematic.pptx (OVERWRITES it: edit the .pptx by hand afterwards)
// One 13.33 x 7.5 in slide; every shape is named (Home > Arrange > Selection Pane).
// (a) the arc between the W-ThO2 cathode and the Cu anode, with the two resolved regions zoomed into
// (b) the cathode region (cathode on the left) and (c) the anode region (anode on the right), not to scale.
const path = require("path");
const pptxgen = require("pptxgenjs");

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.title = "Arc schematic";
const slide = pres.addSlide();
slide.background = { color: "FFFFFF" };
const S = pres.shapes;

const FONT = "Times New Roman";
const C_E = "2A78D6", C_I = "EB6834", INK = "2B2B2B", LOSS = "D13B3A", HEAT = "8B1A1A";
const W_COL = "8C8F94", CU_COL = "B87333", PLASMA = "F6B26B";
const SHEATH = "CFE0F5", LAYER = "FDE3C4", LTE = "FFF4E6", GREY = "9A9A9A";

// ---------- helpers ----------
// Rich text: $..$ italic, _{..} subscript, ^{..} superscript, _[..] italic subscript, \n new line.
function runs(src, base) {
  const out = [];
  const lines = src.split("\n");
  lines.forEach((line, li) => {
    const parts = line.split(/(\$[^$]+\$|_\{[^}]*\}|_\[[^\]]*\]|\^\{[^}]*\})/).filter((p) => p !== "");
    parts.forEach((p, pi) => {
      const o = Object.assign({}, base);
      let t = p;
      if (p.startsWith("$")) { t = p.slice(1, -1); o.italic = true; }
      else if (p.startsWith("_{")) { t = p.slice(2, -1); o.subscript = true; }
      else if (p.startsWith("_[")) { t = p.slice(2, -1); o.subscript = true; o.italic = true; }
      else if (p.startsWith("^{")) { t = p.slice(2, -1); o.superscript = true; }
      if (pi === parts.length - 1 && li < lines.length - 1) o.breakLine = true;
      out.push({ text: t, options: o });
    });
  });
  return out;
}
function text(name, src, x, y, w, h, opt = {}) {
  const base = { fontFace: FONT, fontSize: opt.size || 15, color: opt.color || INK };
  if (opt.bold) base.bold = true;
  slide.addText(runs(src, base), {
    x, y, w, h, isTextBox: true, margin: 0, objectName: name,
    align: opt.align || "center", valign: opt.valign || "middle", lineSpacingMultiple: 0.95,
  });
}
function arrow(name, x1, y1, x2, y2, color, width = 1.5, dash) {
  const line = { color, width, endArrowType: "triangle" };
  if (dash) line.dashType = dash;
  slide.addShape(S.LINE, {
    x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.abs(x2 - x1), h: Math.abs(y2 - y1),
    flipH: x2 < x1, flipV: y2 < y1, line, objectName: name,
  });
}
function seg(name, x1, y1, x2, y2, color, width = 0.75, dash) {
  const line = { color, width };
  if (dash) line.dashType = dash;
  slide.addShape(S.LINE, {
    x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.abs(x2 - x1) || 0.001, h: Math.abs(y2 - y1) || 0.001,
    flipH: (x2 < x1) !== (y2 < y1) && x1 !== x2 && y1 !== y2, line, objectName: name,
  });
}
// Free-form shape through absolute points (inches); closed shapes get a fill, open ones a line only.
function poly(name, pts, { fill, transparency = 0, line, lineWidth = 1, closed = true } = {}) {
  const xs = pts.map((p) => p[0]), ys = pts.map((p) => p[1]);
  const x0 = Math.min(...xs), y0 = Math.min(...ys);
  const w = Math.max(Math.max(...xs) - x0, 0.01), h = Math.max(Math.max(...ys) - y0, 0.01);
  const points = pts.map((p, i) => (i === 0 ? { x: p[0] - x0, y: p[1] - y0, moveTo: true } : { x: p[0] - x0, y: p[1] - y0 }));
  if (closed) points.push({ close: true });
  const o = { x: x0, y: y0, w, h, points, objectName: name };
  o.fill = fill ? { color: fill, transparency } : { type: "none" };
  o.line = line ? { color: line, width: lineWidth } : { type: "none" };
  slide.addShape(S.CUSTOM_GEOMETRY, o);
}
function rect(name, x, y, w, h, fill, line, lineWidth = 1, dash) {
  const o = { x, y, w, h, objectName: name };
  o.fill = fill ? { color: fill } : { type: "none" };
  o.line = line ? { color: line, width: lineWidth, dashType: dash || "solid" } : { type: "none" };
  slide.addShape(S.RECTANGLE, o);
}
const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));

// ---------- (a) the arc ----------
// Geometry of arc_schematic.py (units: its data coordinates), mapped to inches.
const A_SCALE = 0.62, A_CX = 1.95, A_TOP = 0.55;
const ax = (u) => A_CX + A_SCALE * u, ay = (y) => A_TOP + A_SCALE * (10.3 - y);
const YC = 7.3, YA = 1.0;                       // cathode tip, anode surface
text("label (a)", "(a)", 0.25, 0.2, 0.5, 0.35, { size: 17, align: "left" });

const ys = linspace(YA, YC, 120);
const rad = ys.map((y) => { const s = (y - YA) / (YC - YA); return 0.11 + 0.8 * Math.sqrt(1 - s) + 0.25 * Math.pow(1 - s, 6); });
[[1.25, 82, "column glow outer"], [1.1, 70, "column glow"], [1.0, 45, "arc column"], [0.45, 10, "column core"]].forEach(([k, tr, nm]) => {
  const right = ys.map((y, i) => [ax(k * rad[i]), ay(y)]);
  const left = ys.slice().reverse().map((y, i) => [ax(-k * rad[ys.length - 1 - i]), ay(y)]);
  poly(nm, right.concat(left), { fill: nm === "column core" ? "FBE7C6" : PLASMA, transparency: tr });
});
poly("cathode", [[-0.32, 10.2], [0.32, 10.2], [0.32, 8.4], [0.06, YC], [-0.06, YC], [-0.32, 8.4]].map(([u, y]) => [ax(u), ay(y)]),
  { fill: W_COL, line: INK, lineWidth: 0.75 });
text("cathode label", "cathode\n(W–ThO_{2})", ax(0.45), ay(10.0), 1.3, 0.6, { align: "left" });
rect("anode", ax(-1.55), ay(YA), A_SCALE * 3.1, A_SCALE * 0.75, CU_COL, INK, 0.75);
text("anode label", "anode (Cu)", ax(-1.55), ay(YA), A_SCALE * 3.1, A_SCALE * 0.75, { color: "FFFFFF" });
text("column label", "LTE arc column\n(bulk model)", ax(-1.0), ay(4.6), A_SCALE * 2.0, 0.6, { size: 15 });
arrow("current arrow", ax(1.45), ay(3.0), ax(1.45), ay(6.0), INK, 1.25);
text("current label", "$I$", ax(1.52), ay(4.75), 0.25, 0.3, { size: 16, align: "left" });

// resolved regions (thickness exaggerated) and the 1-D x axes
const HC = 0.55;
const cBox = [ax(-0.32), ay(YC), A_SCALE * 0.64, A_SCALE * HC];
const aBox = [ax(-1.25), ay(YA + HC), A_SCALE * 2.5, A_SCALE * HC];
rect("cathode region box", ...cBox, null, C_E, 1.25, "dash");
rect("anode region box", ...aBox, null, C_E, 1.25, "dash");
arrow("x axis cathode", ax(0), ay(YC), ax(0), ay(YC - HC - 0.4), INK, 1);
text("x cathode", "$x$", ax(0.07), ay(YC - HC - 0.25), 0.2, 0.25, { align: "left" });
arrow("x axis anode", ax(0.95), ay(YA), ax(0.95), ay(YA + HC + 0.45), INK, 1);
text("x anode", "$x$", ax(1.02), ay(YA + HC + 0.55), 0.2, 0.25, { align: "left" });
text("cathode region label", "cathode\nregion (b)", ax(-0.42) - 1.1, ay(YC) - 0.05, 1.1, 0.5, { color: C_E, align: "right", size: 14 });
text("anode region label", "anode region (c)", ax(-1.25), ay(YA + HC) - 0.32, 1.6, 0.28, { color: C_E, align: "left", size: 14 });

// ---------- (b), (c) the resolved regions ----------
const X0 = 4.75, X1 = 12.75, WD = X1 - X0;     // panel extent (x = 0 .. L), electrode outside
const US = 0.11, UL = 0.6;                     // sheath edge, end of the non-equilibrium layer
const H = 2.4;

function panel(el, top) {
  const cat = el === "cathode";
  const X = (u) => (cat ? X0 + u * WD : X1 - u * WD);       // u = 0 at the electrode, 1 at x = L
  const Y = (v) => top + H * (1 - v);                       // v = 0 bottom, 1 top of the strip
  const P = cat ? "(b) " : "(c) ";
  const band = (name, u0, u1, col) => rect(P + name, Math.min(X(u0), X(u1)), top, Math.abs(X(u1) - X(u0)), H, col);
  band("sheath", 0, US, SHEATH);
  band("non-equilibrium layer", US, UL, LAYER);
  band("near-LTE plasma", UL, 1, LTE);
  seg(P + "sheath edge", X(US), top, X(US), top + H, INK, 0.75, "sysDot");
  seg(P + "layer edge", X(UL), top, X(UL), top + H, INK, 0.75, "sysDot");
  seg(P + "x = L boundary", X(1), top - 0.05, X(1), top + H + 0.05, INK, 1.5, "dash");
  const ex = cat ? X0 - 0.32 : X1;
  rect(P + (cat ? "cathode" : "anode"), ex, top - 0.06, 0.32, H + 0.12, cat ? W_COL : CU_COL, INK, 0.75);

  text(P + "title", (cat ? "(b) cathode region" : "(c) anode region") + " (not to scale)",
    cat ? X0 - 0.32 : X1 - 5.5 + 0.32, top - 0.45, 5.5, 0.35, { size: 17, align: cat ? "left" : "right" });

  // zone labels
  const lab = (name, u0, u1, s, v, size, color) => text(P + name, s, Math.min(X(u0), X(u1)), Y(v), Math.abs(X(u1) - X(u0)), 0.27, { size, color });
  lab("sheath label", 0, US, "sheath", 0.985, 15);
  lab("sheath width", 0, US, "≲ 1 μm", 0.875, 13);
  lab("sheath charge", 0, US, "$n$_[i] > $n$_[e]", 0.775, 13);
  if (cat) {
    lab("layer label", US, UL, "ionization layer, ~10^{−4} m", 0.985, 15);
    lab("layer detail", US, UL, "$T$_[h] < $T$_[e],  $n$_[e] ≫ $n$_{eq}($T$_[h])", 0.875, 13);
  } else {
    lab("layer label", US, UL, "quasi-neutral anode layer", 0.985, 15);
    lab("layer detail", US, UL, "$T$_[h] < $T$_[e],  $n$_[e] > $n$_{eq}($T$_[h])  (recombination)", 0.875, 13);
  }
  lab("LTE label", UL, 1, "near-LTE plasma", 0.985, 15);
  lab("LTE detail", UL, 1, "$T$_[e] ≈ $T$_[h],  $n$_[e] ≈ $n$_{eq}", 0.875, 13);

  // qualitative temperatures (values of arc_schematic.py, rescaled to the strip)
  const sc = (t) => 0.39 + ((t - 0.24) / (0.56 - 0.24)) * 0.24;
  const us = linspace(0, 1, 160);
  let Te, Th;
  if (cat) {
    Te = us.map((u) => 0.53 + 0.03 / (1 + Math.exp(-(u - 0.6) / 0.08)) - 0.012 * Math.exp(-(((u - 0.3) / 0.06) ** 2)));
    Th = us.map((u) => 0.3 + 0.26 * (1 - Math.exp(-u / 0.17)));
  } else {
    Te = us.map((u) => 0.4 + (0.16 * (1 - Math.exp(-(u + 0.02) / 0.4))) / (1 - Math.exp(-1.02 / 0.4)));
    Th = us.map((u) => 0.24 + (0.32 * (1 - Math.exp(-u / 0.28))) / (1 - Math.exp(-1 / 0.28)));
  }
  poly(P + "T_e curve", us.map((u, i) => [X(u), Y(sc(Te[i]))]), { line: C_E, lineWidth: 2, closed: false });
  poly(P + "T_h curve", us.map((u, i) => [X(u), Y(sc(Th[i]))]), { line: C_I, lineWidth: 2, closed: false });
  const k = cat ? 24 : 75;
  text(P + "T_e label", "$T$_[e]", X(us[k]) - 0.2, Y(sc(Te[k])) - 0.32, 0.4, 0.28, { color: C_E, size: 15 });
  text(P + "T_h label", "$T$_[h]", X(us[k]) - 0.2, Y(sc(Th[k])) + 0.05, 0.4, 0.28, { color: C_I, size: 15 });

  // fluxes at the surface (lower half of the sheath); labels just outside the sheath
  const fl = cat
    ? [[HEAT, "heat flux $q$_[w]", true, 2.5], [C_I, "ions", true, 1.5], [C_E, "plasma e^{−}", true, 1.5], [C_E, "emitted e^{−}", false, 1.5]]
    : [[HEAT, "heat flux $q$_[w]", true, 2.5], [C_E, "e^{−}", true, 1.5], [C_I, "ions", true, 1.5]];
  fl.forEach(([col, s, inward, wd], i) => {
    const v = 0.31 - i * 0.08;
    const a = X(0.012), b = X(US - 0.008);
    if (inward) arrow(P + "flux " + i, b, Y(v), a, Y(v), col, wd); else arrow(P + "flux " + i, a, Y(v), b, Y(v), col, wd);
    const tx = cat ? X(US) + 0.06 : X(US) - 1.66;
    text(P + "flux label " + i, s, tx, Y(v) - 0.12, 1.6, 0.24, { color: col, size: 13, align: cat ? "left" : "right" });
  });

  // radial loss, the whole layer and the plasma (out of the column)
  [0.36, 0.5, 0.64, 0.78, 0.92].forEach((u, i) => arrow(P + "radial loss " + i, X(u), Y(0.25), X(u), Y(0.05), LOSS, 1.25));
  text(P + "radial loss label", "radial loss $Q$_[e], $Q$_[h]", X(0.64) - 1.1, Y(0.36), 2.2, 0.26,
    { color: LOSS, size: 13, align: "center" });

  // boundary data
  const yb = top + H + 0.1;
  const surf = cat
    ? "$x$ = 0:  $T$_{wall} from the heat balance,\nthermionic emission, wall fluxes"
    : "$x$ = 0:  $T$_{wall} prescribed,\nwall fluxes";
  const bulk = "$x$ = $L$:  $T$_[e] = $T$_[h] = $T$_{bulk},  $n$_[e] = $n$_[i] = $n$_{eq}($T$_{bulk}),\n$J$_{app} (from the LTE bulk model)";
  text(P + "surface data", surf, cat ? X0 : X1 - 3.4, yb, 3.4, 0.5, { size: 13, align: cat ? "left" : "right", valign: "top" });
  text(P + "bulk data", bulk, cat ? X1 - 4.6 : X0, yb, 4.6, 0.5, { size: 13, align: cat ? "right" : "left", valign: "top" });
  return { top, X };
}

const pb = panel("cathode", 0.75);
const pc = panel("anode", 4.45);

// zoom lines from the dashed boxes in (a) to the panels
seg("zoom cathode top", cBox[0] + cBox[2], cBox[1], X0 - 0.32, pb.top - 0.06, GREY, 0.75, "dash");
seg("zoom cathode bottom", cBox[0] + cBox[2], cBox[1] + cBox[3], X0 - 0.32, pb.top + H + 0.06, GREY, 0.75, "dash");
seg("zoom anode top", aBox[0] + aBox[2], aBox[1], X0, pc.top, GREY, 0.75, "dash");
seg("zoom anode bottom", aBox[0] + aBox[2], aBox[1] + aBox[3], X0, pc.top + H, GREY, 0.75, "dash");

pres.writeFile({ fileName: path.join(__dirname, "schematic.pptx") }).then((f) => console.log("wrote " + f));
