'use strict';

const PARK_W = 3;
const PARK_H = 2;
const SNAP_DIST = 0.15;
const PALETTE = ['#c74440', '#2d70b3', '#388c46', '#6042a6', '#fa7e19', '#0e8ac8', '#b8860b', '#c2185b'];

const state = {
  segments: [],
  features: [],
  mode: 'draw',
  view: { cx: 1.5, cy: 1.0, ppk: 260 },
  selectedFeature: null,
  selectedSeg: null,
  isolateId: null,
  nextId: 1
};

function round3(v) { return Math.round(v * 1000) / 1000; }

function fmt(v, digits) {
  const d = digits === undefined ? 3 : digits;
  if (!isFinite(v)) return '0';
  let s = v.toFixed(d);
  if (s.indexOf('.') >= 0) s = s.replace(/0+$/, '').replace(/\.$/, '');
  if (s === '-0') s = '0';
  return s;
}

function zeros(n) { return new Array(n).fill(0); }

function solveLinear(M, v) {
  const n = v.length;
  const A = M.map((row, i) => row.concat([v[i]]));
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
    if (Math.abs(A[p][c]) < 1e-12) return null;
    const t = A[c]; A[c] = A[p]; A[p] = t;
    for (let r = 0; r < n; r++) {
      if (r === c) continue;
      const f = A[r][c] / A[c][c];
      if (f === 0) continue;
      for (let k = c; k <= n; k++) A[r][k] -= f * A[c][k];
    }
  }
  return A.map((row, i) => row[n] / A[i][i]);
}

function lstsq(X, ys) {
  const m = X[0].length;
  const M = [];
  for (let i = 0; i < m; i++) M.push(zeros(m));
  const v = zeros(m);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < m; j++) {
      let s = 0;
      for (let r = 0; r < X.length; r++) s += X[r][i] * X[r][j];
      M[i][j] = s;
    }
    let s = 0;
    for (let r = 0; r < X.length; r++) s += X[r][i] * ys[r];
    v[i] = s;
  }
  return solveLinear(M, v);
}

function bestConstantRms(base, xs, ys) {
  const n = xs.length;
  let c = 0;
  for (let i = 0; i < n; i++) c += ys[i] - base(xs[i]);
  c /= n;
  let s = 0;
  for (let i = 0; i < n; i++) { const e = base(xs[i]) + c - ys[i]; s += e * e; }
  return { c, rms: Math.sqrt(s / n) };
}

function resamplePts(raw, n) {
  const pts = raw.map(p => ({ x: p.x, y: p.y }));
  pts.sort((a, b) => a.x - b.x);
  const x1 = pts[0].x, x2 = pts[pts.length - 1].x;
  const out = [];
  const w = (x2 - x1) / (n * 2) + 1e-12;
  for (let i = 0; i < n; i++) {
    const xt = x1 + (x2 - x1) * i / (n - 1);
    let sum = 0, cnt = 0, best = pts[0].y, bd = Infinity;
    for (const p of pts) {
      const d = Math.abs(p.x - xt);
      if (d <= w) { sum += p.y; cnt++; }
      if (d < bd) { bd = d; best = p.y; }
    }
    out.push({ x: xt, y: cnt ? sum / cnt : best });
  }
  return out;
}

function candidateFitters(xs, ys, x1, x2) {
  const spanY = (Math.max(...ys) - Math.min(...ys)) || 0.1;
  const xr = Math.max(x2 - x1, 0.1);
  const cands = [];

  function add(type, label, params, base, penalty) {
    const { rms } = bestConstantRms(base, xs, ys);
    cands.push({ type, label, params, base, rms, score: rms / spanY + penalty });
  }

  for (const deg of [1, 2, 3]) {
    const X = xs.map(x => {
      const row = [];
      for (let d = deg; d >= 1; d--) row.push(Math.pow(x, d));
      row.push(1);
      return row;
    });
    const co = lstsq(X, ys);
    if (!co) continue;
    const cs = co.map(round3);
    const base = x => {
      let s = 0;
      for (let d = 1; d <= deg; d++) s += cs[deg - d] * Math.pow(x, d);
      return s + cs[deg];
    };
    const label = deg === 1 ? 'Linear (polynomial)' : deg === 2 ? 'Quadratic (polynomial)' : 'Cubic (polynomial)';
    add('poly' + deg, label, { cs }, base, deg === 1 ? 0.004 : deg === 2 ? 0.006 : 0.008);
  }

  for (const f of [0, 0.05, 0.1, 0.2]) {
    const h = round3(x1 - f * xr);
    const X = xs.map(x => [Math.sqrt(Math.max(x - h, 0)), 1]);
    const co = lstsq(X, ys);
    if (!co) continue;
    const a = round3(co[0]);
    if (!isFinite(a) || Math.abs(a) < 1e-4) continue;
    const base = x => a * Math.sqrt(Math.max(x - h, 0));
    add('sqrt', 'Square root', { a, h }, base, 0.006);
  }

  for (const f of [0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3]) {
    const h = round3(x1 - f * xr);
    if (x1 - h < 0.02) continue;
    const X = xs.map(x => [1 / (x - h), 1]);
    const co = lstsq(X, ys);
    if (!co) continue;
    const a = round3(co[0]);
    if (!isFinite(a) || Math.abs(a) < 1e-4) continue;
    const base = x => a / (x - h);
    add('hyperbola', 'Rectangular hyperbola', { a, h }, base, 0.006);
  }

  for (const f of [0.02, 0.05, 0.1, 0.2, 0.4]) {
    const h = round3(x1 - f * xr);
    if (x1 - h < 0.02) continue;
    const X = xs.map(x => [Math.log(x - h), 1]);
    const co = lstsq(X, ys);
    if (!co) continue;
    const a = round3(co[0]);
    if (!isFinite(a) || Math.abs(a) < 1e-4) continue;
    const base = x => a * Math.log(x - h);
    add('log', 'Logarithmic', { a, h }, base, 0.008);
  }

  const ymin = Math.min(...ys), ymax = Math.max(...ys);
  const cList = [];
  for (const f of [0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32]) {
    cList.push(ymin - spanY * f);
    cList.push(ymax + spanY * f);
  }
  let expFallback = null;
  const before = cands.length;
  for (const c0 of cList) {
    const d = ys.map(y => y - c0);
    const sign = d[0] > 0 ? 1 : -1;
    if (d.some(v => v * sign <= 0)) continue;
    const X = xs.map(x => [x - x1, 1]);
    const u = d.map(v => Math.log(Math.abs(v)));
    const co = lstsq(X, u);
    if (!co) continue;
    const p = co[0], q = co[1];
    const b = round3(Math.exp(p));
    const a = round3(sign * Math.exp(q));
    if (!isFinite(a) || !isFinite(b) || b <= 0.01 || b >= 50 || Math.abs(b - 1) < 0.02) continue;
    if (Math.abs(a) < 1e-4 || Math.abs(a) > 1e6) continue;
    const base = x => a * Math.pow(b, x - x1);
    if (a > 0) add('exp', 'Exponential', { a, b, h: x1 }, base, 0.008);
    else {
      const { rms } = bestConstantRms(base, xs, ys);
      if (!expFallback || rms < expFallback.rms) expFallback = { type: 'exp', label: 'Exponential', params: { a, b, h: x1 }, base, rms, score: rms / spanY + 0.008 };
    }
  }
  if (cands.length === before && expFallback) cands.push(expFallback);

  if (spanY > 1e-6) {
    const rmsAt = b => {
      const X = xs.map(x => [Math.sin(b * x), Math.cos(b * x), 1]);
      const co = lstsq(X, ys);
      if (!co) return Infinity;
      let s = 0;
      for (let i = 0; i < xs.length; i++) {
        const e = co[0] * Math.sin(b * xs[i]) + co[1] * Math.cos(b * xs[i]) + co[2] - ys[i];
        s += e * e;
      }
      return Math.sqrt(s / xs.length);
    };
    let best = null;
    for (let b = 0.3; b <= 12.0001; b += 0.02) {
      const rms = rmsAt(b);
      if (!best || rms < best.rms) best = { b, rms };
    }
    if (best) {
      let lo = Math.max(0.05, best.b - 0.06), hi = best.b + 0.06;
      for (let it = 0; it < 50; it++) {
        const m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
        if (rmsAt(m1) < rmsAt(m2)) hi = m2; else lo = m1;
      }
      const bBest = (lo + hi) / 2;
      const X = xs.map(x => [Math.sin(bBest * x), Math.cos(bBest * x), 1]);
      const co = lstsq(X, ys);
      const aRaw = Math.hypot(co[0], co[1]);
      if (aRaw >= 0.35 * spanY && aRaw > 1e-4) {
        const phi = Math.atan2(co[1], co[0]);
        const b = round3(bBest);
        const a = round3(aRaw);
        const h = round3(-phi / bBest);
        const base = x => a * Math.sin(b * (x - h));
        add('sine', 'Sine', { a, b, h }, base, 0.01);
      }
    }
  }

  return cands;
}

function coefLatex(a) {
  if (Math.abs(a - 1) < 1e-9) return '';
  if (Math.abs(a + 1) < 1e-9) return '-';
  return fmt(a);
}

function xShiftLatex(h) {
  if (Math.abs(h) < 1e-9) return 'x';
  return h > 0 ? 'x-' + fmt(h) : 'x+' + fmt(-h);
}

function plusLatex(v) {
  if (Math.abs(v) < 1e-9) return '';
  return v > 0 ? '+' + fmt(v) : '-' + fmt(-v);
}

function desmosString(seg) {
  const D = '\\left\\{' + fmt(seg.x1) + '\\le x\\le' + fmt(seg.x2) + '\\right\\}';
  const p = seg.params;
  let body;
  if (seg.type === 'poly1' || seg.type === 'poly2' || seg.type === 'poly3') {
    const cs = p.cs;
    const deg = cs.length - 1;
    const parts = [];
    for (let d = deg; d >= 0; d--) {
      const a = cs[deg - d];
      if (Math.abs(a) < 1e-9) continue;
      const suffix = d === 0 ? '' : d === 1 ? 'x' : 'x^{' + d + '}';
      const coef = (Math.abs(Math.abs(a) - 1) < 1e-9 && d > 0) ? '' : fmt(Math.abs(a));
      parts.push({ neg: a < 0, coef, suffix });
    }
    if (!parts.length) body = '0';
    else body = parts.map((t, i) => (i === 0 ? (t.neg ? '-' : '') : (t.neg ? '-' : '+')) + t.coef + t.suffix).join('');
  } else if (seg.type === 'sqrt') {
    body = coefLatex(p.a) + '\\sqrt{' + xShiftLatex(p.h) + '}' + plusLatex(seg.c);
  } else if (seg.type === 'hyperbola') {
    body = '\\frac{' + fmt(p.a) + '}{' + xShiftLatex(p.h) + '}' + plusLatex(seg.c);
  } else if (seg.type === 'exp') {
    const pw = p.h ? '^{' + xShiftLatex(p.h) + '}' : '^{x}';
    if (Math.abs(p.a - 1) < 1e-9) body = fmt(p.b) + pw;
    else if (Math.abs(p.a + 1) < 1e-9) body = '-' + fmt(p.b) + pw;
    else body = fmt(p.a) + '\\cdot' + fmt(p.b) + pw;
    body += plusLatex(seg.c);
  } else if (seg.type === 'log') {
    body = coefLatex(p.a) + '\\ln\\left(' + xShiftLatex(p.h) + '\\right)' + plusLatex(seg.c);
  } else {
    const bstr = Math.abs(p.b - 1) < 1e-9 ? '' : fmt(p.b);
    body = coefLatex(p.a) + '\\sin\\left(' + bstr + '\\left(' + xShiftLatex(p.h) + '\\right)\\right)' + plusLatex(seg.c);
  }
  return 'y=' + body + D;
}

function htmlNum(v) { return fmt(v).replace(/-/g, '−'); }

function xShiftHTML(h) {
  if (Math.abs(h) < 1e-9) return 'x';
  return h > 0 ? 'x − ' + htmlNum(h) : 'x + ' + htmlNum(-h);
}

function plusHTML(v) {
  if (Math.abs(v) < 1e-9) return '';
  return v > 0 ? ' + ' + htmlNum(v) : ' − ' + htmlNum(-v);
}

function coefHTML(a) {
  if (Math.abs(a - 1) < 1e-9) return '';
  if (Math.abs(a + 1) < 1e-9) return '−';
  return htmlNum(a);
}

function equationHTML(seg) {
  const p = seg.params;
  const D = ' <span class="domain">{' + htmlNum(seg.x1) + ' ≤ x ≤ ' + htmlNum(seg.x2) + '}</span>';
  let body;
  if (seg.type === 'poly1' || seg.type === 'poly2' || seg.type === 'poly3') {
    const cs = p.cs;
    const deg = cs.length - 1;
    const parts = [];
    for (let d = deg; d >= 0; d--) {
      const a = cs[deg - d];
      if (Math.abs(a) < 1e-9) continue;
      const suffix = d === 0 ? '' : d === 1 ? 'x' : 'x<sup>' + d + '</sup>';
      const coef = (Math.abs(Math.abs(a) - 1) < 1e-9 && d > 0) ? '' : htmlNum(Math.abs(a));
      parts.push({ neg: a < 0, coef, suffix });
    }
    if (!parts.length) body = '0';
    else body = parts.map((t, i) => (i === 0 ? (t.neg ? '−' : '') : (t.neg ? ' − ' : ' + ')) + t.coef + t.suffix).join('');
  } else if (seg.type === 'sqrt') {
    body = coefHTML(p.a) + '√<span class="rad">' + xShiftHTML(p.h) + '</span>' + plusHTML(seg.c);
  } else if (seg.type === 'hyperbola') {
    body = '<span class="frac"><span>' + htmlNum(p.a) + '</span><span>' + xShiftHTML(p.h) + '</span></span>' + plusHTML(seg.c);
  } else if (seg.type === 'exp') {
    const pw = p.h ? '<sup>' + xShiftHTML(p.h) + '</sup>' : '<sup>x</sup>';
    if (Math.abs(p.a - 1) < 1e-9) body = htmlNum(p.b) + pw;
    else if (Math.abs(p.a + 1) < 1e-9) body = '−' + htmlNum(p.b) + pw;
    else body = htmlNum(p.a) + ' · ' + htmlNum(p.b) + pw;
    body += plusHTML(seg.c);
  } else if (seg.type === 'log') {
    body = coefHTML(p.a) + ' ln(' + xShiftHTML(p.h) + ')' + plusHTML(seg.c);
  } else {
    const bstr = Math.abs(p.b - 1) < 1e-9 ? '' : htmlNum(p.b);
    body = coefHTML(p.a) + ' sin(' + bstr + '(' + xShiftHTML(p.h) + '))' + plusHTML(seg.c);
  }
  return 'y = ' + body + D;
}

function segmentLength(seg) {
  const n = 400;
  let s = 0;
  let px = seg.x1, py = seg.fn(seg.x1);
  for (let i = 1; i <= n; i++) {
    const x = seg.x1 + (seg.x2 - seg.x1) * i / n;
    const y = seg.fn(x);
    s += Math.hypot(x - px, y - py);
    px = x; py = y;
  }
  return s;
}

function findAnchor(rawPts, segments) {
  const p0 = rawPts[0], p1 = rawPts[rawPts.length - 1];
  if (!segments.length) return { target: round3(p0.y), segEnd: 'left', dist: 0 };
  const opts = [];
  for (const seg of segments) {
    opts.push({ d: Math.hypot(p0.x - seg.x1, p0.y - seg.y1), target: seg.y1, segEnd: 'left' });
    opts.push({ d: Math.hypot(p1.x - seg.x1, p1.y - seg.y1), target: seg.y1, segEnd: 'right' });
    opts.push({ d: Math.hypot(p0.x - seg.x2, p0.y - seg.y2), target: seg.y2, segEnd: 'left' });
    opts.push({ d: Math.hypot(p1.x - seg.x2, p1.y - seg.y2), target: seg.y2, segEnd: 'right' });
  }
  opts.sort((a, b) => a.d - b.d);
  return { target: opts[0].target, segEnd: opts[0].segEnd, dist: opts[0].d };
}

function buildSegmentFromStroke(rawPts, anchor, forceType) {
  const pts = resamplePts(rawPts, 90);
  const xs = pts.map(p => p.x), ys = pts.map(p => p.y);
  const x1 = round3(xs[0]), x2 = round3(xs[xs.length - 1]);
  if (!(x2 - x1 >= 0.12)) {
    return { error: 'That stroke is too small or too vertical. Draw a curve at least 0.12 km wide (left to right).' };
  }
  const cands = candidateFitters(xs, ys, x1, x2);
  if (!cands.length) return { error: 'Could not fit a function to that stroke. Try drawing it again.' };
  let pick = null;
  if (forceType && forceType !== 'auto') {
    const group = forceType === 'poly' ? ['poly1', 'poly2', 'poly3']
      : forceType === 'explog' ? ['exp', 'log']
      : [forceType];
    pick = cands.filter(c => group.indexOf(c.type) >= 0).sort((a, b) => a.score - b.score)[0] || null;
    if (!pick) return { error: 'Could not fit the selected function type to that stroke. Try Auto, or redraw the curve.' };
  } else {
    pick = cands.sort((a, b) => a.score - b.score)[0];
  }
  const base = pick.base;
  const snapped = anchor.dist <= SNAP_DIST;
  const c = snapped
    ? (anchor.segEnd === 'left' ? round3(anchor.target - base(x1)) : round3(anchor.target - base(x2)))
    : 0;
  const fn = x => base(x) + c;
  const seg = {
    id: state.nextId++,
    type: pick.type,
    label: pick.label,
    params: pick.params,
    base,
    c,
    fn,
    x1, x2,
    y1: anchor.segEnd === 'left' ? anchor.target : round3(fn(x1)),
    y2: anchor.segEnd === 'right' ? anchor.target : round3(fn(x2)),
    color: PALETTE[(state.segments.length) % PALETTE.length]
  };
  seg.length = segmentLength(seg);
  seg.desmos = desmosString(seg);
  seg.html = equationHTML(seg);
  return { seg };
}

function catOf(type) {
  if (type === 'poly1' || type === 'poly2' || type === 'poly3') return 'poly';
  if (type === 'sqrt') return 'sqrt';
  if (type === 'hyperbola') return 'hyperbola';
  if (type === 'exp' || type === 'log') return 'explog';
  return 'trig';
}

function segmentBoundsOK(seg) {
  const n = 100;
  for (let i = 0; i <= n; i++) {
    const x = seg.x1 + (seg.x2 - seg.x1) * i / n;
    const y = seg.fn(x);
    if (x < -1e-9 || x > PARK_W + 1e-9 || y < -1e-9 || y > PARK_H + 1e-9) return false;
  }
  return true;
}

let els = null, ctx = null, dpr = 1;
let drawing = false, stroke = [], panning = null, draggingFeature = null, dragMoved = false, draggingSeg = null, segMoved = false, segDragY = 0;

function addStroke(rawPts, forceType) {
  const anchor = findAnchor(rawPts, state.segments);
  const res = buildSegmentFromStroke(rawPts, anchor, forceType);
  if (res.error) return res;
  const seg = res.seg;
  const hitFeature = state.features.find(f => segHitsFeature(seg, f));
  if (hitFeature) {
    return { error: 'That section would cross the ' + hitFeature.label.toLowerCase() + ' — draw around it.' };
  }
  state.segments.push(seg);
  state.segments.sort((a, b) => a.x1 - b.x1);
  if (anchor.dist > SNAP_DIST) {
    toast('Drawn free-hand — nearest endpoint is ' + fmt(anchor.dist) + ' km away (within ' + fmt(SNAP_DIST) + ' km it would have snapped).');
  }
  if (state.segments.length > 1) {
    const i = state.segments.indexOf(seg);
    const prev = state.segments[i - 1];
    if (prev && Math.abs(prev.y2 - seg.y1) > 1e-9) {
      toast('Warning: this section was inserted between existing sections — check continuity in the checklist.');
    }
  }
  if (!segmentBoundsOK(seg)) {
    toast('Warning: part of this section lies outside the 3 km × 2 km park.');
  }
  refresh();
  return { seg };
}

function loadExample() {
  state.segments = [];
  state.isolateId = null;
  const strokes = [];
  let s = [];
  for (let i = 0; i <= 40; i++) {
    const x = 0.02 + (0.68 - 0.02) * i / 40;
    s.push({ x, y: 0.4 * x * x + 0.5 });
  }
  strokes.push({ pts: s, type: 'poly' });
  const y1 = 0.4 * 0.68 * 0.68 + 0.5;
  s = [];
  for (let i = 0; i <= 40; i++) {
    const x = 0.68 + (1.7 - 0.68) * i / 40;
    s.push({ x, y: y1 + 0.55 * Math.sqrt(x - 0.68) });
  }
  strokes.push({ pts: s, type: 'sqrt' });
  const y2 = y1 + 0.55;
  s = [];
  for (let i = 0; i <= 48; i++) {
    const x = 1.7 + (2.3 - 1.7) * i / 48;
    s.push({ x, y: y2 + 0.1 * Math.sin(3 * (x - 1.7)) });
  }
  strokes.push({ pts: s, type: 'sine' });
  const y3 = y2 + 0.1 * Math.sin(3 * 0.6);
  const k4 = y3 - 0.3;
  s = [];
  for (let i = 0; i <= 40; i++) {
    const x = 2.3 + (2.8 - 2.3) * i / 40;
    s.push({ x, y: k4 + 0.15 / (x - 1.8) });
  }
  strokes.push({ pts: s, type: 'hyperbola' });
  const y4 = k4 + 0.15;
  const c5 = y4 - 0.02 * Math.pow(1.6, 2.8);
  s = [];
  for (let i = 0; i <= 30; i++) {
    const x = 2.8 + (3.0 - 2.8) * i / 30;
    s.push({ x, y: c5 + 0.02 * Math.pow(1.6, x) });
  }
  strokes.push({ pts: s, type: 'exp' });
  for (const st of strokes) {
    const anchor = findAnchor(st.pts, state.segments);
    const res = buildSegmentFromStroke(st.pts, anchor, st.type);
    if (res.seg) {
      state.segments.push(res.seg);
      state.segments.sort((a, b) => a.x1 - b.x1);
    }
  }
  refresh();
  toast('Example track loaded — 5 sections, all function types, fully continuous.');
}

function worldFromEvent(e) {
  const rect = els.canvas.getBoundingClientRect();
  const px = e.clientX - rect.left, py = e.clientY - rect.top;
  const v = state.view;
  return { x: (px - rect.width / 2) / v.ppk + v.cx, y: (rect.height / 2 - py) / v.ppk + v.cy };
}

function pointInFeature(f, x, y) {
  if (f.type === 'ellipse' || f.type === 'circle') {
    const dx = (x - f.cx) / f.rx, dy = (y - f.cy) / f.ry;
    return dx * dx + dy * dy <= 1;
  }
  let inside = false;
  const pts = f.pts;
  for (let a = 0, b = pts.length - 1; a < pts.length; b = a++) {
    if ((pts[a].y > y) !== (pts[b].y > y) &&
        x < (pts[b].x - pts[a].x) * (y - pts[a].y) / (pts[b].y - pts[a].y) + pts[a].x) inside = !inside;
  }
  return inside;
}

function featureAt(w) {
  for (let i = state.features.length - 1; i >= 0; i--) {
    if (pointInFeature(state.features[i], w.x, w.y)) return state.features[i];
  }
  return null;
}

function featuresOverlap(a, b) {
  const samplesOf = f => {
    const out = [];
    if (f.type === 'ellipse' || f.type === 'circle') {
      for (let i = 0; i < 24; i++) {
        const t = i / 24 * Math.PI * 2;
        out.push({ x: f.cx + f.rx * Math.cos(t), y: f.cy + f.ry * Math.sin(t) });
      }
      out.push({ x: f.cx, y: f.cy });
    } else {
      for (const p of f.pts) out.push({ x: p.x, y: p.y });
      out.push({
        x: f.pts.reduce((s, p) => s + p.x, 0) / f.pts.length,
        y: f.pts.reduce((s, p) => s + p.y, 0) / f.pts.length
      });
    }
    return out;
  };
  const A = samplesOf(a), B = samplesOf(b);
  return A.some(p => pointInFeature(b, p.x, p.y)) || B.some(p => pointInFeature(a, p.x, p.y));
}

function segHitsFeature(seg, f) {
  const n = 80;
  for (let i = 0; i <= n; i++) {
    const x = seg.x1 + (seg.x2 - seg.x1) * i / n;
    if (pointInFeature(f, x, seg.fn(x))) return true;
  }
  return false;
}

function featureOverlapsAnything(f, exclude) {
  for (const o of state.features) {
    if (o === f || o === exclude) continue;
    if (featuresOverlap(f, o)) return o.label.toLowerCase();
  }
  for (const seg of state.segments) {
    if (segHitsFeature(seg, f)) return 'a track section';
  }
  return null;
}

function segmentAt(w) {
  const v = state.view;
  let best = null, bestD = Infinity;
  for (const seg of state.segments) {
    const n = 120;
    for (let i = 0; i <= n; i++) {
      const x = seg.x1 + (seg.x2 - seg.x1) * i / n;
      const d = Math.hypot((w.x - x) * v.ppk, (w.y - seg.fn(x)) * v.ppk);
      if (d < bestD) { bestD = d; best = seg; }
    }
  }
  return bestD <= 8 ? best : null;
}

function moveFeature(f, dx, dy) {
  f.cx += dx; f.cy += dy;
  if (f.pts) f.pts = f.pts.map(p => ({ x: p.x + dx, y: p.y + dy }));
}

function featureMoved(f, dx, dy) {
  const cand = Object.assign({}, f);
  cand.cx = f.cx + dx;
  cand.cy = f.cy + dy;
  if (f.pts) cand.pts = f.pts.map(p => ({ x: p.x + dx, y: p.y + dy }));
  return cand;
}

function featurePosOK(f, dx, dy) {
  return !featureOverlapsAnything(featureMoved(f, dx, dy), f);
}

function findFreeSpot(f) {
  if (featurePosOK(f, 0, 0)) return { x: f.cx, y: f.cy };
  for (let r = 0.15; r <= 1.25; r += 0.15) {
    for (let a = 0; a < 12; a++) {
      const ang = a / 12 * Math.PI * 2;
      const dx = r * Math.cos(ang), dy = r * Math.sin(ang);
      const nx = f.cx + dx, ny = f.cy + dy;
      if (nx < 0.15 || nx > PARK_W - 0.15 || ny < 0.15 || ny > PARK_H - 0.15) continue;
      if (featurePosOK(f, dx, dy)) return { x: nx, y: ny };
    }
  }
  return null;
}

function segShiftValid(seg, newC) {
  const test = { x1: seg.x1, x2: seg.x2, fn: x => seg.base(x) + newC };
  for (const f of state.features) {
    if (segHitsFeature(test, f)) return false;
  }
  return true;
}

function snappedShift(seg, newC) {
  let best = null;
  for (const o of state.segments) {
    if (o === seg) continue;
    for (const ex of [seg.x1, seg.x2]) {
      for (const pair of [[o.x1, o.y1], [o.x2, o.y2]]) {
        if (Math.abs(pair[0] - ex) >= 0.02) continue;
        const dy = pair[1] - (seg.base(ex) + newC);
        if (Math.abs(dy) <= 0.08 && (best === null || Math.abs(dy) < Math.abs(best))) best = dy;
      }
    }
  }
  return best === null ? newC : newC + best;
}

function applySegShift(seg, newC) {
  seg.c = newC;
  seg.fn = x => seg.base(x) + seg.c;
  seg.y1 = round3(seg.fn(seg.x1));
  seg.y2 = round3(seg.fn(seg.x2));
}

function draw() {
  const canvas = els.canvas;
  const W = canvas.clientWidth, H = canvas.clientHeight;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, W, H);
  const v = state.view;
  const sx = x => (x - v.cx) * v.ppk + W / 2;
  const sy = y => H / 2 - (y - v.cy) * v.ppk;
  const wx = px => (px - W / 2) / v.ppk + v.cx;
  const wy = py => (H / 2 - py) / v.ppk + v.cy;

  const steps = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 5, 10];
  const major = steps.find(s => s * v.ppk >= 55) || 10;
  const minor = major / 5;
  const x0 = wx(0), x1 = wx(W), y0 = wy(H), y1v = wy(0);

  ctx.lineWidth = 1;
  ctx.strokeStyle = '#f0f2f5';
  ctx.beginPath();
  for (let x = Math.ceil(x0 / minor) * minor; x <= x1; x += minor) { ctx.moveTo(sx(x), 0); ctx.lineTo(sx(x), H); }
  for (let y = Math.ceil(y0 / minor) * minor; y <= y1v; y += minor) { ctx.moveTo(0, sy(y)); ctx.lineTo(W, sy(y)); }
  ctx.stroke();

  ctx.strokeStyle = '#dce1e7';
  ctx.beginPath();
  for (let x = Math.ceil(x0 / major) * major; x <= x1; x += major) { ctx.moveTo(sx(x), 0); ctx.lineTo(sx(x), H); }
  for (let y = Math.ceil(y0 / major) * major; y <= y1v; y += major) { ctx.moveTo(0, sy(y)); ctx.lineTo(W, sy(y)); }
  ctx.stroke();

  ctx.fillStyle = 'rgba(144, 190, 109, 0.10)';
  ctx.fillRect(sx(0), sy(PARK_H), PARK_W * v.ppk, PARK_H * v.ppk);
  ctx.strokeStyle = '#7fb069';
  ctx.lineWidth = 1.5;
  ctx.strokeRect(sx(0), sy(PARK_H), PARK_W * v.ppk, PARK_H * v.ppk);

  ctx.strokeStyle = '#9aa4ae';
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(sx(0), 0); ctx.lineTo(sx(0), H);
  ctx.moveTo(0, sy(0)); ctx.lineTo(W, sy(0));
  ctx.stroke();

  ctx.fillStyle = '#6b7680';
  ctx.font = '12px Arial';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  const yAxisVisible = sx(0) >= 0 && sx(0) <= W;
  const xAxisVisible = sy(0) >= 0 && sy(0) <= H;
  for (let x = Math.ceil(x0 / major) * major; x <= x1; x += major) {
    if (Math.abs(x) < 1e-9) continue;
    const ly = xAxisVisible ? sy(0) + 5 : (wy(0) > 0 ? H - 18 : 4);
    ctx.fillText(fmt(x, 2), sx(x), ly);
  }
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  for (let y = Math.ceil(y0 / major) * major; y <= y1v; y += major) {
    if (Math.abs(y) < 1e-9) continue;
    const lx = yAxisVisible ? sx(0) - 6 : (wx(0) > 0 ? 30 : W - 6);
    ctx.fillText(fmt(y, 2), lx, sy(y));
  }
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText('km', W - 26, (xAxisVisible ? sy(0) : H - 20) + 4);

  ctx.font = 'bold 12px Arial';
  ctx.fillStyle = '#4a7c3f';
  ctx.fillText('Park boundary: 3 km × 2 km', sx(0) + 6, sy(PARK_H) + 5);

  for (const f of state.features) {
    ctx.beginPath();
    if (f.type === 'ellipse' || f.type === 'circle') {
      ctx.ellipse(sx(f.cx), sy(f.cy), f.rx * v.ppk, f.ry * v.ppk, 0, 0, Math.PI * 2);
    } else {
      f.pts.forEach((p, i) => i ? ctx.lineTo(sx(p.x), sy(p.y)) : ctx.moveTo(sx(p.x), sy(p.y)));
      ctx.closePath();
    }
    ctx.fillStyle = f.fill;
    ctx.fill();
    ctx.strokeStyle = f.stroke;
    ctx.lineWidth = 1.5;
    ctx.stroke();
    if (state.selectedFeature === f) {
      ctx.setLineDash([5, 4]);
      ctx.strokeStyle = '#333';
      ctx.lineWidth = 1.5;
      ctx.stroke();
      ctx.setLineDash([]);
    }
    let lx = sx(f.cx), ly = sy(f.cy);
    if (f.type === 'polygon') {
      const c = f.pts.reduce((a, p) => ({ x: a.x + p.x / f.pts.length, y: a.y + p.y / f.pts.length }), { x: 0, y: 0 });
      lx = sx(c.x); ly = sy(c.y);
    }
    ctx.font = 'bold 12px Arial';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.lineWidth = 3;
    ctx.strokeStyle = 'rgba(255,255,255,0.9)';
    ctx.strokeText(f.label, lx, ly);
    ctx.fillStyle = '#333';
    ctx.fillText(f.label, lx, ly);
  }

  const isolated = state.isolateId;
  for (const seg of state.segments) {
    const dim = isolated && seg.id !== isolated;
    ctx.beginPath();
    const n = 200;
    for (let i = 0; i <= n; i++) {
      const x = seg.x1 + (seg.x2 - seg.x1) * i / n;
      const X = sx(x), Y = sy(seg.fn(x));
      if (i) ctx.lineTo(X, Y); else ctx.moveTo(X, Y);
    }
    if (state.selectedSeg === seg) {
      ctx.globalAlpha = dim ? 0.15 : 0.25;
      ctx.strokeStyle = '#333';
      ctx.lineWidth = 9;
      ctx.stroke();
    }
    ctx.strokeStyle = seg.color;
    ctx.globalAlpha = dim ? 0.15 : 1;
    ctx.lineWidth = dim ? 2 : 3;
    ctx.stroke();
  }

  ctx.globalAlpha = 1;
  const sorted = state.segments.slice().sort((a, b) => a.x1 - b.x1);
  const endpoints = [];
  sorted.forEach((seg, i) => {
    if (i === 0) endpoints.push({ x: seg.x1, y: seg.y1 });
    endpoints.push({ x: seg.x2, y: seg.y2 });
  });
  ctx.font = '11px Arial';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'bottom';
  for (const p of endpoints) {
    const dim = isolated && !sorted.some(sg => sg.id === isolated && ((Math.abs(sg.x1 - p.x) < 1e-9 && Math.abs(sg.y1 - p.y) < 1e-9) || (Math.abs(sg.x2 - p.x) < 1e-9 && Math.abs(sg.y2 - p.y) < 1e-9)));
    ctx.globalAlpha = dim ? 0.15 : 1;
    ctx.beginPath();
    ctx.arc(sx(p.x), sy(p.y), 4, 0, Math.PI * 2);
    ctx.fillStyle = '#fff';
    ctx.fill();
    ctx.lineWidth = 2;
    ctx.strokeStyle = '#333';
    ctx.stroke();
    const label = '(' + fmt(p.x) + ', ' + fmt(p.y) + ')';
    ctx.lineWidth = 3;
    ctx.strokeStyle = 'rgba(255,255,255,0.9)';
    ctx.strokeText(label, sx(p.x) + 7, sy(p.y) - 6);
    ctx.fillStyle = '#333';
    ctx.fillText(label, sx(p.x) + 7, sy(p.y) - 6);
  }
  ctx.globalAlpha = 1;

  if (drawing && stroke.length > 1) {
    ctx.setLineDash([6, 4]);
    ctx.strokeStyle = '#888';
    ctx.lineWidth = 2;
    ctx.beginPath();
    stroke.forEach((p, i) => {
      const X = sx(p.x), Y = sy(p.y);
      if (i) ctx.lineTo(X, Y); else ctx.moveTo(X, Y);
    });
    ctx.stroke();
    ctx.setLineDash([]);
  }
}

function refresh() {
  renderSidebar();
  draw();
}

function renderSidebar() {
  const list = els.segList;
  list.innerHTML = '';
  const sorted = state.segments.slice().sort((a, b) => a.x1 - b.x1);
  sorted.forEach(seg => {
    const div = document.createElement('div');
    div.className = 'seg-item' + (state.isolateId === seg.id ? ' isolated' : '');
    div.style.borderLeftColor = seg.color;
    const head = document.createElement('div');
    head.className = 'seg-head';
    const dot = document.createElement('span');
    dot.className = 'seg-dot';
    dot.style.background = seg.color;
    const name = document.createElement('span');
    name.className = 'seg-name';
    name.textContent = 'Section ' + (sorted.indexOf(seg) + 1) + ' — ' + seg.label;
    const copy = document.createElement('button');
    copy.className = 'seg-btn';
    copy.textContent = 'Copy';
    copy.addEventListener('click', e => { e.stopPropagation(); copyText(seg.desmos); });
    const del = document.createElement('button');
    del.className = 'seg-btn';
    del.textContent = '✕';
    del.addEventListener('click', e => {
      e.stopPropagation();
      state.segments = state.segments.filter(s => s.id !== seg.id);
      if (state.isolateId === seg.id) state.isolateId = null;
      if (state.selectedSeg === seg) state.selectedSeg = null;
      refresh();
    });
    head.appendChild(dot);
    head.appendChild(name);
    head.appendChild(copy);
    head.appendChild(del);
    const eq = document.createElement('div');
    eq.className = 'eq';
    eq.innerHTML = seg.html;
    const meta = document.createElement('div');
    meta.className = 'seg-meta';
    meta.textContent = 'Length: ' + fmt(seg.length) + ' km · endpoints (' + fmt(seg.x1) + ', ' + fmt(seg.y1) + ') → (' + fmt(seg.x2) + ', ' + fmt(seg.y2) + ')';
    div.appendChild(head);
    div.appendChild(eq);
    div.appendChild(meta);
    div.addEventListener('click', () => {
      state.isolateId = state.isolateId === seg.id ? null : seg.id;
      refresh();
    });
    list.appendChild(div);
  });
  if (!sorted.length) {
    list.innerHTML = '<div class="feat-hint">No sections yet. Left-drag on the graph to draw one — it will be converted into an equation automatically.</div>';
  }
  renderChecklist();
}

function renderChecklist() {
  const el = els.checklist;
  const sorted = state.segments.slice().sort((a, b) => a.x1 - b.x1);
  const cats = new Set(sorted.map(s => catOf(s.type)));
  const rows = [
    ['poly', 'Polynomial function (Study Area 1)'],
    ['sqrt', 'Square root function (Study Area 1)'],
    ['hyperbola', 'Rectangular hyperbola (Study Area 1)'],
    ['explog', 'Exponential or logarithmic (Study Area 2)'],
    ['trig', 'Sine or cosine function (Study Area 2)']
  ];
  let html = '';
  for (const [key, text] of rows) {
    const ok = cats.has(key);
    html += '<div class="check-row"><span class="' + (ok ? 'check-ok">✓' : 'check-no">✗') + '</span><span>' + text + '</span></div>';
  }
  let cont = true;
  for (let i = 1; i < sorted.length; i++) {
    if (Math.abs(sorted[i - 1].y2 - sorted[i].y1) > 1e-9) cont = false;
  }
  const contText = sorted.length < 2 ? 'add more sections' : (cont ? 'all junctions meet exactly' : 'mismatch at a junction');
  html += '<div class="check-row"><span class="' + (sorted.length >= 2 && cont ? 'check-ok">✓' : 'check-no">✗') + '</span><span>Continuity at every junction — ' + contText + '</span></div>';
  const inBounds = sorted.length > 0 && sorted.every(segmentBoundsOK);
  html += '<div class="check-row"><span class="' + (inBounds ? 'check-ok">✓' : 'check-no">✗') + '</span><span>Track inside the 3 km × 2 km park, first quadrant</span></div>';
  html += '<div class="check-row"><span class="' + (state.features.length >= 3 ? 'check-ok">✓' : 'check-no">✗') + '</span><span>Lake, playground and café placed (shapes: ' + state.features.length + ')</span></div>';
  const total = sorted.reduce((a, s) => a + s.length, 0);
  els.totalLen.textContent = sorted.length
    ? 'Total track length: ' + fmt(total) + ' km' + (total >= 2 ? '' : ' — aim for a reasonable running distance (≈ 3–6 km)')
    : '';
  el.innerHTML = html;
}

function copyText(text) {
  const done = () => toast('Copied — paste straight into Desmos.');
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(done, () => fallbackCopy(text, done));
  } else {
    fallbackCopy(text, done);
  }
}

function fallbackCopy(text, done) {
  const ta = document.createElement('textarea');
  ta.value = text;
  ta.style.position = 'fixed';
  ta.style.opacity = '0';
  document.body.appendChild(ta);
  ta.select();
  try { document.execCommand('copy'); done(); } catch (e) { toast('Copy failed — select the equation text manually.'); }
  document.body.removeChild(ta);
}

let toastTimer = null;
function toast(msg) {
  els.toast.textContent = msg;
  els.toast.classList.add('show');
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => els.toast.classList.remove('show'), 3200);
}

function resetView() {
  const W = els.canvas.clientWidth, H = els.canvas.clientHeight;
  state.view.ppk = Math.min(W / (PARK_W + 0.7), H / (PARK_H + 0.7));
  state.view.cx = PARK_W / 2;
  state.view.cy = PARK_H / 2;
  draw();
}

function resize() {
  const canvas = els.canvas;
  dpr = window.devicePixelRatio || 1;
  const W = canvas.clientWidth, H = canvas.clientHeight;
  canvas.width = Math.round(W * dpr);
  canvas.height = Math.round(H * dpr);
  draw();
}

function makeFeature(kind, at) {
  const p = at || { x: PARK_W / 2, y: PARK_H / 2 };
  if (kind === 'lake') return { id: state.nextId++, kind, label: 'Lake', type: 'ellipse', cx: p.x, cy: p.y, rx: 0.4, ry: 0.28, fill: 'rgba(70,130,180,0.30)', stroke: '#2d70b3' };
  if (kind === 'cafe') return { id: state.nextId++, kind, label: 'Café', type: 'circle', cx: p.x, cy: p.y, rx: 0.14, ry: 0.14, fill: 'rgba(139,69,19,0.28)', stroke: '#8b4513' };
  return {
    id: state.nextId++, kind, label: 'Playground', type: 'polygon',
    cx: p.x, cy: p.y,
    pts: [{ x: p.x - 0.3, y: p.y - 0.15 }, { x: p.x + 0.3, y: p.y - 0.15 }, { x: p.x + 0.3, y: p.y + 0.15 }, { x: p.x - 0.3, y: p.y + 0.15 }],
    fill: 'rgba(250,160,30,0.25)', stroke: '#e08000'
  };
}

function undo() {
  if (!state.segments.length) return;
  const sorted = state.segments.slice().sort((a, b) => a.x1 - b.x1);
  const last = sorted[sorted.length - 1];
  state.segments = state.segments.filter(s => s !== last);
  if (state.isolateId === last.id) state.isolateId = null;
  if (state.selectedSeg === last) state.selectedSeg = null;
  refresh();
}

function init() {
  els = {
    canvas: document.getElementById('graph'),
    segList: document.getElementById('segList'),
    checklist: document.getElementById('checklist'),
    totalLen: document.getElementById('totalLen'),
    toast: document.getElementById('toast'),
    forceType: document.getElementById('forceType')
  };
  ctx = els.canvas.getContext('2d');

  state.features = [
    makeFeature('lake', { x: 2.05, y: 1.0 }),
    makeFeature('playground', { x: 0.5, y: 1.7 }),
    makeFeature('cafe', { x: 0.55, y: 0.45 })
  ];

  els.canvas.addEventListener('contextmenu', e => e.preventDefault());

  els.canvas.addEventListener('mousedown', e => {
    if (e.button === 2 || e.button === 1) {
      panning = { x: e.clientX, y: e.clientY };
      els.canvas.classList.add('panning');
      e.preventDefault();
      return;
    }
    if (e.button !== 0) return;
    const w = worldFromEvent(e);
    if (state.mode === 'edit') {
      const f = featureAt(w);
      state.selectedFeature = f;
      state.selectedSeg = null;
      if (f) {
        draggingFeature = f;
        draggingFeature._off = { x: w.x - f.cx, y: w.y - f.cy };
      } else {
        const seg = segmentAt(w);
        state.selectedSeg = seg;
        if (seg) {
          draggingSeg = seg;
          segDragY = w.y;
        }
      }
      draw();
      return;
    }
    state.selectedFeature = null;
    state.selectedSeg = null;
    drawing = true;
    stroke = [w];
    draw();
  });

  window.addEventListener('mousemove', e => {
    if (panning) {
      state.view.cx -= (e.clientX - panning.x) / state.view.ppk;
      state.view.cy += (e.clientY - panning.y) / state.view.ppk;
      panning = { x: e.clientX, y: e.clientY };
      draw();
      return;
    }
    if (draggingFeature) {
      const w = worldFromEvent(e);
      const f = draggingFeature;
      const nx = w.x - f._off.x, ny = w.y - f._off.y;
      const dx = nx - f.cx, dy = ny - f.cy;
      if (dx || dy) dragMoved = true;
      let applied = false;
      if ((dx || dy) && featurePosOK(f, dx, dy)) {
        moveFeature(f, dx, dy);
        applied = true;
      } else if (dx && featurePosOK(f, dx, 0)) {
        moveFeature(f, dx, 0);
        applied = true;
      } else if (dy && featurePosOK(f, 0, dy)) {
        moveFeature(f, 0, dy);
        applied = true;
      }
      if ((dx || dy) && !applied) {
        const hit = featureOverlapsAnything(featureMoved(f, dx, dy), f);
        if (hit) toast('Blocked — the ' + f.label.toLowerCase() + ' would overlap ' + hit + '.');
      }
      draw();
      return;
    }
    if (draggingSeg) {
      const w = worldFromEvent(e);
      const dy = w.y - segDragY;
      if (dy) {
        segMoved = true;
        let target = snappedShift(draggingSeg, draggingSeg.c + dy);
        if (!segShiftValid(draggingSeg, target)) target = draggingSeg.c + dy;
        if (segShiftValid(draggingSeg, target)) {
          applySegShift(draggingSeg, target);
          segDragY = w.y;
        } else {
          toast('Blocked — the section would overlap a park feature.');
        }
        draw();
      }
      return;
    }
    if (drawing) {
      const w = worldFromEvent(e);
      const last = stroke[stroke.length - 1];
      if (Math.hypot(w.x - last.x, w.y - last.y) * state.view.ppk > 2 && stroke.length < 2000) {
        stroke.push(w);
        draw();
      }
    }
  });

  window.addEventListener('mouseup', () => {
    if (panning) {
      panning = null;
      els.canvas.classList.remove('panning');
    }
    if (draggingFeature) {
      if (dragMoved) {
        state.selectedFeature = null;
        draw();
      }
      draggingFeature = null;
      dragMoved = false;
    }
    if (draggingSeg) {
      if (segMoved) {
        const s = draggingSeg;
        s.length = segmentLength(s);
        s.desmos = desmosString(s);
        s.html = equationHTML(s);
        state.selectedSeg = null;
        refresh();
      }
      draggingSeg = null;
      segMoved = false;
    }
    if (drawing) {
      drawing = false;
      const pts = stroke;
      stroke = [];
      if (pts.length >= 5) {
        const res = addStroke(pts, els.forceType.value);
        if (res.error) toast(res.error);
      } else {
        draw();
      }
    }
  });

  els.canvas.addEventListener('wheel', e => {
    e.preventDefault();
    const rect = els.canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left, my = e.clientY - rect.top;
    const v = state.view;
    const bx = (mx - rect.width / 2) / v.ppk + v.cx;
    const by = (rect.height / 2 - my) / v.ppk + v.cy;
    const factor = Math.exp(-e.deltaY * 0.0012);
    v.ppk = Math.max(20, Math.min(2000, v.ppk * factor));
    const ax = (mx - rect.width / 2) / v.ppk + v.cx;
    const ay = (rect.height / 2 - my) / v.ppk + v.cy;
    v.cx += bx - ax;
    v.cy += by - ay;
    draw();
  }, { passive: false });

  window.addEventListener('keydown', e => {
    const ae = document.activeElement;
    const tag = ae && ae.tagName;
    if (tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA' || (ae && ae.isContentEditable)) return;
    if (e.key === 'Delete' && state.selectedFeature) {
      const f = state.selectedFeature;
      state.features = state.features.filter(x => x !== f);
      state.selectedFeature = null;
      toast(f.label + ' deleted');
      draw();
      e.preventDefault();
    } else if (e.key === 'Delete' && state.selectedSeg) {
      const s = state.selectedSeg;
      state.segments = state.segments.filter(x => x !== s);
      if (state.isolateId === s.id) state.isolateId = null;
      state.selectedSeg = null;
      toast(s.label + ' deleted');
      refresh();
      e.preventDefault();
    }
    if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'z') {
      e.preventDefault();
      undo();
    }
  });

  function setMode(m) {
    state.mode = m;
    document.getElementById('modeDraw').classList.toggle('active', m === 'draw');
    document.getElementById('modeEdit').classList.toggle('active', m === 'edit');
    els.canvas.classList.toggle('editing', m === 'edit');
    if (m === 'edit') toast('Edit mode: drag features, or click a section to move it. Delete removes the selection.');
    else toast('Draw mode: left-drag to draw a track section.');
  }
  document.getElementById('modeDraw').addEventListener('click', () => setMode('draw'));
  document.getElementById('modeEdit').addEventListener('click', () => setMode('edit'));

  document.querySelectorAll('[data-add]').forEach(btn => {
    btn.addEventListener('click', () => {
      const f = makeFeature(btn.dataset.add);
      const spot = findFreeSpot(f);
      if (!spot) {
        toast('No room — a ' + f.label.toLowerCase() + ' would overlap something everywhere.');
        return;
      }
      moveFeature(f, spot.x - f.cx, spot.y - f.cy);
      state.features.push(f);
      draw();
    });
  });

  document.getElementById('btnUndo').addEventListener('click', undo);
  document.getElementById('btnClear').addEventListener('click', () => {
    state.segments = [];
    state.isolateId = null;
    state.selectedSeg = null;
    refresh();
  });
  document.getElementById('btnExample').addEventListener('click', loadExample);
  document.getElementById('btnCopyAll').addEventListener('click', () => {
    const sorted = state.segments.slice().sort((a, b) => a.x1 - b.x1);
    if (!sorted.length) { toast('Nothing to copy yet — draw some track first.'); return; }
    copyText(sorted.map(s => s.desmos).join('\n'));
  });
  document.getElementById('btnPng').addEventListener('click', () => {
    els.canvas.toBlob(blob => {
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = 'park-track.png';
      a.click();
      URL.revokeObjectURL(a.href);
    });
  });
  document.getElementById('zin').addEventListener('click', () => { state.view.ppk *= 1.25; draw(); });
  document.getElementById('zout').addEventListener('click', () => { state.view.ppk /= 1.25; draw(); });
  document.getElementById('zreset').addEventListener('click', resetView);

  window.addEventListener('resize', resize);
  resize();
  resetView();
  refresh();
}

if (typeof document !== 'undefined' && document.getElementById) {
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    round3, fmt, resamplePts, candidateFitters, findAnchor, buildSegmentFromStroke,
    segmentLength, segmentBoundsOK, catOf, desmosString, equationHTML, state,
    PARK_W, PARK_H
  };
}
