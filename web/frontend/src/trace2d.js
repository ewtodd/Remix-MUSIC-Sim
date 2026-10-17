// 2D trace view: read `events_MeV` from the served .root via JSROOT and render
// the per-strip dE curve on a canvas. Reconstructed from the bulk tree (the
// `traces/` dir is only written under run.update=1, which fails headless).
//
// Data model (verified):
//   events_MeV:  330 events (events_per_strip * n_strips)
//    - LeftdE[16], RightdE[16]: L/R readout per strip 1..16, interleaved
//    - Strip0dE, Strip17dE, Cathode: scalars per event
//    - strip s trace = Strip0dE (s=0) | LeftdE[s-1]+RightdE[s-1] (1..16) | Strip17dE (s=17)

import { openFile, treeDraw } from "jsroot";

const N_MID = 16; // strips 1..16 have L/R electrodes

// Pull a branch into a Float array (scalar) or Flattened-by-event array.
async function pull(handle, expr) {
  const res = await treeDraw(handle, {
    expr,
    graph: true,
    arr_limit: 10000000,
    want_hist: false,
  });
  return res && res.fY ? res.fY : new Float64Array(0);
}

// Gather the full events_MeV column data + MC event mapping in one call.
export async function loadTraceData(rootHandle) {
  const ev = await rootHandle.readObject("events_MeV");
  const [LeftdE, RightdE, Strip0dE, Strip17dE, Cathode] = await Promise.all([
    pull(ev, "LeftdE"),
    pull(ev, "RightdE"),
    pull(ev, "Strip0dE"),
    pull(ev, "Strip17dE"),
    pull(ev, "Cathode"),
  ]);

  const nEvents = Cathode.length;
  const events = [];
  for (let e = 0; e < nEvents; e++) {
    const strip = [Strip0dE[e], ...new Array(N_MID), Strip17dE[e]];
    // For strips 1..16, L and R arrays are flattened by event (16 each).
    for (let s = 0; s < N_MID; s++) {
      strip[1 + s] = LeftdE[e * N_MID + s] + RightdE[e * N_MID + s];
    }
    events.push({ strip, cathode: Cathode[e] });
  }

  // Map event index -> reaction strip from MC (optional; for annotations).
  let reactionStrip = null;
  let totalEntries = null;
  try {
    const mc = await rootHandle.readObject("MC");
    const rs = await pull(mc, "reaction_strip");
    const n = await pull(mc, "n_steps");
    reactionStrip = rs;
    totalEntries = n.length;
  } catch { /* MC optional */ }

  return { nEvents, events, reactionStrip, totalEntries };
}

// Render one event's per-strip dE curve onto ctx (HTMLCanvasElement 2D).
export function drawTrace(ctx, strip, cathode, opts = {}) {
  const { width = 640, height = 360, bg = "#0b0f18" } = opts;
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = bg;
  ctx.fillRect(0, 0, width, height);

  const n = strip.length; // 18 strips (0..17)
  const pad = { l: 46, r: 14, t: 16, b: 30 };
  const plotW = width - pad.l - pad.r;
  const plotH = height - pad.t - pad.b;
  const xAt = (i) => pad.l + (i / (n - 1)) * plotW;
  const maxY = Math.max(...strip, 0) * 1.1 || 1;
  const yAt = (v) => pad.t + plotH - (v / maxY) * plotH;

  // Axes
  ctx.strokeStyle = "#2c3852";
  ctx.fillStyle = "#93a2bd";
  ctx.font = "11px system-ui";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t);
  ctx.lineTo(pad.l, pad.t + plotH);
  ctx.lineTo(pad.l + plotW, pad.t + plotH);
  ctx.stroke();

  // Y labels
  for (let gy = 0; gy <= 4; gy++) {
    const v = (maxY / 4) * gy;
    const y = yAt(v);
    ctx.fillText(v.toFixed(2), 8, y + 4);
    ctx.beginPath();
    ctx.moveTo(pad.l, y);
    ctx.lineTo(pad.l + plotW, y);
    ctx.strokeStyle = "#1a2233";
    ctx.stroke();
  }

  // Trace line + points
  ctx.strokeStyle = "#7ad0ff";
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let i = 0; i < n; i++) {
    const x = xAt(i), y = yAt(strip[i]);
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();
  ctx.fillStyle = "#4c8dff";
  for (let i = 0; i < n; i++) {
    ctx.beginPath();
    ctx.arc(xAt(i), yAt(strip[i]), 3, 0, Math.PI * 2);
    ctx.fill();
  }

  // X labels (strip index) + title
  ctx.fillStyle = "#93a2bd";
  for (let i = 0; i < n; i++) {
    ctx.fillText(String(i), xAt(i) - 6, height - 10);
  }
  ctx.fillStyle = "#e6ecf5";
  ctx.fillText(`Strip dE (MeV)   cathode ${cathode.toFixed(2)} MeV`, 12, 12);
}
