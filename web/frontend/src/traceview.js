// Trace view controller: fetch the served .root, load events_MeV via JSROOT,
// and draw the trace for the selected event. Handles the slider + normalize
// toggle. Called from main.js after a run reaches `done`.

import { openResultRoot } from "./rootio.js";
import { loadTraceData, drawTrace } from "./trace2d.js";
import { computeGains, applyGains } from "./normalize.js";

const $$ = (sel) => document.querySelector(sel);

let traceData = null;
let gains = null;

export async function showTrace(runId) {
  const vis = $$("#results");
  const canvas = $$("#trace-canvas");
  const ctx = canvas.getContext("2d");
  ctx.fillStyle = "#0b0f18";
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  try {
    // Reveal the trace body; front-and-center.
    vis.querySelector(".vis-placeholder").hidden = true;
    vis.querySelector(".vis-body").hidden = false;

    // Fetch the result .root and open it with JSROOT in the browser.
    const root = await openResultRoot(runId);

    traceData = await loadTraceData(root);
    const { nEvents, reactionStrip } = traceData;

    if (!nEvents) {
      throw new Error("no events in events_MeV");
    }

    // Build per-strip gains from the unreacted (beam) events in this run.
    const beamTraces = [];
    const rs = reactionStrip;
    for (let e = 0; e < nEvents && rs; e++) {
      if (rs[e] === -1) beamTraces.push(traceData.events[e].strip);
    }
    gains = computeGains(beamTraces, 18);

    const slider = $$("#event-slider");
    slider.max = nEvents - 1;
    slider.value = 0;
    $$("#event-label").textContent = "0";
    $$("#trace-meta").textContent =
      `n=${nEvents} | beam(react_strip=-1)=${beamTraces.length}` +
      (gains ? ` | gains applied` : " | no beam events (gains=1)");

    const draw = () => {
      const idx = Number(slider.value);
      let strip = traceData.events[idx].strip;
      if ($$("#normalize-view").checked && gains) strip = applyGains(strip, gains);
      drawTrace(ctx, strip, traceData.events[idx].cathode);
    };
    slider.addEventListener("input", draw);
    $$("#normalize-view").addEventListener("change", draw);
    draw();
  } catch (err) {
    console.error("traceview error:", err);
    $$("#trace-meta").textContent = `trace error: ${err.message}`;
  }
}
