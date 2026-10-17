// Client-side beam-normalization orchestration.
//
// When the normalize toggle is on, the reaction control needs per-strip gains
// (gain[s] = 1 / mean(beam per-strip total)) computed from a strip = -1 beam
// reference run's events_MeV, injected as [detector] normalize_per_strip +
// strip_gains. This module:
//   1. builds the beam-reference control fields (same beam/gas; strip = -1),
//   2. runs it via POST /api/run and polls to done,
//   3. reads the beam run's events_MeV with JSROOT,
//   4. computes gain[s] = 1 / mean(sum(LeftdE[s],RightdE[s]) over beam events)
//      and returns them as an 18-element array.
//
// Mirrors the ~/MUSIC treatment (SimBeamGains per-strip beam gains) and the
// client-side jsroot projection — no uproot / backend change.

import { treeDraw } from "jsroot";
import { openResultRoot } from "./rootio.js";

const N_SEG = 16; // segmented strips 1..16

/**
 * Per-strip beam gains from a set of beam (unreacted, strip = -1) traces.
 * gain[s] = 1 / mean(beam per-strip total). `beamTraces` is an array of
 * per-strip arrays (18 values, strips 0..17). Mirrors the ~/MUSIC treatment so
 * the 2D trace view can normalize without a dedicated reference run.
 */
export function computeGains(beamTraces, nStrips = 18) {
  if (!beamTraces || !beamTraces.length) return null;
  const sums = new Array(nStrips).fill(0);
  for (const t of beamTraces) {
    for (let s = 0; s < nStrips && s < t.length; s++) sums[s] += t[s];
  }
  const gains = new Array(nStrips).fill(1.0);
  for (let s = 0; s < nStrips; s++) {
    const mean = sums[s] / beamTraces.length;
    if (mean > 0) gains[s] = 1.0 / mean;
  }
  return gains;
}

/** Apply per-strip gains to a trace (18 values, strips 0..17). */
export function applyGains(strip, gains) {
  return strip.map((v, s) => v * (gains && gains[s] ? gains[s] : 1));
}

/** POST /api/run and resolve the run_id immediately. */
export async function startRun(fields, normalize = false) {
  const res = await fetch("/api/run", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ fields, normalize, output: "beam.root" }),
  });
  if (!res.ok) throw new Error(`start beam run failed (${res.status})`);
  return (await res.json()).run_id;
}

/** Poll a run until done/error; reject on error. */
export async function waitRun(runId) {
  for (;;) {
    const res = await fetch(`/api/run/${runId}`);
    const d = await res.json();
    if (d.status === "done") return d;
    if (d.status === "error") throw new Error(d.error || "beam run error");
    await new Promise((r) => setTimeout(r, 800));
  }
}

/** Read the beam run's events_MeV and compute per-strip gains (18 values). */
export async function computeBeamGains(runId) {
  const root = await openResultRoot(runId);
  const ev = await root.readObject("events_MeV");
  const L = await treeDraw(ev, { expr: "LeftdE", graph: true, arr_limit: 1e9 });
  const R = await treeDraw(ev, { expr: "RightdE", graph: true, arr_limit: 1e9 });

  const w = N_SEG;
  const n = Math.floor(L.fY.length / w);
  // per-strip total accumulated over all beam events, segmented strip 1..16.
  const segTotal = new Array(N_SEG).fill(0);
  for (let e = 0; e < n; e++) {
    const b = e * w;
    for (let s = 0; s < w; s++) segTotal[s] += L.fY[b + s] + R.fY[b + s];
  }
  // gain[s] = 1 / mean(beam per-strip total). Strips that saw no beam (mean 0)
  // get gain 1.0 (leave unchanged). Index 0 and 17 are the full-width strips,
  // which are not in the L/R arrays; default them to 1.0.
  const gains = new Array(18).fill(1.0);
  for (let s = 0; s < N_SEG; s++) {
    const mean = segTotal[s] / Math.max(n, 1);
    gains[s + 1] = mean > 0 ? 1.0 / mean : 1.0;
  }
  return gains;
}

/**
 * Build the beam-reference fields (strip = -1, no reaction step) from the
 * reaction control's beam/gas/target. events_per_strip is small.
 */
export function buildBeamRefFields(control) {
  const c = JSON.parse(JSON.stringify(control)); // deep copy
  delete c.reaction;
  c.detector = { ...(c.detector || {}) };
  delete c.detector.strip_first;
  delete c.detector.strip_last;
  if ("strip" in c.detector) delete c.detector.strip;
  c.detector.strip = -1;
  // Beam reference needs a small number of events.
  c.run = { ...(c.run || {}), events_per_strip: 40, seed: 777 };
  return c;
}
