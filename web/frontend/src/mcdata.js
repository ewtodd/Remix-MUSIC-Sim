// Client-side projection of the MC tree into per-event records for the 3D path
// view. Uses JSROOT treeDraw + graph:true to pull branch values, then zips them
// by entry index. No uproot / backend round-trip.

import { treeDraw } from "jsroot";
import { openResultRoot } from "./rootio.js";

const MC_BRANCHES = [
  "n_steps", "event_index", "reaction_strip", "reaction_channel",
  "vertex_x", "vertex_y", "vertex_z",
  "beam_stop_x", "beam_stop_y", "beam_stop_z", "beam_termination",
  "evap_stop_x", "evap_stop_y", "evap_stop_z", "evap_termination",
  "residue_stop_x", "residue_stop_y", "residue_stop_z", "residue_stop_strip", "residue_termination",
  "beam_energy_gas", "beam_energy_exit",
];

export async function loadMC(runId) {
  const root = await openResultRoot(runId);
  const mc = await root.readObject("MC");
  const data = {};
  for (const name of MC_BRANCHES) {
    const res = await treeDraw(mc, { expr: name, graph: true, arr_limit: 10000000, want_hist: false });
    data[name] = res.fY;
  }
  const N = (data.n_steps || []).length;
  const events = [];
  for (let i = 0; i < N; i++) {
    const ev = { index: i };
    for (const name of MC_BRANCHES) ev[name] = data[name][i];
    events.push(ev);
  }
  return events;
}
