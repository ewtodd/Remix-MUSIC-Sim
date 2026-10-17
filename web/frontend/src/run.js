// Run flow: POST /api/run, then poll GET /api/run/{id} until done, surfacing
// status + log tail. When done, the client is ready to fetch the result .root
// (used by the Phase 3 trace view).

import { state, controlsToJson } from "./controls.js";
import { computeBeamGains, buildBeamRefFields, startRun, waitRun } from "./normalize.js";

const $$ = (sel) => document.querySelector(sel);

let pollTimer = null;

export async function submitRun() {
  const body = {
    fields: controlsToJson(),
    normalize: state.normalize,
    output: "result.root",
  };

  // If normalization is on, first run a strip = -1 beam reference to derive the
  // per-strip gains, then inject them into the reaction control. This keeps the
  // whole normalization feature client-side (jsroot) — the backend stays thin.
  if (state.normalize) {
    setStatus("Normalizing: running beam reference…", "pending");
    try {
      const beamFields = buildBeamRefFields(body.fields);
      const beamRunId = await startRun(beamFields);
      await waitRun(beamRunId);
      const gains = await computeBeamGains(beamRunId);
      body.fields = body.fields || {};
      body.fields.detector = body.fields.detector || {};
      body.fields.detector.normalize_per_strip = true;
      body.fields.detector.strip_gains = gains;
      setStatus("Normalizing: beam reference done — running reaction…", "pending");
    } catch (err) {
      setStatus(`Normalization failed: ${err.message}`, "error");
      return null;
    }
  }

  const statusEl = $$("#run-status");
  setStatus("Starting…", "pending");

  let runId;
  try {
    const res = await fetch("/api/run", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    if (!res.ok) {
      const err = await res.text();
      throw new Error(`POST /api/run failed (${res.status}): ${err}`);
    }
    const data = await res.json();
    runId = data.run_id;
  } catch (err) {
    setStatus(String(err.message || err), "error");
    return null;
  }

  setStatus(`Run ${runId}: running…`, "running");
  $$("#cancel-btn").hidden = false;
  showLog(true);

  pollTimer = setInterval(() => pollRun(runId), 1000);
  return runId;
}

async function pollRun(runId) {
  try {
    const res = await fetch(`/api/run/${runId}`);
    if (!res.ok) throw new Error(`status ${res.status}`);
    const data = await res.json();

    if (data.log_tail) {
      const log = $$("#run-log");
      log.textContent = data.log_tail;
      log.scrollTop = log.scrollHeight;
    }

    if (data.status === "done") {
      clearInterval(pollTimer);
      setStatus(`Run ${runId}: done`, "done");
      $$("#cancel-btn").hidden = true;
      state.onDone && state.onDone(runId);
    } else if (data.status === "error") {
      clearInterval(pollTimer);
      setStatus(`Run ${runId}: error — ${data.error || "see log"}`, "error");
      $$("#cancel-btn").hidden = true;
    } else {
      // running — phase may be "cache" while the one-time tables are built.
      const phase = data.phase === "cache" ? " (generating tables…)" : "";
      setStatus(`Run ${runId}: ${data.status}${phase}`, "running");
    }
  } catch (err) {
    clearInterval(pollTimer);
    setStatus(`Polling error: ${err.message}`, "error");
  }
}

function setStatus(text, kind) {
  const el = $$("#run-status");
  el.textContent = text;
  el.className = "status " + kind;
}

function showLog(show) {
  $$("#run-log").hidden = !show;
}
