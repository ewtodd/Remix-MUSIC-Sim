// App shell: load the reaction/control list, select a control, render the
// curated form, and wire up the run button.

import { renderControls, state } from "./controls.js";
import { submitRun } from "./run.js";

const $$ = (sel) => document.querySelector(sel);

async function loadReactions() {
  const res = await fetch("/api/reactions");
  if (!res.ok) throw new Error(`GET /api/reactions failed (${res.status})`);
  const data = await res.json();
  return data.reactions || [];
}

async function loadControl(reactionId, name) {
  const res = await fetch(`/api/control/${encodeURIComponent(reactionId)}/${encodeURIComponent(name)}`);
  if (!res.ok) throw new Error(`GET /api/control failed (${res.status})`);
  return res.json();
}

function populateReactionSelect(reactions) {
  const sel = $$("#reaction-select");
  sel.innerHTML = "";
  for (const r of reactions) {
    const opt = document.createElement("option");
    opt.value = r.id;
    opt.textContent = r.label;
    sel.appendChild(opt);
  }
  sel.addEventListener("change", () => populateControlSelect(reactions));
  populateControlSelect(reactions);
}

function populateControlSelect(reactions) {
  const rsel = $$("#reaction-select");
  const reaction = reactions.find((r) => r.id === rsel.value);
  const csel = $$("#control-select");
  csel.innerHTML = "";
  if (!reaction) return;
  for (const c of reaction.controls) {
    const opt = document.createElement("option");
    opt.value = c.name;
    opt.textContent = c.label;
    csel.appendChild(opt);
  }
}

async function onLoadControl() {
  const reactionId = $$("#reaction-select").value;
  const name = $$("#control-select").value;
  if (!reactionId || !name) return;
  try {
    const info = await loadControl(reactionId, name);
    renderControls(info);
  } catch (err) {
    alert(`Failed to load control: ${err.message}`);
  }
}

async function boot() {
  try {
    const reactions = await loadReactions();
    populateReactionSelect(reactions);
  } catch (err) {
    alert(`Failed to load reactions: ${err.message}`);
  }

  // Collapsible control drawer — visuals stay front-and-center.
  $$("#toggle-drawer").addEventListener("click", () => {
    $$(".layout").classList.toggle("drawer-collapsed");
    $$("#toggle-drawer").textContent =
      $$(".layout").classList.contains("drawer-collapsed") ? "Show controls" : "Controls";
  });

  $$("#load-control-btn").addEventListener("click", onLoadControl);
  $$("#run-btn").addEventListener("click", () => submitRun());
  $$("#cancel-btn").addEventListener("click", () => {
    // The backend doesn't support cancelling yet; hide for now.
    // TODO(Phase 5/6): add a stop endpoint.
  });

  state.onDone = (runId) => {
    console.log("Run done", runId);
    // Phase 3: 2D trace view.
    import("./traceview.js").then((m) => m.showTrace(runId));
    // Phase 4: 3D path view.
    import("./pathview.js").then((m) => m.show3D(runId));
  };
}

boot();
