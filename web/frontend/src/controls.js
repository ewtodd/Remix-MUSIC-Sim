// Curated control form + raw TOML editor, kept in sync, plus the normalize
// toggle. The curated fields are a subset of the schema; the raw TOML covers
// everything else. When the user edits either side, the other follows.

/**
 * Update the raw TOML textarea from the current curated state is deferred:
 * editing raw TOML is authoritative. Editing a curated field re-serializes from
 * the curated dict. We keep a single source of truth (`state.control`), a dict
 * that mirrors the control schema.
 */

export const state = {
  control: null,      // parsed control dict (from /api/control)
  normalize: false,   // per-strip normalization toggle
};

const $$ = (sel) => document.querySelector(sel);

/** Render the curated form + raw TOML editor into #controls. */
export function renderControls(controlInfo) {
  // controlInfo: { toml, curated } from GET /api/control.
  state.control = parseToml(controlInfo.toml);
  const root = $$("#controls");
  root.innerHTML = "";

  // ---- Curated form ----
  const form = document.createElement("div");
  form.className = "curated";

  const curated = controlInfo.curated || {};
  form.appendChild(fieldset("Gas", gasFields(curated.gas)));
  form.appendChild(fieldset("Beam", beamFields(curated.beam)));
  form.appendChild(fieldset("Target", targetFields(curated.target)));
  form.appendChild(fieldset("Detector", detectorFields(curated.detector)));
  form.appendChild(fieldset("Run", runFields(curated.run)));

  // Normalize toggle.
  const norm = document.createElement("label");
  norm.className = "toggle";
  norm.innerHTML = `<input type="checkbox" id="normalize-toggle" />
    <span>Normalize traces to 1 a.u. per strip (beam reference run)</span>`;
  norm.querySelector("input").checked = state.normalize;
  norm.querySelector("input").addEventListener("change", (e) => {
    state.normalize = e.target.checked;
  });
  form.appendChild(norm);

  root.appendChild(form);

  // ---- Raw TOML editor ----
  const raw = document.createElement("div");
  raw.className = "raw";
  const rawLabel = document.createElement("h3");
  rawLabel.textContent = "Raw TOML (authoritative on edit)";
  const ta = document.createElement("textarea");
  ta.id = "raw-toml";
  ta.spellcheck = false;
  ta.value = dumpToml(state.control);
  ta.addEventListener("input", (e) => {
    try {
      state.control = parseToml(e.target.value);
      markRawValid(ta, true);
      state.onChange && state.onChange();
    } catch (err) {
      markRawValid(ta, false, err.message);
    }
  });
  raw.appendChild(rawLabel);
  raw.appendChild(ta);
  root.appendChild(raw);

  // Wire curated field edits to re-serialize the raw editor.
  root.addEventListener("input", (e) => {
    const inp = e.target;
    if (inp.closest(".curated") && inp.dataset.key) {
      setPath(state.control, inp.dataset.key, coerce(inp));
      ta.value = dumpToml(state.control);
      markRawValid(ta, true);
      state.onChange && state.onChange();
    }
  });
}

export function controlsToJson() {
  // The curated form edits `state.control` which already includes the overlay
  // shape. Send the whole thing as fields for the run.
  return state.control;
}

// ---------------------------------------------------------------------------
// Field builders
// ---------------------------------------------------------------------------

function fieldset(title, fields) {
  const fs = document.createElement("div");
  fs.className = "fieldset";
  const h = document.createElement("h3");
  h.textContent = title;
  fs.appendChild(h);
  for (const f of fields) {
    fs.appendChild(inputRow(f));
  }
  return fs;
}

function inputRow(f) {
  const label = document.createElement("label");
  label.className = "field";
  label.dataset.key = f.key;
  const span = document.createElement("span");
  span.textContent = f.label;
  let inp;
  if (f.type === "number") {
    inp = document.createElement("input");
    inp.type = "number";
    inp.step = "any";
    inp.value = f.value ?? "";
  } else {
    inp = document.createElement("input");
    inp.type = "text";
    inp.value = f.value ?? "";
  }
  inp.addEventListener("input", () => {
    label.dataset.key && (label.value = inp.value);
  });
  label.appendChild(span);
  label.appendChild(inp);
  // store the key on the input so the delegated listener can read it
  inp.dataset.key = f.key;
  inp.dataset.value = f.value ?? "";
  return label;
}

function gasFields(g) {
  g = g || {};
  return [
    { key: "gas.species", label: "Species", value: g.species },
    { key: "gas.pressure", label: "Pressure (Torr)", type: "number", value: g.pressure },
    { key: "gas.temperature", label: "Temperature (K)", type: "number", value: g.temperature },
  ];
}

function beamFields(b) {
  b = b || {};
  return [
    { key: "beam.species", label: "Species", value: b.species },
    { key: "beam.energy", label: "Energy (MeV)", type: "number", value: b.energy },
    { key: "beam.energy_fwhm", label: "Energy FWHM (MeV)", type: "number", value: b.energy_fwhm },
    { key: "beam.dedx_scale", label: "dE/dx scale", type: "number", value: b.dedx_scale },
  ];
}

function targetFields(t) {
  t = t || {};
  return [
    { key: "target.species", label: "Species", value: t.species },
    { key: "target.compound", label: "Compound", value: t.compound },
  ];
}

function detectorFields(d) {
  d = d || {};
  return [
    { key: "detector.eloss_bins", label: "eloss bins", type: "number", value: d.eloss_bins },
    { key: "detector.max_eloss", label: "max eloss (MeV)", type: "number", value: d.max_eloss },
    { key: "detector.strip", label: "strip (-1 = unreacted)", type: "number", value: d.strip },
    { key: "detector.strip_first", label: "strip_first", type: "number", value: d.strip_first },
    { key: "detector.strip_last", label: "strip_last", type: "number", value: d.strip_last },
  ];
}

function runFields(r) {
  r = r || {};
  return [
    { key: "run.events_per_strip", label: "events / strip", type: "number", value: r.events_per_strip },
    { key: "run.seed", label: "seed", type: "number", value: r.seed },
    { key: "run.threads", label: "threads", type: "number", value: r.threads },
  ];
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function setPath(obj, dotted, value) {
  const parts = dotted.split(".");
  let cur = obj;
  for (let i = 0; i < parts.length - 1; i++) {
    if (!cur[parts[i]] || typeof cur[parts[i]] !== "object") cur[parts[i]] = {};
    cur = cur[parts[i]];
  }
  cur[parts[parts.length - 1]] = value;
}

function coerce(inp) {
  if (inp.type === "number") {
    const v = Number(inp.value);
    return Number.isFinite(v) ? v : 0;
  }
  return inp.value;
}

function markRawValid(ta, ok, msg) {
  ta.classList.toggle("ok", ok);
  ta.classList.toggle("err", !ok);
  ta.title = ok ? "" : (msg || "invalid TOML");
}

// ---- TOML serialize/parse (small subset) ----
// The backend uses tomli_w; here we do a minimal, schema-shaped serializer so
// the browser shows a faithful editable TOML. Full fidelity is not required for
// the editor's round-trip because the backend re-parses the authored text.
import { dumpToml, parseToml } from "./toml.js";
