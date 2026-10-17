// 3D path view controller. Loads MC records (jsroot), builds the three.js
// scene, and wires the event slider, product checkboxes, and overlay toggle.

import { loadMC } from "./mcdata.js";
import { setProducts, setMode, mount } from "./path3d.js";

const $$ = (sel) => document.querySelector(sel);

export async function show3D(runId) {
  const panel = $$("#results3d");
  // Reveal the 3D body; front-and-center.
  panel.querySelector(".vis-placeholder").hidden = true;
  panel.querySelector(".vis-body").hidden = false;
  const container = $$("#path3d-container");
  const meta = $$("#ev3d-meta");
  const label = $$("#ev3d-label");
  const slider = $$("#ev3d-slider");

  // Product toggles.
  const toggles = $$("#ev3d-products");
  toggles.innerHTML = "";
  const active = new Set(["beam", "evap", "residue"]);
  for (const [id, color, labelText] of [["beam", "#ffffff", "Beam"], ["evap", "#00ffff", "Evap"], ["residue", "#ff3366", "Residue"]]) {
    const lab = document.createElement("label");
    lab.className = "inline product";
    lab.innerHTML = `<input type="checkbox" data-product="${id}" checked /><span style="color:${color}">${labelText}</span>`;
    lab.querySelector("input").addEventListener("change", (e) => {
      e.target.checked ? active.add(id) : active.delete(id);
      setProducts([...active]);
      redraw();
    });
    toggles.appendChild(lab);
  }

  let events = [];
  let current = 0;

  function redraw() {
    setMode($$("#ev3d-overlay").checked ? "overlay" : "single");
    setProducts([...active]);
    mount(container, events, current);
    label.textContent = $$("#ev3d-overlay").checked ? "all" : String(current);
  }

  slider.addEventListener("input", () => {
    current = Number(slider.value);
    label.textContent = String(current);
    redraw();
  });
  $$("#ev3d-overlay").addEventListener("change", () => {
    slider.disabled = $$("#ev3d-overlay").checked;
    redraw();
  });

  meta.textContent = "loading MC records…";
  try {
    events = await loadMC(runId);
    slider.max = String(Math.max(0, events.length - 1));
    const reacted = events.filter((e) => e.beam_termination === 1).length;
    meta.textContent = `n=${events.length} (${reacted} reacted)`;
    redraw();
  } catch (err) {
    console.error("3D view error:", err);
    meta.textContent = `3D error: ${err.message}`;
  }
}
