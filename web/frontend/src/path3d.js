// 3D path view using three.js. Renders the MUSIC anode stack (from the
// detector geometry mapping) plus the coarse per-event particle paths:
//   beam:      (0,0,0) -> vertex (reacted) OR -> beam_stop (unreacted)
//   evap:      vertex -> evap_stop
//   residue:   vertex -> residue_stop
// Supports single-event and multi-event overlay modes, product toggles, and
// orbit controls.

import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { ROWS, ROW_Z0, N_STRIPS } from "./detectorgeom.js";

const PRODUCTS = [
  { id: "beam",    color: 0xffffff, label: "Beam" },
  { id: "evap",    color: 0x00ffff, label: "Evap" },
  { id: "residue", color: 0xff3366, label: "Residue" },
];
const PRODUCT_COLOR = Object.fromEntries(PRODUCTS.map(p => [p.id, p.color]));

let renderer = null, scene = null, camera = null, controls = null;
let trafficState = { products: ["beam","evap","residue"], mode: "single" };
let renderLoop = false;

function makeRenderer(container) {
  if (renderer) return;
  // Cap the pixel ratio: full devicePixelRatio (often 2-3x on HiDPI) massively
  // inflates the fragment cost of the translucent detector stack and hundreds
  // of overlay tubes. 1.5 is a good balance of crispness and drag speed.
  const px = Math.min(window.devicePixelRatio || 1, 1.5);
  renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: "high-performance" });
  renderer.setPixelRatio(px);
  renderer.setSize(container.clientWidth, container.clientHeight);
  container.appendChild(renderer.domElement);
  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x0f1420);
  camera = new THREE.PerspectiveCamera(55, container.clientWidth / container.clientHeight, 0.1, 2000);
  camera.position.set(40, 30, 55);
  camera.lookAt(0, 0, 16);
  controls = new OrbitControls(camera, renderer.domElement);
  controls.target.set(0, 0, 16);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.update();
  // OrbitControls dispatches a 'change' event whenever the camera moves (incl.
  // each damped frame). We re-render only then, so an idle scene costs nothing.
  controls.addEventListener("change", () => { needsRender = true; });
  scene.add(new THREE.AmbientLight(0xffffff, 0.85));
  scene.add(new THREE.DirectionalLight(0xffffff, 0.6));
  startRenderLoop();
}

let needsRender = true;

// Continuous render loop. Damping requires controls.update() every frame (the
// camera eases toward its target over several frames). We still only run the
// expensive renderer.render() when a 'change' event flagged a real camera move,
// so an idle scene does no GPU work and dragging stays responsive.
function startRenderLoop() {
  if (renderLoop) return;
  renderLoop = true;
  const loop = () => {
    requestAnimationFrame(loop);
    controls.update(); // applies damping; fires 'change' only while moving
    if (needsRender) {
      needsRender = false;
      renderer.render(scene, camera);
    }
  };
  loop();
}

// Build the detector stack (translucent anode boxes + dead caps + axes).
function buildDetector(group) {
  const height = 10;
  for (let i = 0; i < ROWS.length; i++) {
    const r = ROWS[i];
    const z0 = ROW_Z0[i];
    const isDead = r.id < 0;
    const opacity = isDead ? 0.25 : 0.16;
    const mat = new THREE.MeshLambertMaterial({
      color: r.color, transparent: true, opacity, side: THREE.DoubleSide, depthWrite: false,
    });
    // Two column boxes (L/R) if the strip is split; otherwise one full box.
    const cols = (r.id >= 1 && r.id <= 16) ? 2 : 1;
    for (let c = 0; c < cols; c++) {
      const dx = (cols === 2) ? (c === 0 ? r.dxR : r.dxL) : r.dxL;
      const geo = new THREE.BoxGeometry(dx, height, r.dz);
      const mesh = new THREE.Mesh(geo, mat);
      const total = r.dxL + r.dxR;
      const xc = -total / 2;
      const x = (cols === 2) ? (c === 0 ? xc + r.dxR / 2 : xc + r.dxR + r.dxL / 2) : xc + total / 2;
      mesh.position.set(x, 0, z0 + r.dz / 2);
      group.add(mesh);
    }
    // Edges for a crisp grid.
    const edge = new THREE.BoxGeometry(r.dxL + r.dxR, height, r.dz);
    const edgeMat = new THREE.LineBasicMaterial({ color: 0x2c3852 });
    const eg = new THREE.LineSegments(new THREE.EdgesGeometry(edge), edgeMat);
    const xc = -(r.dxL + r.dxR) / 2;
    eg.position.set(xc + (r.dxL + r.dxR) / 2, 0, z0 + r.dz / 2);
    group.add(eg);
  }
  // Axes.
  const axis = (len, color, dir) => {
    const geo = new THREE.CylinderGeometry(0.2, 0.2, len, 8);
    const m = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ color }));
    m.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0), dir);
    m.position.copy(dir.clone().multiplyScalar(len/2));
    group.add(m);
  };
  axis(16, 0xff5555, new THREE.Vector3(1,0,0));   // +x red
  axis(16, 0x55ff55, new THREE.Vector3(0,1,0));   // +y green
  axis(50, 0x5588ff, new THREE.Vector3(0,0,1));   // +z blue (beam dir)
}

function tube(group, p0, p1, color, radius = 0.5, opacity = 0.95) {
  const a = new THREE.Vector3(...p0), b = new THREE.Vector3(...p1);
  const dir = b.clone().sub(a);
  const len = dir.length();
  if (len < 1e-4) return;
  const geo = new THREE.CylinderGeometry(radius, radius, len, 12);
  const mat = new THREE.MeshBasicMaterial({ color, transparent: true, opacity, depthTest: false });
  const mesh = new THREE.Mesh(geo, mat);
  mesh.renderOrder = 10;
  mesh.position.copy(a.clone().add(b).multiplyScalar(0.5));
  mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0), dir.normalize());
  return mesh;
}

// Build the coarse path for one event from its MC record.
export function eventPath(ev) {
  const vx = ev.vertex_x, vy = ev.vertex_y, vz = ev.vertex_z;
  const paths = [];
  const beamStop = [ev.beam_stop_x, ev.beam_stop_y, ev.beam_stop_z];
  const reacted = ev.beam_termination === 1; // absorbed at the vertex
  const vertex = [vx, vy, vz];
  // Beam: entrance (0,0,0) -> vertex (reacted) or -> beam_stop (unreacted).
  if (reacted) {
    paths.push({ product: "beam", pts: [[0,0,0], vertex] });
    paths.push({ product: "evap", pts: [vertex, [ev.evap_stop_x, ev.evap_stop_y, ev.evap_stop_z]] });
    paths.push({ product: "residue", pts: [vertex, [ev.residue_stop_x, ev.residue_stop_y, ev.residue_stop_z]] });
  } else {
    paths.push({ product: "beam", pts: [[0,0,0], beamStop] });
  }
  return { vertex, reacted, paths };
}

export function setProducts(products) { trafficState.products = products; }
export function setMode(mode) { trafficState.mode = mode; }

export function show3D(mcEvents, container, current = 0) {
  makeRenderer(container);
  // Re-append the renderer canvas if it was detached (e.g. container cleared).
  if (!renderer.domElement.parentNode) container.appendChild(renderer.domElement);
  // Reset scene children (keep background).
  while (scene.children.length) scene.remove(scene.children[0]);
  buildDetector(scene);

  const tracks = new THREE.Group();
  scene.add(tracks);

  const drawOne = (ev) => {
    const { paths } = eventPath(ev);
    for (const p of paths) {
      if (!trafficState.products.includes(p.product)) continue;
      const t = tube(tracks, p.pts[0], p.pts[1], PRODUCT_COLOR[p.product]);
      if (t) tracks.add(t);
    }
  };

  if (trafficState.mode === "overlay") {
    mcEvents.forEach(drawOne);
  } else if (mcEvents[current]) {
    drawOne(mcEvents[current]);
  }

  renderer.setSize(container.clientWidth, container.clientHeight);
  // Signal the render loop to draw the (re)built scene.
  needsRender = true;
}

// Load MC records and hand them to the viewer.
export function mount(container, mcEvents, current = 0) {
  // Don't wipe the container's canvas — we re-append it in show3D if needed.
  show3D(mcEvents, container, current);
}

export function dispose() {
  if (renderer) {
    renderer.dispose();
    renderer = null; scene = null; camera = null; controls = null;
  }
}
