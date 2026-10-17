// Validate three.js detector + coarse MC path geometry builds without error (headless).
import * as THREE from "three";

// Mirrors Geometry.cpp: 20 anode rows, 2 columns. Strip pitch ~1.578 cm along z.
// strip 0 & 17 full-width; strips 1..16 split L/R. Height along y.
const ANODE_ROWS = 20;
const PITCH = 1.578;         // cm along z per row
const DET_WIDTH = 10.0;      // cm in x (full width)
const DET_HEIGHT = 10.0;     // cm in y
const HALF_W = DET_WIDTH / 2, HALF_H = DET_HEIGHT / 2;

function stripZ(i) { return (i - ANODE_ROWS / 2) * PITCH; }  // center z of row i

// Build one anode as a thin box (half-width in x for split strips).
function anodeBox(i, col) {
  const isFull = (i === 0 || i === 17);
  const w = isFull ? DET_WIDTH : DET_WIDTH / 2;           // x extent
  const x = isFull ? 0 : (col === 0 ? -HALF_W / 2 : HALF_W / 2);
  const z = stripZ(i);
  const geo = new THREE.BoxGeometry(w, DET_HEIGHT, PITCH * 0.95);
  const mesh = new THREE.Mesh(geo);
  mesh.position.set(x, 0, z);
  return mesh;
}

const scene = new THREE.Scene();
for (let i = 0; i < ANODE_ROWS; i++) {
  for (let c = 0; c < 2; c++) scene.add(anodeBox(i, c));
}

// Coarse MC path: beam enters along z at x=y=0, reacts at vertex, then
// evap + residue diverge to their stop points.
function pathLine(p0, p1, color) {
  const g = new THREE.BufferGeometry().setFromPoints([
    new THREE.Vector3(...p0), new THREE.Vector3(...p1)
  ]);
  return new THREE.Line(g, new THREE.LineBasicMaterial({ color }));
}

const vertex = [0.0, 0.0, stripZ(8)];
const beamStart = [0.0, 0.0, stripZ(0)];
const evapStop = [1.2, 0.3, stripZ(12)];
const residueStop = [-0.4, -0.2, stripZ(10)];

scene.add(pathLine(beamStart, vertex, 0x888888));   // beam (gray)
scene.add(pathLine(vertex, evapStop, 0x4444ff));    // evap particle
scene.add(pathLine(vertex, residueStop, 0xff4444)); // residue

let nAnode = 0, nLines = 0;
scene.traverse(o => { if (o.isMesh) nAnode++; if (o.isLine) nLines++; });
console.log("Detector anode meshes:", nAnode, "| path lines:", nLines);

// Validate positions / math
const first = scene.children[0].position;
console.log("row0 mesh x,y,z:", first.x, first.y, first.z.toFixed(2));
console.log("stripZ(8) =", stripZ(8).toFixed(2), "cm;  stripZ(17) =", stripZ(17).toFixed(2), "cm");

// Attempt an offscreen render (works if a GL context is available; otherwise report).
let rendered = false;
try {
  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(800, 600);
  const camera = new THREE.PerspectiveCamera(60, 800 / 600, 0.1, 1000);
  camera.position.set(20, 15, 25); camera.lookAt(0, 0, 0);
  scene.add(new THREE.AmbientLight(0xffffff, 1));
  renderer.render(scene, camera);
  rendered = true;
  console.log("Headless WebGL render: OK");
} catch (e) {
  console.log("Headless WebGL unavailable (expected on headless):", e.message.split("\n")[0]);
}
console.log("RESULT:", rendered ? "render-ok" : "geometry-math-ok (no GL)");
