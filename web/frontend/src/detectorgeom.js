// Detector geometry, mirroring the hardcoded mapping in `Geometry.cpp`
// (LoadHardcodedAnodeGeometry, 2019/11/04 anode geometry).
//
// 20 rows × 2 columns. Strip IDs 0 and 17 are single-column full-width strips;
// 1..16 are split into left/right halves; rows 0 and 19 are dead-layer caps.
// Units are centimetres. `+z` is the beam direction (downstream).

export const ANODE_ROWS = 20;
export const ANODE_COLS = 2;
export const N_STRIPS = 18; // readout strips 0..17

// Row definitions: { id, name, dxL, dxR, dz, color }.
// id = readout strip id (-1 upstream dead, -2 downstream dead). Colors are the
// ROOT kColor* values from the C++ (916/425/609/...) — used to tint the 3D mesh.
export const ROWS = (() => {
  const rows = [];
  const row = (id, name, dxL, dxR, dz, color) => ({ id, name, dxL, dxR, dz, color });
  // Upstream dead layer (full width).
  rows.push(row(-1, "DeadUS", 9, 9, 3.590, 920));
  // Strip 0 (full width).
  rows.push(row(0, "S0", 9, 9, 1.578, 416));
  // Strips 1..16 split L/R. Odd s → R short(4)/L long(5); even s → L short(4)/R long(5).
  for (let s = 1; s <= 16; s++) {
    const dxR = (s % 2 === 1) ? 4 : 5;
    const dxL = 9 - dxR;
    rows.push(row(s, `S${s}`, dxL, dxR, 1.578, 425));
  }
  // Strip 17 (full width).
  rows.push(row(17, "S17", 9, 9, 1.578, 416));
  // Downstream dead layer.
  rows.push(row(-2, "DeadDS", 9, 9, 3.522, 920));
  return rows;
})();

// Accumulate z at the START of each row (cm). Row 0 starts at z=0.
export const ROW_Z0 = (() => {
  const z = [];
  let acc = 0;
  for (const r of ROWS) {
    z.push(acc);
    acc += r.dz;
  }
  return z;
})();

// Total detector depth along beam axis.
export const ANODE_DEPTH = ROWS.reduce((a, r) => a + r.dz, 0);

// The x-centre of column `col` (0=R at -x, 1=L at +x) for a given row.
// x is centred: full active length per row = dxL + dxR (9 cm for split strips).
export function colX(rowIdx, col) {
  const r = ROWS[rowIdx];
  const total = r.dxL + r.dxR;
  const xc = -total / 2;
  return (col === 0) ? xc + r.dxR / 2 : xc + r.dxR + r.dxL / 2;
}

// Strip id -> z-centre (cm) of that strip's row.
export function stripZ(id) {
  for (let i = 0; i < ROWS.length; i++) {
    if (ROWS[i].id === id) return ROW_Z0[i] + ROWS[i].dz / 2;
  }
  return null;
}

// Strip id -> {z0, z1, rowIdx, dxL, dxR, height}.
export function stripInfo(id) {
  for (let i = 0; i < ROWS.length; i++) {
    if (ROWS[i].id === id) {
      const r = ROWS[i];
      return {
        id,
        name: r.name,
        rowIdx: i,
        z0: ROW_Z0[i],
        z1: ROW_Z0[i] + r.dz,
        zc: ROW_Z0[i] + r.dz / 2,
        dxL: r.dxL,
        dxR: r.dxR,
        height: 10,
      };
    }
  }
  return null;
}
