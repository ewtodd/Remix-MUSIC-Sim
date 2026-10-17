// Open a served .root file with JSROOT from a raw byte buffer.
//
// jsroot's openFile() accepts an ArrayBuffer directly (no blob URL), which is
// the robust path: it avoids object-URL lifecycle/range-request failures
// ("can't access property getUint8, this.arr is null", ERR_REQUEST_RANGE_
// NOT_SATISFIABLE) that happen when jsroot lazily fetches ranges from a blob
// that was already revoked or that a CDN/browser serves without ranges.
//
// Use this everywhere a .root is fetched from /api.
import { openFile } from "jsroot";

export async function openResultRoot(runId) {
  const res = await fetch(`/api/run/${runId}/result.root`);
  if (!res.ok) throw new Error(`result.root ${res.status}`);
  const buf = await res.arrayBuffer();
  return openFile(buf);
}
