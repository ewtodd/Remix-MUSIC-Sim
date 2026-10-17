import { openFile } from "jsroot";

const rootPath = process.argv[2] || "/home/e-work/Remix-MUSIC-Sim/web/runs/_spike_sim/out/spike.root";
const f = await openFile(rootPath);
console.log("opened:", rootPath);

// Top-level keys
const top = f.fKeys || [];
console.log("top-level keys:", top.map(k => `${k.fName} (${k.fclassName})`).join(", "));

async function dirKeys(name) {
  // A TDirectoryFile is read back as an object; use readObject to get it
  const dir = await f.readObject(name);
  if (!dir) return [];
  const kk = dir.fKeys || [];
  return kk.map(k => `${k.fName} (${k.fclassName})`);
}

// MC tree
try {
  const mc = await f.readObject("MC");
  console.log("\n=== MC tree ===");
  console.log("typename:", mc?._typename);
  if (mc) {
    const br = mc.fBranches || [];
    console.log("branches:", br.map(b => b.fName).join(", "));
  }
} catch (e) {
  console.log("MC readObject error:", e.message);
}

// traces directory
try {
  console.log("\n=== traces/ directory ===");
  const tk = await dirKeys("traces");
  console.log("trace keys count:", tk.length);
  console.log("first 8:", tk.slice(0, 8).join(", "));
} catch (e) {
  console.log("traces readObject error:", e.message);
}

// metadata directory
try {
  console.log("\n=== metadata/ directory ===");
  const mk = await dirKeys("metadata");
  console.log("metadata keys count:", mk.length);
  console.log("first 10:", mk.slice(0, 10).join(", "));
} catch (e) {
  console.log("metadata readObject error:", e.message);
}
