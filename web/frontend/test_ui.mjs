import puppeteer from "puppeteer-core";

const CHROME = "/nix/store/ky199ji3a5drs6rm3wr5zidamx9dhz56-chromium-155.0.8059.39/bin/chromium";
const URL = "http://127.0.0.1:5180/";

const browser = await puppeteer.launch({
  executablePath: CHROME,
  headless: "new",
  args: ["--no-sandbox", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"],
});
const page = await browser.newPage();
await page.setViewport({ width: 1280, height: 900 });

const consoleMsgs = [];
page.on("console", (m) => consoleMsgs.push(`[${m.type()}] ${m.text()}`));
page.on("pageerror", (e) => consoleMsgs.push(`[pageerror] ${e.message}`));

await page.goto(URL, { waitUntil: "networkidle2", timeout: 20000 });

// 1. reaction dropdown populated?
const reactionOptions = await page.$$eval("#reaction-select option", (o) => o.map((x) => x.value));
console.log("reaction options:", reactionOptions.length);

// 2. select the 37Cl_alpha_n reaction
await page.select("#reaction-select", "37Cl_alpha_n");
const controlOptions = await page.$$eval("#control-select option", (o) => o.map((x) => x.value));
console.log("control options:", controlOptions);

// 3. load a control
await page.select("#control-select", "37Cl_alpha_n_bulk.toml");
await page.click("#load-control-btn");
await page.waitForFunction(() => {
  const v = document.querySelector('input[data-key="beam.energy"]')?.value;
  return v && v !== "";
}, { timeout: 5000 });
const beamEnergy = await page.$eval('input[data-key="beam.energy"]', (el) => el.value);
const gasSpecies = await page.$eval('input[data-key="gas.species"]', (el) => el.value);
console.log("curated form loaded: beam.energy =", beamEnergy, "| gas.species =", gasSpecies);
const rawLen = await page.$eval("#raw-toml", (el) => el.value.length);
console.log("raw toml length:", rawLen);

// 4. Set a small event count for speed, then run.
await page.$eval('input[data-key="run.events_per_strip"]', (el) => {
  el.value = "20";
  el.dispatchEvent(new Event("input"));
});
await page.click("#run-btn");
await page.waitForFunction(() => {
  const s = document.querySelector("#run-status")?.textContent || "";
  return s.includes("done") || s.includes("error");
}, { timeout: 90000 });
const statusText = await page.$eval("#run-status", (el) => el.textContent);
console.log("run status:", statusText);

console.log("--- console messages ---");
for (const m of consoleMsgs) console.log(m);

await browser.close();
console.log("DONE");
