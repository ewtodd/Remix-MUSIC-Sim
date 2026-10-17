import puppeteer from "puppeteer-core";
const CHROME = "/nix/store/ky199ji3a5drs6rm3wr5zidamx9dhz56-chromium-155.0.8059.39/bin/chromium";
const browser = await puppeteer.launch({
  executablePath: CHROME,
  headless: "new",
  args: ["--no-sandbox", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"],
});
const msgs = [];
const page = await browser.newPage();
await page.setViewport({ width: 1280, height: 900 });
page.on("console", (m) => msgs.push(`[${m.type()}] ${m.text()}`));
page.on("pageerror", (e) => msgs.push(`[pageerror] ${e.message}`));

await page.goto("http://127.0.0.1:5180/", { waitUntil: "networkidle2", timeout: 25000 });

await page.select("#reaction-select", "37Cl_alpha_n");
await page.select("#control-select", "37Cl_alpha_n_bulk.toml");
await page.click("#load-control-btn");
await page.waitForFunction(() => document.querySelector('input[data-key="beam.energy"]')?.value, { timeout: 8000 });
console.log("control loaded");

// Throttle events for speed.
await page.$eval('input[data-key="run.events_per_strip"]', (el) => {
  el.value = "40"; el.dispatchEvent(new Event("input", { bubbles: true }));
});

await page.click("#run-btn");
await page.waitForFunction(() => {
  const s = document.querySelector("#run-status")?.textContent || "";
  return s.includes("done");
}, { timeout: 90000 });
console.log("RUN done:", await page.$eval("#run-status", (el) => el.textContent));

// 3D view should appear.
try {
  await page.waitForFunction(() => {
    const r = document.querySelector("#results3d");
    return r && !r.hidden && (document.querySelector("#ev3d-meta")?.textContent || "").includes("n=");
  }, { timeout: 60000 });
} catch (e) {
  console.log("3D WAIT TIMEOUT.");
  console.log("  results3d hidden:", await page.$eval("#results3d", (el) => el.hidden).catch(()=>"(none)"));
  console.log("  ev3d-meta:", await page.$eval("#ev3d-meta", (el) => el.textContent).catch(()=>"(none)"));
  console.log("  path3d-container childCount:", await page.$eval("#path3d-container", (el) => el.childElementCount).catch(()=>"(none)"));
  for (const m of msgs) console.log(m);
  await browser.close(); process.exit(1);
}
const meta = await page.$eval("#ev3d-meta", (el) => el.textContent);
console.log("3D META:", meta);
const sliderMax = await page.$eval("#ev3d-slider", (el) => el.max);
const productCount = await page.$$eval("#ev3d-products input", (els) => els.length);
console.log("3D slider max:", sliderMax, "| product toggles:", productCount);

// Non-background pixel count in the 3D canvas.
async function nonBg() {
  return page.evaluate(() => {
    const c = document.querySelector("#path3d-container canvas");
    if (!c) return -1;
    const gl = c.getContext("webgl2") || c.getContext("webgl");
    if (!gl) return -2;
    const w = c.width, h = c.height;
    const px = new Uint8Array(w * h * 4);
    gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, px);
    let n = 0;
    for (let i = 0; i < px.length; i += 4) {
      // background is ~#0f1420 (15,20,32); count clearly different pixels
      if (Math.abs(px[i]-15) > 30 || Math.abs(px[i+1]-20) > 30 || Math.abs(px[i+2]-32) > 30) n++;
    }
    return n;
  });
}
console.log("3D non-bg pixels:", await nonBg());

await page.screenshot({ path: "path3d_view.png", fullPage: true });
for (const m of msgs) console.log(m);
console.log("DONE");
await browser.close();
