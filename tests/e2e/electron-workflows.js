// FLASHApp scenario for electron/scripts/e2e.js: runs FLASHDeconv and FLASHTnT
// (DecoyDatabase -> FLASHDeconv -> FLASHTnT -> FLASHAnnotator) on the bundled
// example data inside the Electron window, then opens each Viewer and checks
// that the Vue component draws its plots and tables. FLASHQuant is view-only,
// so its Viewer is checked on the example workspace.
const RUN_TIMEOUT = 45 * 60 * 1000;
const VIEW_TIMEOUT = 3 * 60 * 1000;

async function runWorkflow(t, workflowPage, inputs) {
  await t.open(workflowPage);
  await t.clickTab("Configure");
  for (const [label, file] of inputs) await t.choose(label, file);
  await t.shot(`${workflowPage}-configured`);
  await t.clickTab("Run");
  await t.clickButton("Start Workflow");
  await t.waitForText("Workflow completed successfully.", {
    timeout: RUN_TIMEOUT,
    failOn: ["Errors occurred, check log file."],
  });
  await t.shot(`${workflowPage}-finished`);
}

async function checkViewer(t, viewerPage, experiment, extraSelector) {
  await t.open(viewerPage);
  if (experiment) await t.choose("choose experiment", experiment);
  await t.waitForComponent(".js-plotly-plot, .tabulator");
  // The grid draws its cells one by one (WebGL in software can be slow), so
  // count across all component frames until plots and tables are both there.
  const count = async (sel) => {
    let n = 0;
    for (const f of t.componentFrames()) n += await f.locator(sel).count().catch(() => 0);
    return n;
  };
  const text = async () => (await Promise.all(t.componentFrames().map((f) => f.locator("body").innerText().catch(() => "")))).join("\n");
  const started = Date.now();
  let plots = 0;
  let tables = 0;
  while (Date.now() - started < VIEW_TIMEOUT) {
    await t.assertNoException(viewerPage);
    if ((await text()).includes("WebGL is not supported")) throw new Error(`${viewerPage}: WebGL plots did not render`);
    plots = await count(".js-plotly-plot");
    tables = await count(".tabulator");
    if (plots > 0 && tables > 0 && (!extraSelector || (await count(extraSelector)) > 0)) break;
    await t.page.waitForTimeout(2000);
  }
  t.log(`${viewerPage}: ${plots} plots, ${tables} tables in ${t.componentFrames().length} component frame(s)`);
  await t.shot(`${viewerPage}-${experiment || "default"}`);
  if (plots === 0 || tables === 0) {
    t.log((await text()).slice(0, 2000));
    throw new Error(`${viewerPage}: expected plots and tables, got ${plots} plots and ${tables} tables`);
  }
  if (extraSelector && (await count(extraSelector)) === 0) throw new Error(`${viewerPage}: ${extraSelector} missing`);
}

module.exports = async (t) => {
  await runWorkflow(t, "FLASHDeconvWorkflow", [["mzML-files", "example_fd.mzML"]]);
  await checkViewer(t, "FLASHDeconvViewer", "example_fd");

  await runWorkflow(t, "FLASHTnTWorkflow", [["mzML-files", "GapA_5458mz_HCD.mzML"]]);
  const log = await t.page.locator("body").innerText();
  if (!log.includes("FLASHAnnotator")) throw new Error("FLASHTnT run did not reach FLASHAnnotator");
  await checkViewer(t, "FLASHTnTViewer", "GapA_5458mz_HCD");

  await checkViewer(t, "FLASHQuantViewer", null);
};
