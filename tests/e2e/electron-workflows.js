// FLASHApp scenario for electron/scripts/e2e.js: runs the workflow twice on the
// bundled example data inside the Electron window, deconvolution only
// (FLASHDeconv) and with identification (DecoyDatabase -> FLASHDeconv ->
// FLASHTnT -> FLASHAnnotator), then opens the Viewer on each run and checks
// that the Vue component draws its plots and tables.
const RUN_TIMEOUT = 45 * 60 * 1000;
const VIEW_TIMEOUT = 3 * 60 * 1000;

async function setCheckbox(t, label, checked) {
  const box = t.page.getByRole("checkbox", { name: label }).first();
  if ((await box.isChecked()) !== checked) {
    // Streamlit draws the input at 0x0 behind its styled box, so even a forced
    // click on it fails ("outside of the viewport"); click its label instead.
    await box.locator("xpath=ancestor::label[1]").click();
    await t.settle();
  }
}

async function runWorkflow(t, workflowPage, inputs, identify) {
  // One page per step: Upload data, Describe experiment, Run analysis
  await t.open(`${workflowPage}Configure`);
  for (const [label, file] of inputs) await t.choose(label, file);
  await setCheckbox(t, "Identify proteoforms", identify);
  await t.shot(`${workflowPage}-${identify ? "identify" : "deconvolve"}-configured`);
  await t.open(`${workflowPage}Run`);
  // No settle() here: while the workflow runs, the status fragment reruns the
  // page, which can keep it from going quiet for 120 s. waitForText waits.
  await t.page.getByRole("button", { name: "Start Workflow" }).first().click();
  // While the new run is running, its page reruns every second and can still
  // show the previous run's "completed" message as stale elements, so wait
  // for the run to start (Stop Workflow appears) and to end (it goes away)
  // before reading the result.
  const stop = t.page.getByRole("button", { name: "Stop Workflow" }).first();
  await stop.waitFor({ timeout: 60000 }).catch(() => t.log("no running state seen within 60 s"));
  await stop.waitFor({ state: "detached", timeout: RUN_TIMEOUT });
  await t.waitForText("Workflow completed successfully.", {
    timeout: RUN_TIMEOUT,
    failOn: ["Errors occurred, check log file."],
  });
  await t.shot(`${workflowPage}-${identify ? "identify" : "deconvolve"}-finished`);
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
  await runWorkflow(t, "TopDown", [["MS data", "Deconvolution-only_ProteinMix_LC-MS.mzML"]], false);
  await checkViewer(t, "TopDownViewer", "Deconvolution-only_ProteinMix_LC-MS");

  await runWorkflow(t, "TopDown", [["MS data", "Native-identified_GapA_HCD.mzML"]], true);
  const log = await t.page.locator("body").innerText();
  if (!log.includes("FLASHAnnotator")) {
    t.log(log.slice(-4000));
    throw new Error("The identification did not reach FLASHAnnotator");
  }
  await checkViewer(t, "TopDownViewer", "Native-identified_GapA_HCD");
};
