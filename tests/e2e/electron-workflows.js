// FLASHApp scenario for electron/scripts/e2e.js: runs FLASHDeconv and FLASHTnT
// (DecoyDatabase -> FLASHDeconv -> FLASHTnT -> FLASHAnnotator) on the bundled
// example data inside the Electron window, then opens each Viewer and checks
// that the Vue component draws its plots and tables. FLASHQuant is view-only,
// so its Viewer is checked on the example workspace.
const RUN_TIMEOUT = 45 * 60 * 1000;

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
  const frame = await t.waitForComponent(".js-plotly-plot, .tabulator");
  // Let every cell of the grid finish drawing before judging it.
  await t.page.waitForTimeout(5000);
  const text = await frame.locator("body").innerText();
  if (text.includes("WebGL is not supported")) throw new Error(`${viewerPage}: WebGL plots did not render`);
  const plots = await frame.locator(".js-plotly-plot").count();
  const tables = await frame.locator(".tabulator").count();
  t.log(`${viewerPage}: ${plots} plots, ${tables} tables`);
  if (plots === 0 || tables === 0) throw new Error(`${viewerPage}: expected plots and tables, got ${plots} plots and ${tables} tables`);
  if (extraSelector && (await frame.locator(extraSelector).count()) === 0)
    throw new Error(`${viewerPage}: ${extraSelector} missing`);
  await t.shot(`${viewerPage}-${experiment || "default"}`);
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
