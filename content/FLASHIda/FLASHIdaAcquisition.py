from pathlib import Path

import streamlit as st

from src import flashida
from src.common.common import local_workspaces_dir, page_setup, save_params
from src.workflow._log_status import classify_log_outcome


params = page_setup()
if not flashida.available():
    st.info('FLASHIda is a developer feature of the desktop app.')
    st.stop()

from src.Workflow import IdaWorkflow  # noqa: E402
from src.render import ida_runs  # noqa: E402

wf = IdaWorkflow()
running = wf.get_workflow_status()['running']
wf_params = wf.parameter_manager.get_parameters_from_json()

st.title('FLASHIda acquisition')
st.caption('🛠️ Developer feature · runs on the instrument PC')

configured = all(str(wf_params.get(k) or '').strip() for k in ('executable', 'raw-files'))
with st.expander('Instrument PC', expanded=not configured and not running):
    if running:
        st.caption('Stop the acquisition to change these.')
        st.text_input('FLASHIda executable', str(wf_params.get('executable') or ''), disabled=True)
        st.text_input('Raw file folder', str(wf_params.get('raw-files') or ''), disabled=True)
    else:
        wf.configure()
        wf.parameter_manager.save_parameters()
        wf_params = wf.parameter_manager.get_parameters_from_json()

raw_dir = str(wf_params.get('raw-files') or '').strip()

st.subheader('Runs')
st.caption('One row per acquisition: the name is the raw file name the sequence writes '
           '(without .raw). FLASHIda is started only for raw files named here, each with '
           'its own method.')


def start_bar() -> None:
    """Run analysis, or Stop while running. Drawn in the fragment of the runs
    table, so its checks follow every edit of the table without a page reload."""
    now_running = wf.get_workflow_status()['running']
    if now_running != running:
        # The acquisition started or ended since the page was drawn
        st.rerun()
    st.write('')
    with st.container(border=True):
        a, b = st.columns([5, 1.4], vertical_alignment='center')
        if running:
            status = flashida.read_status(wf.workflow_dir)
            started = status.get('started')
            a.markdown(
                (f"**Running since {ida_runs._when(started)[:5]}.** " if started else '**Running.** ')
                + 'The table shows what FLASHIda does for each run; the **Viewer** follows it live.'
            )
            if b.button('⏹️ Stop', use_container_width=True):
                wf.stop_workflow()
                st.rerun()
        else:
            current = wf.parameter_manager.get_parameters_from_json()
            problems = ida_runs.start_problems(wf.workflow_dir, current)
            a.markdown(
                '\n'.join(f'- {p}' for p in problems) if problems else
                f"Watches `{str(current.get('raw-files') or '').strip()}` for the newest sequence "
                'and the raw files it writes, and starts FLASHIda for each run of the table.'
            )
            if b.button('🚀 Run analysis', type='primary', use_container_width=True,
                        disabled=bool(problems)):
                wf.start_workflow()
                st.rerun()


ida_runs.runs_table(wf.workflow_dir, local_workspaces_dir(st.session_state.settings), raw_dir, running,
                    footer=start_bar)

log = Path(wf.workflow_dir, 'logs', 'minimal.log')
if log.exists():
    content = log.read_text(encoding='utf-8', errors='replace')
    lines = content.splitlines()
    # A start that ended at once (Stop was not pressed) says why, next to the button
    failed = not running and classify_log_outcome(content) == 'error'
    if failed:
        reason = next((l.strip() for l in reversed(lines)
                       if l.strip() and 'Workflow execution reported failure' not in l), '')
        reason = reason.removeprefix('ERROR:').strip()
        st.error(f'The last acquisition stopped: {reason}' if reason else
                 'The last acquisition stopped with an error, see the log.')
    with st.expander('Log', expanded=failed):
        st.code('\n'.join(lines[-200:]), language=None)

save_params(params)
