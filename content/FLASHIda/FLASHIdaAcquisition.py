import time
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
from src import flashida_folder  # noqa: E402
from src.render import ida_runs, ida_sections as sections  # noqa: E402

wf = IdaWorkflow()
running = wf.get_workflow_status()['running']
wf_params = wf.parameter_manager.get_parameters_from_json()

st.title('FLASHIda acquisition')
st.caption('🛠️ Developer feature · runs on the instrument PC')

configured = all(str(wf_params.get(k) or '').strip() for k in ('executable', 'raw-files', 'method-folder'))
with st.expander('Instrument PC', expanded=not configured and not running):
    if running:
        st.caption('Stop the acquisition to change these.')
        st.text_input('FLASHIda executable', str(wf_params.get('executable') or ''), disabled=True)
        st.text_input('Raw file folder', str(wf_params.get('raw-files') or ''), disabled=True)
        st.text_input('Method folder', str(wf_params.get('method-folder') or ''), disabled=True)
    else:
        wf.configure()
        wf.parameter_manager.save_parameters()
        wf_params = wf.parameter_manager.get_parameters_from_json()

raw_dir = str(wf_params.get('raw-files') or '').strip()
method_dir = str(wf_params.get('method-folder') or '').strip()

# Sequence: the runs table, FLASHIda for each raw file as it appears. Single
# run: FLASHIda at once for one raw file the instrument is writing already.
SEQUENCE, SINGLE = 'Sequence', 'Single run'
running_single = running and bool(flashida.read_status(wf.workflow_dir).get('single'))
if running:
    mode = SINGLE if running_single else SEQUENCE
    # Kept for the switch, which is hidden while running: the page stays in this mode afterwards
    st.session_state['ida-mode'] = mode
else:
    mode = st.segmented_control(
        'Mode', [SEQUENCE, SINGLE], key='ida-mode', label_visibility='collapsed',
        default=None if 'ida-mode' in st.session_state else SEQUENCE,
        help='**Sequence**: FLASHIda for each run of the table as its raw file appears. '
             '**Single run**: FLASHIda now for one raw file the instrument is writing already.',
    ) or SEQUENCE

if mode == SEQUENCE:
    st.subheader('Runs')
    st.caption('One row per acquisition: the name is the raw file name the sequence writes '
               '(without .raw). FLASHIda is started only for raw files named here, each with '
               'the method picked from the method folder.')


def start_bar() -> None:
    """Run analysis, or Stop while running. Drawn in the fragment of the runs
    table, so its checks follow every edit of the table without a page reload."""
    now_running = wf.get_workflow_status()['running']
    if now_running != running:
        # The acquisition started or ended since the page was drawn
        st.rerun()
    st.write('')
    ida_runs.etd_notice(wf.workflow_dir, method_dir)
    with st.container(border=True):
        a, b = st.columns([5, 1.4], vertical_alignment='center')
        if running:
            status = flashida.read_status(wf.workflow_dir)
            started = status.get('started')
            a.markdown(
                (f"**Running since {ida_runs._when(started)[:5]}.** " if started else '**Running.** ')
                + 'The table shows what FLASHIda does for each run; the **Viewer** follows it live.'
            )
            if st.session_state.get('ida-confirm-stop'):
                # A stray click must not end an acquisition: Stop asks first
                st.warning('**Stop the acquisition?** FLASHIda is stopped now and not started for '
                           'the runs still waiting.', icon='⏹️')
                y, n, _ = st.columns([1.4, 1.4, 4])
                if y.button('⏹️ Yes, stop', type='primary', use_container_width=True):
                    st.session_state.pop('ida-confirm-stop', None)
                    wf.stop_workflow()
                    st.rerun()
                if n.button('Keep running', use_container_width=True):
                    st.session_state.pop('ida-confirm-stop', None)
                    st.rerun(scope='fragment')
            elif b.button('⏹️ Stop', use_container_width=True):
                st.session_state['ida-confirm-stop'] = True
                st.rerun(scope='fragment')
        else:
            current = wf.parameter_manager.get_parameters_from_json()
            problems = ida_runs.start_problems(wf.workflow_dir, current)
            # A method that fails FLASHApp's check only warns: the bundled schema
            # may be older than the FLASHIda on this PC, so the user can override.
            warnings = ida_runs.method_warnings(wf.workflow_dir, str(current.get('method-folder') or ''))
            a.markdown(
                '\n'.join(f'- {p}' for p in problems + warnings) if problems or warnings else
                f"Watches `{str(current.get('raw-files') or '').strip()}` for the newest sequence "
                'and the raw files it writes, and starts FLASHIda for each run of the table.'
            )
            override = bool(warnings) and not problems
            etd = [] if problems else ida_runs.etd_runs(wf.workflow_dir, method_dir)
            if st.session_state.get('ida-confirm-etd') and etd:
                # ETD is switched on in Tune, not by a method: Run asks first
                st.warning('**Is ETD set up in Tune?** '
                           f"{'This method acquires' if len(etd) == 1 else 'These methods acquire'} with ETD, "
                           'which FLASHIda cannot switch on itself. Without it the ETD scans fail.\n\n'
                           + ida_runs.etd_lines(etd), icon='⚡')
                y, n, _ = st.columns([2.2, 1.2, 3.6])
                if y.button('⚡ Yes, ETD is set up: start', type='primary', use_container_width=True):
                    st.session_state.pop('ida-confirm-etd', None)
                    wf.start_workflow()
                    st.rerun()
                if n.button('Cancel', use_container_width=True):
                    st.session_state.pop('ida-confirm-etd', None)
                    st.rerun(scope='fragment')
            elif b.button('⚠️ Run anyway' if override else '🚀 Run analysis',
                          type='secondary' if override else 'primary', use_container_width=True,
                          disabled=bool(problems),
                          help="Starts with the methods as they are, overriding FLASHApp's method check. "
                               'FLASHIda checks each method itself when its run starts.' if override else None):
                if etd:
                    st.session_state['ida-confirm-etd'] = True
                    st.rerun(scope='fragment')
                else:
                    wf.start_workflow()
                    st.rerun()
            else:
                # A question left open while the table changed is not asked again by itself
                st.session_state.pop('ida-confirm-etd', None)

@st.fragment(run_every=2 if running else 5)
def single_run_panel() -> None:
    """One raw file the instrument is writing already, with a method of the
    method folder: FLASHIda starts for it at once. Re-drawn every few seconds,
    so a raw file the instrument starts to write appears without a reload (the
    picks are kept by their keys)."""
    if wf.get_workflow_status()['running'] != running:
        st.rerun()
    st.subheader('Single run')
    if running:
        status = flashida.read_status(wf.workflow_dir)
        single = status.get('single') or {}
        started = status.get('started')
        with st.container(border=True):
            st.markdown(f"**FLASHIda is running for `{single.get('raw', '?')}.raw`** with the method "
                        f"**{single.get('method_file', '?')}**"
                        + (f", since {ida_runs._when(started)[:5]}." if started else '.')
                        + ' The **Viewer** follows it live.')
            if st.session_state.get('ida-confirm-stop'):
                st.warning('**Stop FLASHIda?** It stops steering this acquisition now.', icon='⏹️')
                y, n, _ = st.columns([1.4, 1.4, 4])
                if y.button('⏹️ Yes, stop', type='primary', use_container_width=True):
                    st.session_state.pop('ida-confirm-stop', None)
                    wf.stop_workflow()
                    st.rerun()
                if n.button('Keep running', use_container_width=True):
                    st.session_state.pop('ida-confirm-stop', None)
                    st.rerun(scope='fragment')
            elif st.button('⏹️ Stop'):
                st.session_state['ida-confirm-stop'] = True
                st.rerun(scope='fragment')
        return

    st.caption('For a raw file the instrument is writing already: FLASHIda starts at once with the '
               'method picked here and steers the rest of that acquisition. No sequence is watched. '
               'A raw file the instrument starts to write appears here by itself; for one that is not '
               'written yet, add it as a run in **Sequence**.')
    exe = Path(str(wf_params.get('executable') or '').strip())
    problems = []
    if exe.suffix.lower() != '.exe' or not exe.is_file():
        problems.append('Set the FLASHIda executable under Instrument PC.')
    raws = flashida.raw_files_newest_first(raw_dir) if raw_dir and Path(raw_dir).is_dir() else []
    if not raws:
        problems.append('No raw file in the raw file folder (Instrument PC) yet; it appears here '
                        'once the instrument starts writing it.')
    methods, favs = ida_runs._folder_methods(method_dir)
    if not methods:
        problems.append('Set the method folder under Instrument PC.')
    now = time.time()
    names = [r[0] for r in raws]
    growing = {n for n, t in raws if now - t < 60}
    raw = st.selectbox(
        'Raw file', names, index=0 if names else None, key='ida-single-raw', disabled=not names,
        format_func=lambda n: f'{n}  ·  being written' if n in growing else n,
        help='The raw files of the raw file folder, newest first. One changed in the last minute is '
             'marked: the instrument is writing it.')
    method = ida_runs._method_pick('Method', methods, favs, 'ida-single-method',
                                   ida_runs._default_method(methods, favs, []),
                                   disabled=not methods) if methods else None
    for p in problems:
        st.caption(f'⚠️ {p}')
    try:
        etd = sections.etd_scans(flashida_folder.read(method_dir, method)) if method else []
    except ValueError as e:
        st.error(str(e))
        etd, method = [], None
    if etd:
        st.warning(f"**ETD: set up ETD in Tune before FLASHIda starts.** The method acquires with ETD "
                   f"({', '.join(etd)}), which FLASHIda cannot switch on itself.", icon='⚡')

    def start() -> None:
        flashida.write_single_run(wf.workflow_dir, raw, method)
        wf.start_workflow()
        st.rerun()

    if st.session_state.get('ida-confirm-etd-single') and etd:
        st.warning('**Is ETD set up in Tune?** Without it the ETD scans fail.', icon='⚡')
        y, n, _ = st.columns([2.2, 1.2, 3.6])
        if y.button('⚡ Yes, ETD is set up: start', type='primary', use_container_width=True):
            st.session_state.pop('ida-confirm-etd-single', None)
            start()
        if n.button('Cancel', use_container_width=True, key='ida-single-cancel'):
            st.session_state.pop('ida-confirm-etd-single', None)
            st.rerun(scope='fragment')
    elif st.button('🚀 Start FLASHIda', type='primary', disabled=bool(problems) or not raw or not method):
        if etd:
            st.session_state['ida-confirm-etd-single'] = True
            st.rerun(scope='fragment')
        else:
            start()


if mode == SINGLE:
    single_run_panel()
else:
    ida_runs.runs_table(wf.workflow_dir, local_workspaces_dir(st.session_state.settings), raw_dir, running,
                        footer=start_bar, folder=method_dir)

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
