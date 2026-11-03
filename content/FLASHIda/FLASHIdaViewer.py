from pathlib import Path

import streamlit as st

from src import flashida
from src.common.common import page_setup
from src.parse import ida_log

page_setup()
if not flashida.available():
    st.info('FLASHIda is a developer feature of the desktop app.')
    st.stop()

from src.Workflow import IdaWorkflow  # noqa: E402
from src.render.ida_live import show_run  # noqa: E402

wf = IdaWorkflow()

st.title('FLASHIda Viewer')
st.caption(
    '🛠️ Developer feature. Follows the logs of a FLASHIda run while it acquires: '
    'the precursors it chose from each survey scan, the scans it commanded for '
    'them and what came of them. Click a precursor to see its scans.'
)

known = flashida.log_folders(wf.workflow_dir)
FOLDER_KEY = 'flashida-viewer-log-folder'
if FOLDER_KEY not in st.session_state:
    st.session_state[FOLDER_KEY] = known[0] if known else ''

c1, c2 = st.columns([3, 1])
with c1:
    st.text_input(
        'Log folder', key=FOLDER_KEY,
        help="FLASHIda's runtime.log_dir (the folder holding the run folders), or one run folder. "
             'The Acquisition page fills it in when it starts FLASHIda.',
    )
    if len(known) > 1:
        picked = st.selectbox(
            'Log folders FLASHIda was started with', ['-'] + known,
            key='flashida-viewer-known-folder',
        )
        if picked != '-' and picked != st.session_state[FOLDER_KEY]:
            st.session_state[FOLDER_KEY] = picked
            st.rerun()
with c2:
    live = st.toggle('Live', value=True, key='flashida-viewer-live',
                     help='Re-read the logs every few seconds.')
    interval = st.select_slider('Every', options=[2, 5, 10, 30], value=5,
                                format_func=lambda s: f'{s} s', key='flashida-viewer-interval',
                                disabled=not live)

folder = Path(st.session_state[FOLDER_KEY].strip()) if st.session_state[FOLDER_KEY].strip() else None
if folder is None:
    st.info('Enter the folder FLASHIda writes its logs to.')
    st.stop()
if not folder.is_dir():
    st.warning(f"'{folder}' is not a folder.")
    st.stop()

FOLLOW = '⏩ Newest run'
cache_root = Path(wf.workflow_dir, 'cache', 'insight-live')


@st.fragment(run_every=interval if live else None)
def runs_view():
    runs = ida_log.find_runs(folder)
    if not runs:
        st.info('No FLASHIda runs in this folder yet. A run appears as soon as FLASHIda starts.')
        return
    names = [FOLLOW] + [r.name for r in runs]
    if st.session_state.get('flashida-viewer-run') not in names:
        st.session_state['flashida-viewer-run'] = FOLLOW
    choice = st.selectbox('Run', names, key='flashida-viewer-run',
                          help='Newest run follows the run FLASHIda is writing.')
    run = runs[0] if choice == FOLLOW else runs[names.index(choice) - 1]
    if choice == FOLLOW:
        st.caption(f'Showing **{run.name}**')
    show_run(run, cache_root, 'viewer')


runs_view()
