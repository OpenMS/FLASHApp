"""The History of a result, on Export results (a section of each card) and
Explore results (the **History** button of the experiment row).

It shows how the result came about (``history.timeline``), the tool settings
that differ from the defaults, the processing report of the start that made
it, and **Start again with these settings**, which opens Describe experiment
with the result's files and settings.
"""

import json
from pathlib import Path

import pandas as pd
import streamlit as st

from src import history
from src.workflow.report_ui import report_buttons

WORKFLOW_DIR = 'flashtnt'
CONFIGURE_PAGE = Path('content', 'TopDown', 'TopDownConfigure.py')


def start_again(workspace, file_manager, dataset_id):
    """Open Describe experiment with the files and settings ``dataset_id`` was made with."""
    from src.Workflow import TopDownWorkflow  # the page set up the workspace already
    from src.render.missing_analysis import input_file, run_settings

    runs = history.compared_runs(file_manager, dataset_id) if history.is_comparison(
        file_manager, dataset_id) else [dataset_id]
    params, _ = run_settings(workspace, file_manager, runs[0] if runs else dataset_id)
    try:
        # The settings of the start that made it, kept with its report
        kept = json.loads(Path(history.report_dir(file_manager, dataset_id), 'params.json').read_text(encoding='utf-8'))
        if isinstance(kept, dict):
            params = kept
    except (OSError, ValueError):
        pass
    params = {k: v for k, v in params.items() if k != '_follow_up'}
    files = [f for f in (input_file(workspace, r) for r in runs) if f is not None]
    if files:
        params['mzML-files'] = list(dict.fromkeys(files))
    TopDownWorkflow().start_from(
        params, {}, name=file_manager.get_display_name(dataset_id), page=CONFIGURE_PAGE,
        note='to start again',
    )


def history_section(file_manager, dataset_id, workspace, key):
    """How ``dataset_id`` came about, its settings and its processing report."""
    events = history.timeline(file_manager, dataset_id)
    st.markdown('**How this result came about**')
    if events:
        st.markdown('\n'.join(f'{i}. {history.event_markdown(e)}' for i, e in enumerate(events, start=1)))
    else:
        st.caption('Nothing was recorded about how this result was made (added with Add tool results).')
    if not history.is_comparison(file_manager, dataset_id):
        comparisons = history.comparisons_of(file_manager, dataset_id)
        if comparisons:
            names = ', '.join(f'**{file_manager.get_display_name(c)}**' for c in comparisons)
            st.caption(f'Compared in {names}.')

    rows = history.changed_settings(file_manager, dataset_id, Path(workspace, WORKFLOW_DIR))
    st.markdown('**Settings that differ from the defaults**')
    if rows:
        st.dataframe(pd.DataFrame(rows, columns=['Tool', 'Parameter', 'Default', 'This result']),
                     hide_index=True, width='stretch')
    else:
        st.caption('Every tool ran with its default settings.')

    with st.container(horizontal=True, vertical_alignment='center'):
        if history.has_report(file_manager, dataset_id):
            report_buttons(history.report_dir(file_manager, dataset_id), key=f'{key}-report',
                           workflow=WORKFLOW_DIR)
        if st.button('Start again with these settings', icon=':material/replay:', type='primary',
                     key=f'{key}-again',
                     help='Open Describe experiment with these files and settings. '
                          'The new result is kept next to this one.'):
            start_again(workspace, file_manager, dataset_id)
    if not history.has_report(file_manager, dataset_id):
        st.caption('No processing report was kept for this result: it was made before reports were kept '
                   'with each result.')


@st.dialog('History', width='large')
def history_dialog(file_manager, dataset_id, workspace):
    st.markdown(f'#### {file_manager.get_display_name(dataset_id)}')
    history_section(file_manager, dataset_id, workspace, key=f'history-dialog-{dataset_id}')
    st.caption('The history of every result is also on **Export results**, on its card.')


def history_button(file_manager, dataset_id, workspace):
    """**History** of the result shown, for the experiment row of Explore results."""
    if st.button('History', icon=':material/history:', key='viewer_history',
                 help='How this result came about, its settings and its processing report'):
        history_dialog(file_manager, dataset_id, workspace)
