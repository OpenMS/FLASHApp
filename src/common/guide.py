"""Page furniture that walks users through a workflow.

Each workflow page opens with ``page_intro`` (what the page is for and, folded
away, the background a newcomer needs) and ends with ``next_steps`` (where to
go from here), so a workflow reads as one path: Upload data, Describe
experiment, Run analysis, Viewer.
"""

import streamlit as st

# Page paths per workflow, so links survive a page being renamed in one place only
PAGES = {
    'topdown': {
        'upload': ('content/TopDown/TopDownUpload.py', 'Upload data', '📁'),
        'configure': ('content/TopDown/TopDownConfigure.py', 'Describe experiment', '⚙️'),
        'run': ('content/TopDown/TopDownRun.py', 'Run analysis', '🚀'),
        'sequence': ('content/TopDown/TopDownSequenceInput.py', 'Sequence Input', '🧵'),
        'viewer': ('content/TopDown/TopDownViewer.py', 'Viewer', '👀'),
    },
}


def page_link(tool, key, label=None):
    path, default_label, icon = PAGES[tool][key]
    try:
        st.page_link(path, label=label or default_label, icon=icon)
    except KeyError:
        # A page loaded on its own (AppTest.from_file) has no navigation to link into
        st.markdown(f'{icon} {label or default_label}')


def page_intro(lead, background=None, background_title='Background'):
    """The page's purpose in plain words, with optional background folded away."""
    st.markdown(lead)
    if background:
        with st.expander(f'📖 {background_title}'):
            st.markdown(background)


def next_steps(tool, *steps):
    """Footer linking to the pages that usually come next, each with a reason."""
    st.markdown('---')
    st.markdown('**Where to go next**')
    for key, why in steps:
        cols = st.columns([2, 5])
        with cols[0]:
            page_link(tool, key)
        cols[1].markdown(why)
