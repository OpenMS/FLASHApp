# Needed as polars doesnt work well with forking on unix systems
import multiprocessing as mp
if mp.get_start_method(allow_none=True) != "spawn":
    mp.set_start_method("spawn", force=True)

import streamlit as st
from pathlib import Path
import json
# For some reason the windows version only works if this is imported here
import pyopenms

if "settings" not in st.session_state:
        with open("settings.json", "r") as f:
            st.session_state.settings = json.load(f)

if __name__ == '__main__':
    # A local install exports results to a folder instead of downloading them.
    export = "Export" if not st.session_state.settings["online_deployment"] else "Download"
    export_icon = "📤" if export == "Export" else "⬇️"
    pages = {
        "FLASHApp" : [
            st.Page(Path("content", "quickstart.py"), title="Quickstart", icon="👋"),
            #st.Page(Path("content", "user_guide.md").read_text(encoding="utf-8"), unsafe_allow_html=True),
            st.Page(Path("content", "user_guide.py"), title="User Guide", icon="📖"),

        ],
        # One workflow: FLASHDeconv for every run, FLASHTnT when a database is given
        "⚡️ Top-down analysis" : [
            st.Page(Path("content", "TopDown", "TopDownUpload.py"), title="Upload data", icon="📁"),
            st.Page(Path("content", "TopDown", "TopDownConfigure.py"), title="Describe experiment", icon="⚙️"),
            st.Page(Path("content", "TopDown", "TopDownRun.py"), title="Run analysis", icon="🚀"),
            st.Page(Path("content", "TopDown", "TopDownSequenceInput.py"), title="Sequence Input", icon="🧵"),
            st.Page(Path("content", "TopDown", "TopDownViewer.py"), title="Viewer", icon="👀"),
            st.Page(Path("content", "TopDown", "TopDownDownload.py"), title=export, icon=export_icon),
        ],
    }

    # FLASHIda: developer feature of the local app (⚙️ Settings → Developer password)
    from src import flashida
    if flashida.available():
        pages["🧠 FLASHIda"] = [
            st.Page(Path("content", "FLASHIda", "FLASHIdaAcquisition.py"), title="Acquisition", icon="⚙️"),
            st.Page(Path("content", "FLASHIda", "FLASHIdaViewer.py"), title="Viewer", icon="👀"),
        ]

    pg = st.navigation(pages, expanded=True)
    pg.run()
