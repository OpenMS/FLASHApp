import os

import streamlit as st
import streamlit.components.v1 as st_components


# Create a _RELEASE constant. We'll set this to False while we're developing
# the component, and True when we're ready to package and distribute it.
_RELEASE = False


_component_func = None
def get_component_function():
    global _component_func, _RELEASE

    if '_component_func' not in st.session_state:

        if not _RELEASE:
            st.session_state['_component_func'] = st_components.declare_component(
                "flash_viewer_grid",
                url="http://localhost:5173",
        )
        else:
            parent_dir = os.path.dirname(os.path.abspath(__file__))
            build_dir = os.path.join(parent_dir, '..', '..', "js-component", "dist")
            st.session_state['_component_func'] = st_components.declare_component("flash_viewer_grid", path=build_dir)
            
    return st.session_state['_component_func']


class FlashViewerComponent:
    componentArgs = None

    def __init__(self, component_args):
        self.componentArgs = component_args


class PlotlyHeatmap:
    title = None
    showLegend = None

    def __init__(self, title, show_legend=False):
        self.title = title
        self.show_legend = show_legend
        self.componentName = "PlotlyHeatmap"


class Tabulator:
    def __init__(self, table_type):
        if table_type == 'ScanTable':
            self.title = 'Scan Table'
            self.componentName = "TabulatorScanTable"
        elif table_type == 'MassTable':
            self.title = 'Mass Table'
            self.componentName = "TabulatorMassTable"
        elif table_type == 'ProteinTable':
            self.title = 'Protein Table'
            self.componentName = "TabulatorProteinTable"
        elif table_type == 'TagTable':
            self.title = 'Tag Table'
            self.componentName = "TabulatorTagTable"


class PlotlyLineplot:
    def __init__(self, title):
        self.title = title
        self.componentName = "PlotlyLineplot"

class FDRPlotly:
    def __init__(self, title):
        self.title = title
        self.componentName = "FDRPlotly"

class PlotlyLineplotTagger:
    def __init__(self, title):
        self.title = title
        self.componentName = "PlotlyLineplotTagger"


class Plotly3Dplot:
    def __init__(self, title):
        self.title = title
        self.componentName = "Plotly3Dplot"


class SequenceView:
    def __init__(self, title):
        self.title = title
        self.componentName = 'SequenceView'


class InternalFragmentMap:
    def __init__(self, title):
        self.title = title
        self.componentName = 'InternalFragmentMap'


class FLASHQuant:
    def __init__(self):
        self.title = 'QuantVis'
        self.componentName = 'FLASHQuantView'
