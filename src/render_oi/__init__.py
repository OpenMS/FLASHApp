"""OpenMS-Insight rendering engine for FLASHApp (migration).

Additive replacement for the monolithic ``src/render`` grid: each visualization
is an individual OpenMS-Insight component composed with native Streamlit layout
and a per-experiment StateManager. See :mod:`src.render_oi.deconv_viewer`.
"""

from .deconv_viewer import build_component, render_experiment
from .quant_viewer import build_quant_components, render_experiment_quant
from .tnt_viewer import build_component_tnt, render_experiment_tnt

__all__ = [
    "build_component",
    "render_experiment",
    "build_component_tnt",
    "render_experiment_tnt",
    "build_quant_components",
    "render_experiment_quant",
]
