"""FLASHIda: intelligent data acquisition for top-down proteomics.

FLASHIda runs on the instrument PC next to a Thermo Orbitrap. The app starts it
for every raw file the acquisition writes (src/Workflow.py:IdaWorkflow) and
runs its simulator on an mzML (IdaSimulatorWorkflow). Both are shown only in
the local app, and only once developer mode is unlocked (``available()``).

Everything but ``available()`` is Streamlit-free, so the workflow process can
use it.
"""

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pyopenms as oms

# The method file FLASHIda reads: section -> [(parameter, default)]. A list
# default is a list parameter, entered one value per line.
METHOD = {
    'MS1': [
        ('TopN', 3), ('Duration', 67), ('Analyzer', 'Orbitrap'),
        ('FirstMass', 400), ('LastMass', 2000), ('OrbitrapResolution', 120000),
        ('AGCTarget', 800000), ('MaxIT', 50), ('Microscans', 1),
        ('DataType', 'Centroid'), ('RFLens', 30), ('SourceCID', 0),
    ],
    'MS2': [
        ('Analyzer', 'Orbitrap'), ('FirstMass', 200), ('OrbitrapResolution', 60000),
        ('AGCTarget', 500000), ('MaxIT', 118), ('Microscans', 1),
        ('DataType', 'Centroid'), ('Activation', 'HCD'), ('CollisionEnergy', 29),
    ],
    'IDA': [
        ('MaxMs2CountPerMs1', 4), ('QScoreThreshold', 0.2), ('TQScoreThreshold', 0.99),
        ('MinCharge', 4), ('MaxCharge', 50), ('MinMass', 500), ('MaxMass', 50000),
        ('Tolerances', [10.0, 10.0]), ('TargetLogs', []), ('RTWindow', 180),
        ('TargetMode', 0), ('UseFAIMS', False), ('UseCVQScore', False),
        ('CycleTime', 180), ('CVValues', [-10.0, -30.0, -40.0, -50.0, -60.0]),
    ],
}
# Element type of the list parameters
LIST_TYPES = {'Tolerances': 'double', 'CVValues': 'double', 'TargetLogs': 'string'}


def method_key(section: str, name: str) -> str:
    """Workflow parameter key of a method parameter."""
    return f'method-{section}-{name}'


def method_widget_default(value):
    """Widget default of a method parameter: lists as one value per line."""
    if isinstance(value, list):
        return '\n'.join(str(v) for v in value)
    return value


def write_method_xml(params: dict, path) -> None:
    """Write FLASHIda's method.xml from the workflow parameters (defaults for
    the ones not set)."""
    root = ET.Element('MethodParameters')
    for section, entries in METHOD.items():
        node = ET.SubElement(root, section)
        for name, default in entries:
            value = params.get(method_key(section, name), method_widget_default(default))
            tag = ET.SubElement(node, name)
            if name in LIST_TYPES:
                for item in str(value).splitlines():
                    if item.strip():
                        ET.SubElement(tag, LIST_TYPES[name]).text = item.strip()
            elif isinstance(value, bool):
                tag.text = str(value).lower()
            else:
                tag.text = str(value)
    tree = ET.ElementTree(root)
    ET.indent(tree, space='  ')
    tree.write(str(path), encoding='utf-8', xml_declaration=True)


# Raw files the acquisition writes: <anything>FLASHIda_<method>[_<n>...].raw;
# n == 2 means the secondary FLASHIda.
_RAW = re.compile(r'.*FLASHIda_([^.]+)\.raw$', re.IGNORECASE)


def parse_raw_name(file_name: str):
    """(method, secondary) of a FLASHIda raw file name, or None."""
    match = _RAW.match(file_name)
    if match is None:
        return None
    method, _, suffix = match.group(1).partition('_')
    secondary = suffix.split('_')[0] == '2' if suffix else False
    return method, secondary


def find_raws(folder) -> dict:
    """{raw file name: (method, secondary)} of the FLASHIda raw files in a folder."""
    raws = {}
    for path in Path(folder).iterdir():
        if path.is_file():
            parsed = parse_raw_name(path.name)
            if parsed is not None:
                raws[path.name] = parsed
    return raws


def write_simulator_input(in_mzml, out_txt) -> int:
    """The MS1 spectra of an mzML in the simulator's text format
    ("Spec<TAB>rt" then "mz<TAB>intensity" per peak). Returns the number of
    spectra written."""
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(in_mzml), exp)
    count = 0
    with open(out_txt, 'w') as f:
        for spectrum in exp.getSpectra():
            if spectrum.getMSLevel() > 1:
                continue
            count += 1
            f.write(f'Spec\t{spectrum.getRT()}\n')
            for mz, intensity in zip(*spectrum.get_peaks()):
                f.write(f'{mz}\t{intensity}\n')
    return count


def available() -> bool:
    """FLASHIda pages are shown in the local app once developer mode is on."""
    import streamlit as st
    from src.common.developer import is_developer

    settings = st.session_state.get('settings') or {}
    return not settings.get('online_deployment', True) and is_developer()
