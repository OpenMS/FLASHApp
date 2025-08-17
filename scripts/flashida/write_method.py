import json
import sys

import xml.etree.ElementTree as ET

############################
# default paramter values #
###########################
#
# Mandatory keys for each parameter
# key: a unique identifier
# value: the default value
#
# Optional keys for each parameter
# name: the name of the parameter
# hide: don't show the parameter in the parameter section (e.g. for input/output files)
# options: a list of valid options for the parameter
# min: the minimum value for the parameter (int and float)
# max: the maximum value for the parameter (int and float)
# step_size: the step size for the parameter (int and float)
# help: a description of the parameter
# widget_type: the type of widget to use for the parameter (default: auto)
# advanced: whether or not the parameter is advanced (default: False)

DEFAULTS = [
    {"key": "in", "value": [], "help": "Input files for Python Script.", "hide": True},
    {'key': 'TopN', 'name': 'TopN', 'value': 3},
    {'key': 'Duration', 'name': 'Duration', 'value': 67},
    'MS1',
    {'key': 'Analyzer', 'name': 'Analyzer', 'value': 'Orbitrap'},
    {'key': 'FirstMass', 'name': 'FirstMass', 'value': 400},
    {'key': 'LastMass', 'name': 'LastMass', 'value': 2000},
    {'key': 'OrbitrapResolution', 'name': 'OrbitrapResolution', 'value': 120000},
    {'key': 'AGCTarget', 'name': 'AGCTarget', 'value': 800000},
    {'key': 'MaxIT', 'name': 'MaxIT', 'value': 50},
    {'key': 'Microscans', 'name': 'Microscans', 'value': 1},
    {'key': 'DataType', 'name': 'DataType', 'value': 'Centroid'},
    {'key': 'RFLens', 'name': 'RFLens', 'value': 30},
    {'key': 'SourceCID', 'name': 'SourceCID', 'value': 0},
    'MS2',
    {'key': 'Analyzer', 'name': 'Analyzer', 'value': 'Orbitrap'},
    {'key': 'FirstMass', 'name': 'FirstMass', 'value': 200},
    {'key': 'OrbitrapResolution', 'name': 'OrbitrapResolution', 'value': 60000},
    {'key': 'AGCTarget', 'name': 'AGCTarget', 'value': 500000},
    {'key': 'MaxIT', 'name': 'MaxIT', 'value': 118},
    {'key': 'Microscans', 'name': 'Microscans', 'value': 1},
    {'key': 'DataType', 'name': 'DataType', 'value': 'Centroid'},
    {'key': 'Activation', 'name': 'Activation', 'value': 'HCD'},
    {'key': 'CollisionEnergy', 'name': 'CollisionEnergy', 'value': 29},
    'IDA',
    {'key': 'MaxMs2CountPerMs1', 'name': 'MaxMs2CountPerMs1','value': 4},
    {'key': 'QScoreThreshold', 'name': 'QScoreThreshold', 'value': 0.2},
    {'key': 'TQScoreThreshold', 'name': 'TQScoreThreshold', 'value': 0.99},
    {'key': 'MinCharge', 'name': 'MinCharge', 'value': 4},
    {'key': 'MaxCharge', 'name': 'MaxCharge', 'value': 50},
    {'key': 'MinMass', 'name': 'MinMass', 'value': 500},
    {'key': 'MaxMass', 'name': 'MaxMass', 'value': 50000},
    {'key': 'Tolerances', 'name': 'Tolerances', 'value': [10.0, 10.0]},
    {'key': 'TargetLogs', 'name': 'TargetLogs', 'value': [r'C:\Users\KyowonJeong\Desktop\FLASHIdaTmp\test1.log']},
    {'key': 'RTWindow', 'name': 'RTWindow', 'value': 180},
    {'key': 'TargetMode', 'name': 'TargetMode', 'value': 0},
    {'key': 'UseFAIMS', 'name': 'UseFAIMS', 'value': False},
    {'key': 'UseCVQScore', 'name': 'UseCVQScore', 'value': False},
    {'key': 'CycleTime', 'name': 'CycleTime', 'value': 180},
    {'key': 'CVValues', 'name': 'CVValues', 'value': [-10.0, -30.0, -40.0, -50.0, -60.0]},
]

def get_params():
    if len(sys.argv) > 1:
        with open(sys.argv[1], "r") as f:
            return json.load(f)
    else:
        return {}

if __name__ == "__main__":
    params = get_params()

    input_xml = params.pop('input_xml')

    # Create method.xml
    root = ET.Element("MethodParameters")
    subsections = {}
    for key, value in params.items():
        if ':' in key:
            section, param = key.split(':')
            if section not in subsections:
                subsections[section] = ET.SubElement(root, section)
            section = subsections[section]
            tag = ET.SubElement(section, param)
        else:
            tag = ET.SubElement(root, key)
        # Handle list inputs
        if key in ['IDA:Tolerances', 'IDA:CVValues']:
            for item in value.split('\n'):
                element = ET.SubElement(tag, 'double')
                element.text = item.strip()
        elif key in ['IDA:TargetLogs']:
            for item in value.split('\n'):
                element = ET.SubElement(tag, 'string')
                element.text = item.strip().replace('\\\\', '\\')
        elif isinstance(value, bool):
            tag.text = str(value).lower()
        else:
            tag.text = str(value)
    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")
    tree.write(input_xml, encoding="utf-8", xml_declaration=True)