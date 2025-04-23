import json

import streamlit as st

from src.common.common import page_setup, v_space, save_params
from src.workflow.FileManager import FileManager
from pathlib import Path

COMPONENT_OPTIONS=[
    'MS1 raw heatmap',
    'MS1 deconvolved heatmap',
    'Scan table',
    'Deconvolved spectrum (Scan table needed)',
    'Raw spectrum (Scan table needed)',
    'Mass table (Scan table needed)',
    '3D S/N plot (Mass table needed)',
    'Score Distribution Plot'
    # "Sequence view" and "Internal fragment map" is added when "input_sequence" is submitted
]

COMPONENT_NAMES=[
    'ms1_raw_heatmap',
    'ms1_deconv_heat_map',
    'scan_table',
    'deconv_spectrum',
    'anno_spectrum',
    'mass_table',
    '3D_SN_plot',
    'fdr_plot',
    # "sequence view" and "internal fragment map" added when "input_sequence" is submitted
]

# Setup cache access
file_manager = FileManager(
    st.session_state["workspace"],
    Path(st.session_state['workspace'], 'flashdeconv', 'cache')
)

def set_layout(layout, side_by_side=False):
    file_manager.store_data('layout', 'layout', 
        {
            'layout': layout,
            'side_by_side': side_by_side
        }
    )

def get_layout():
    # Check if layout has been set
    if not file_manager.result_exists('layout', 'layout'):
        return None
    # fetch layout from cache
    layout = file_manager.get_results('layout', 'layout')['layout']

    return layout['layout'], layout['side_by_side'] 

def resetSettingsToDefault(num_of_exp=1):
    st.session_state["layout_setting"] = [[['']]] # 1D: experiment, 2D: row, 3D: column, element=component name
    st.session_state["num_of_experiment_to_show"] = num_of_exp
    for index in range(1, num_of_exp):
        st.session_state.layout_setting.append([['']])
    if file_manager.result_exists('layout', 'layout'):
        file_manager.remove_results('layout')
    st.session_state["edit_mode"] = True


def containerForNewComponent(exp_index, row_index, col_index):

    def isThisComponentUnique(new_component_option):
        if any(col for row in st.session_state.layout_setting[exp_index] for col in row if col==new_component_option):
            st.session_state["component_error_message"] = 'Duplicated component!'
            return False
        else:
            return True

    def addNewComponent():
        new_component_option = 'SelectNewComponent%d%d%d'%(exp_index, row_index, col_index)
        if isThisComponentUnique(st.session_state[new_component_option]):
            st.session_state.layout_setting[exp_index][row_index][col_index] = st.session_state[new_component_option]

    # new component
    st.selectbox("New component to add", ['Select...'] + COMPONENT_OPTIONS,
                 key='SelectNewComponent%d%d%d'%(exp_index, row_index, col_index),
                 on_change=addNewComponent,
                 placeholder='Select...',
                 )


def layoutEditorPerExperiment(exp_index):
    layout_info = st.session_state.layout_setting[exp_index]

    for row_index, row in enumerate(layout_info):
        st_cols = st.columns(len(row)+1 if  len(row)<3 else len(row))
        for col_index, col in enumerate(row):
            if not col: # if empty, add newComponent container
                with st_cols[col_index].container():
                    containerForNewComponent(exp_index, row_index, col_index)
            else:
                with st_cols[col_index]:
                    c1, c2 = st.columns([5, 1])
                    c1.info(col)
                    if c2.button("x", key='DelButton%d%d%d'%(exp_index, row_index, col_index), type='primary'):
                        layout_info[row_index].pop(col_index)
                        st.rerun()

        # new column button
        if len(row) < 3: # limit for #column is 3
            if st_cols[-1].button("***+***", key='NewColumnButton%d%d'%(exp_index, row_index)):
                layout_info[row_index].append('')
                st.rerun()

    # new row button
    if st.button("***+***", key='NewRowButton%d'%exp_index):
        layout_info.append([''])
        st.rerun()


def validateSubmittedLayout(input_layout=None):
    layout_setting = input_layout if input_layout is not None else st.session_state.layout_setting

    # check if submitted layout is empty
    if not any(col for exp in layout_setting for row in exp for col in row if col):
        return 'Empty input'

    # check if submitted layout contains "needed" components
    for exp in layout_setting:
        submitted_components = [col for row in exp for col in row if col]
        required_components = [comp.split('(')[1].split('needed')[0].rstrip() for comp in submitted_components if 'needed' in comp]
        if required_components:
            for required in required_components:
                required_exist = False
                for submitted in submitted_components:
                    if submitted.startswith(required):
                        required_exist = True
                if not required_exist:
                    return 'Required component is missing'
    return ''


def getTrimmedLayoutSetting():
    trimmed_layout_setting = []
    for exp in st.session_state.layout_setting:
        rows = []
        for row in exp:
            cols = []
            for col in row:
                if col:
                    cols.append(COMPONENT_NAMES[COMPONENT_OPTIONS.index(col)])
            if cols:
                rows.append(cols)
        if rows:
            trimmed_layout_setting.append(rows)
    return trimmed_layout_setting


def handleEditAndSaveButtons():
    # if "Edit" button was clicked,
    if "edit_btn_clicked" in st.session_state and st.session_state["edit_btn_clicked"]:
        st.session_state["edit_mode"] = True
        # reset variables based on saved layout setting
        st.session_state["num_of_experiment_to_show"] = len(get_layout()[0]) if get_layout() is not None else 1
        st.session_state["layout_setting"] = [[[COMPONENT_OPTIONS[COMPONENT_NAMES.index(col)]
                                                for col in row if col]
                                               for row in exp if row]
                                              for exp in get_layout()[0]]

    # if "Save" button was clicked,
    if "layout_saved" in st.session_state and st.session_state["layout_saved"]:
        got_error = validateSubmittedLayout()
        st.session_state['save_btn_error_message'] = got_error # to show error msg at the end
        if not got_error:
            # get only submitted info from "layout_setting"
            set_layout(getTrimmedLayoutSetting(), side_by_side=st.session_state['side_by_side_view'])
            st.session_state["edit_mode"] = False


def handleSettingButtons():
    if "reset_btn_clicked" in st.session_state and st.session_state.reset_btn_clicked:
        resetSettingsToDefault()

    if "uploaded_json_file" in st.session_state and st.session_state.uploaded_json_file is not None:
        uploaded_layout = json.load(st.session_state.uploaded_json_file)
        validated = validateSubmittedLayout(uploaded_layout)
        if validated!='':
            st.session_state["component_error_message"] = validated
        else:
            st.session_state.layout_setting = [[[COMPONENT_OPTIONS[COMPONENT_NAMES.index(col)]
                                                 for col in row if col]
                                                for row in exp if row]
                                               for exp in uploaded_layout]
            st.session_state.num_of_experiment_to_show = len(uploaded_layout)


def setSequenceView():
    if 'input_sequence' in st.session_state and st.session_state.input_sequence:
        global COMPONENT_OPTIONS
        COMPONENT_OPTIONS = COMPONENT_OPTIONS + ['Sequence view (Mass table needed)',
                                                 'Internal fragment map (Mass table needed)']
        global COMPONENT_NAMES
        COMPONENT_NAMES = COMPONENT_NAMES + ['sequence_view', 'internal_fragment_map']


# page initialization
params = page_setup()

# when sequence is submitted, add "Sequence View" as a component option
setSequenceView()

# handles "onclick" of buttons
if st.session_state.get("edit_mode") is None:
    st.session_state["edit_mode"] = True
handleSettingButtons()
handleEditAndSaveButtons()

# initialize setting information
if "layout_setting" not in st.session_state:
    resetSettingsToDefault()
# the "num_of_experiment_to_show" changed
elif "num_of_experiment_to_show" in st.session_state and \
        len(st.session_state.layout_setting) != st.session_state.num_of_experiment_to_show:
    resetSettingsToDefault(st.session_state.num_of_experiment_to_show)

### title and setting buttons
c1, c2, c3, c4, c5 = st.columns([6, 1, 1, 1, 1])
c1.title("Layout Manager")

# side-by-side view option for 2 experiments
if 'side_by_side_view' not in st.session_state:
    st.session_state['side_by_side_view'] = False
if (
    ('num_of_experiment_to_show' in st.session_state
     and st.session_state.num_of_experiment_to_show == 2)
    or 
    (not st.session_state.edit_mode
     and (get_layout() is not None and len(get_layout()[0]) == 2))
):
    v_space(1, c2)
    st.session_state['side_by_side_view'] = c2.checkbox(
        "Side-by-Side View", value=st.session_state['side_by_side_view'],
        help="If checked, experiments will be shown side-by-side",
        disabled=(not st.session_state.edit_mode)
    )

# Load existing layout setting file
v_space(1, c3)
c3.button("Load Setting", key="load_btn_clicked")

# Save current layout setting (only after "Saved" button)
v_space(1, c4)
c4.download_button(
    label="Save Setting",
    data=json.dumps(getTrimmedLayoutSetting()),
    file_name='FLASHViewer_layout_settings.json',
    mime='json',
    disabled=(validateSubmittedLayout()!=''),
)

# Reset settings to default
v_space(1, c5)
c5.button("Reset Setting", key="reset_btn_clicked")

### space for File Uploader, when "Load Setting" button is clicked
if "load_btn_clicked" in st.session_state and st.session_state.load_btn_clicked:
    st.file_uploader("Choose a json file", type="json", key="uploaded_json_file")

### Main part
if (not st.session_state.edit_mode) and (get_layout() is not None):
    # show saved-mode
    for index_of_experiment in range(len(get_layout()[0])):
        layout_info_per_experiment = get_layout()[0][index_of_experiment]
        with st.expander("Experiment #%d"%(index_of_experiment+1), expanded=True):
            for row_index, row in enumerate(layout_info_per_experiment):
                st_cols = st.columns(len(row))
                for col_index, col in enumerate(row):
                    st_cols[col_index].info(COMPONENT_OPTIONS[COMPONENT_NAMES.index(col)])
else:
    # show edit-mode
    st.selectbox("**#Experiments to view at once**", [1, 2, 3, 4, 5],
                 key="num_of_experiment_to_show",
    )

    for index_of_experiment in range(st.session_state.num_of_experiment_to_show):
        with st.expander("Experiment #%d"%(index_of_experiment+1)):
            layoutEditorPerExperiment(index_of_experiment)

### buttons for edit/save
_, edit_btn_col, save_btn_col = st.columns([9, 1, 1])
edit_btn_col.button("Edit", key="edit_btn_clicked", disabled=st.session_state.edit_mode)
save_btn_col.button("Save", key="layout_saved", disabled=(not st.session_state.edit_mode))

### showing error/success message
if "save_btn_error_message" in st.session_state and st.session_state.layout_saved:
    error_message = st.session_state["save_btn_error_message"]
    if error_message:
        st.error('Error: '+error_message, icon="🚨")
    else:
        st.success('Layouts Saved', icon="✔️")
if "component_error_message" in st.session_state and st.session_state.component_error_message:
    st.error('Error: ' + st.session_state.component_error_message, icon="🚨")
    del st.session_state["component_error_message"]

### TIPs (TODO: Add image)
st.info("""
**💡 Tips**

- If nothing is set, the default layout will be used in the **👀 Viewer** page

- Don't forget to click "save" on the bottom-right corner to save your setting
""")

save_params(params)
