import re
import json
import time
import multiprocessing

import streamlit as st
import pyopenms as oms

from time import sleep
from pathlib import Path
from os import makedirs, listdir
from shutil import copyfile, rmtree
from os.path import join, splitext, basename, exists, dirname, abspath

from src.parse.tnt import parseTnT
from src.parse.deconv import parseDeconv
from src.parse.ida import parseIda
from src.workflow.WorkflowManager import WorkflowManager

DEFAULT_THREADS = 8

class TagWorkflow(WorkflowManager):

    def __init__(self) -> None:
        # Initialize the parent class with the workflow name.
        super().__init__("FLASHTnT", st.session_state["workspace"])
        self.tool_name = 'FLASHTaggerViewer'


    def upload(self)-> None:
        t = st.tabs(["MS data", "Database"])
        with t[0]:
            example_data = [f'example-data/flashtagger/example_spectrum_{n}.mzML' for n in ['aqpz', 'antibody']]
            self.ui.upload_widget(key="mzML-files", name="MS data", file_types="mzML", fallback=example_data)
        with t[1]:
            self.ui.upload_widget(key="fasta-file", name="Database", file_types="fasta",
                                  fallback='example-data/flashtagger/example_database.fasta')


    @st.fragment
    def configure(self) -> None:
        # Input File Selection
        self.ui.select_input_file("mzML-files", multiple=True)
        self.ui.select_input_file("fasta-file", multiple=False)

        # Number of threads cannot be selected in online mode
        if st.session_state.location != "online":
            self.ui.input_widget(
                'threads', name='threads', default=multiprocessing.cpu_count(),
                help='The number of threads that should be used to run the tools.'
            )

        # Decoy database size toggle
        self.ui.input_widget(
            'few_proteins', name='Do you expect <100 Proteins?', widget_type='checkbox', default=True,
            help='If set, the decoy database will be 10 times larger than the target database for better FDR estimation resolution. This increases the runtime significantly.'
        )

        # Create tabs for different analysis steps
        t = st.tabs(
            ["**FLASHDeconv**", "**FLASHTnT**"]
        )
        with t[0]:
            # FLASHDeconv Configuration
            self.ui.input_TOPP(
                'FLASHDeconv',
                exclude_parameters = [
                    'ida_log',
                    'write_detail', 'report_FDR', 'quant_method',
                    'mass_error_ppm', 'min_sample_rate', 'min_trace_length',
                    'max_trace_length', 'min_cos', 'type', 'isotope_correction',
                    'reporter_mz_tol', 'only_fully_quantified'
                ],
                display_subsections=True,
                display_tool_name=False
            )
        with t[1]:
            # FLASHTnT Configuration
            self.ui.input_TOPP('FLASHTnT', display_subsections=True, display_tool_name=False)


    def execution(self) -> None:
        # Get input files
        try:      
            in_mzmls = self.file_manager.get_files(self.params["mzML-files"])
        except ValueError:
            self.logger.log('Please select at least one mzML file.')  
            return
        try: 
            database = self.file_manager.get_files(self.params["fasta-file"])
        except ValueError:
            self.logger.log('Please select a database.')  
            return
        
        # Make sure output directory exists
        base_path = dirname(self.workflow_dir)
        
        # Define output directory
        if 'threads' in self.executor.parameter_manager.get_parameters_from_json():
            threads = self.executor.parameter_manager.get_parameters_from_json()['threads']
        else:
            threads = DEFAULT_THREADS

        # Process files in sequence
        for in_mzml in in_mzmls:
            
            # Generate output folder
            current_base = splitext(basename(in_mzml))[0]
            current_time = time.strftime("%Y%m%d-%H%M%S")
            dataset_id = '%s_%s'%(current_base, current_time)
            folder_path = join(base_path, 'FLASHTaggerOutput', '%s_%s'%(current_base, current_time))
            if exists(folder_path):
                rmtree(folder_path)
            makedirs(folder_path)

            self.logger.log(f'Processing {current_base}:')

            # Define output files
            out_tsv = join(folder_path, f'out.tsv')
            out_deconv = join(folder_path, f'out_deconv.mzML')
            out_anno = join(folder_path, f'anno_annotated.mzML')
            out_spec1 = join(folder_path, f'spec1.tsv')
            out_spec2 = join(folder_path, f'spec2.tsv')
            out_spec3 = join(folder_path, f'spec3.tsv')
            out_spec4 = join(folder_path, f'spec4.tsv')
            out_quant = join(folder_path, f'quant.tsv')
            out_msalign1 = join(folder_path, f'toppic_ms1.msalign')
            out_msalign2 = join(folder_path, f'toppic_ms2.msalign')
            out_feature1 = join(folder_path, f'toppic_ms1.feature')
            out_feature2 = join(folder_path, f'toppic_ms2.feature')

            out_prsm = join(folder_path, f'prsms.tsv')
            out_db = join(folder_path, f'database.fasta')
            out_tag = join(folder_path, f'tags.tsv')
            out_protein = join(folder_path, f'protein.tsv')

            # Check if a decoy database needs to be generated
            tagger_params = self.executor.parameter_manager.get_parameters_from_json()['FLASHTnT']
            if ((tagger_params.get('prsm_fdr', 1) < 1) or (tagger_params.get('pro_fdr', 1) < 1)):
                # If few proteins are present increase decoy size
                if self.executor.parameter_manager.get_parameters_from_json()['few_proteins']:
                    ratio = 10
                else:
                    ratio = 1
                
                self.logger.log(f"-> Creating decoy database with target:decoy ratio 1:{ratio}...")

                # Run decoy database
                self.executor.run_topp(
                    'DecoyDatabase',
                    {
                        'in' : [database[0]],
                        'out' : [out_db],
                    },
                    custom_params = {
                        'method' : 'shuffle',
                        'shuffle_decoy_ratio' : ratio,
                        'enzyme' : 'no cleavage',
                    }
                )
            else:
                # If no decoy database is needed the database file is copied as is
                copyfile(database[0], out_db)
            
            self.logger.log(f"-> Running FLASHDeconv...")

            # Run FLASHDeconv (1/2)
            self.executor.run_topp(
                'FLASHDeconv',
                input_output={
                    'in' : [in_mzml],
                    'out' : [out_tsv],
                    'out_spec1' : [out_spec1],
                    'out_spec2' : [out_spec2],
                    'out_spec3' : [out_spec3],
                    'out_spec4' : [out_spec4],
                    'out_mzml' : [out_deconv],
                    'out_quant' : [out_quant],
                    'out_annotated_mzml' : [out_anno],
                    'out_msalign1' : [out_msalign1],
                    'out_msalign2' : [out_msalign2],
                    'out_feature1' : [out_feature1],
                    'out_feature2' : [out_feature2],
                },
                custom_params = {
                    'threads' : threads
                }
            )

            self.logger.log(f"-> Running FLASHTnT...")

            # Run FLASHTnT (2/2)
            self.executor.run_topp(
                'FLASHTnT',
                input_output={
                    'in' : [out_deconv],
                    'fasta' : [out_db],
                    'out_tag' :  [out_tag],
                    'out_pro' :  [out_protein],
                    'out_prsm' : [out_prsm]
                },
                custom_params = {
                    'threads' : threads
                }
            )

            self.logger.log(f"-> Processing Results...")

            # Store all files
            for file in listdir(folder_path):
                self.file_manager.store_file(
                    dataset_id, str(file).replace('.', '_'), 
                    Path(folder_path, file), file_name=file
                )
            # Store Settings
            FDsettings = self.executor.parameter_manager.get_parameters_from_json()['FLASHDeconv']
            self.file_manager.store_data(
                dataset_id, 'FD_parameters', FDsettings
            )
            json_file = Path(folder_path, 'FD_parameters.json')
            with open(json_file, 'w') as f:
                json.dump(FDsettings, f)
            self.file_manager.store_file(
                dataset_id, 'FD_parameters_json', json_file, 
                file_name='FD_parameters.json'
            )
            FTnTsettings = self.executor.parameter_manager.get_parameters_from_json()['FLASHTnT']
            self.file_manager.store_data(
                dataset_id, 'FTnT_parameters', FTnTsettings
            )
            json_file = Path(folder_path, 'FTnT_parameters.json')
            with open(json_file, 'w') as f:
                json.dump(FTnTsettings, f)
            self.file_manager.store_file(
                dataset_id, 'FTnT_parameters_json', json_file, 
                file_name='FTnT_parameters.json'
            )
            # Fetch results
            results = self.file_manager.get_results(
                dataset_id, [
                    'out_deconv_mzML', 'anno_annotated_mzML',
                    'tags_tsv', 'protein_tsv'
                ]
            )
            out_tsv_ms1 = None
            if self.file_manager.result_exists(dataset_id, 'spec1_tsv'):
                out_tsv_ms1 = self.file_manager.get_results(
                    dataset_id, ['spec1_tsv']
                )['spec1_tsv']
            out_tsv_ms2 = None
            if self.file_manager.result_exists(dataset_id, 'spec2_tsv'):
                out_tsv_ms2 = self.file_manager.get_results(
                    dataset_id, ['spec2_tsv']
                )['spec2_tsv']
            # Parse data
            parseDeconv(
                self.file_manager, dataset_id,
                results['out_deconv_mzML'], results['anno_annotated_mzML'], 
                out_tsv_ms1, out_tsv_ms2, logger=self.logger
            )
            parseTnT(
                self.file_manager, dataset_id,
                results['out_deconv_mzML'], results['anno_annotated_mzML'], 
                results['tags_tsv'], results['protein_tsv'], logger=self.logger
            )

            # Remove temporary folder
            rmtree(folder_path)



class DeconvWorkflow(WorkflowManager):

    def __init__(self) -> None:
        # Initialize the parent class with the workflow name.
        super().__init__("FLASHDeconv", st.session_state["workspace"])
        self.tool_name = 'FLASHDeconvViewer'


    def upload(self)-> None:
        self.ui.upload_widget(key="mzML-files", name="MS data", file_types="mzML",
                              fallback=['example-data/flashdeconv/example_fd.mzML'])


    def configure(self) -> None:
        # Input File Selection
        self.ui.select_input_file("mzML-files", multiple=True)

        # Number of threads cannot be selected in online mode
        if st.session_state.location != "online":
            self.ui.input_widget(
                'threads', name='threads', default=multiprocessing.cpu_count(),
                help='The number of threads that should be used to run the tools.'
            )


        # FLASHDeconv Configuration
        self.ui.input_TOPP(
            'FLASHDeconv', exclude_parameters = ['ida_log'], display_subsections=True,
        )


    def execution(self) -> None:
        # Get input files
        try:
            in_mzmls = self.file_manager.get_files(self.params["mzML-files"])
        except ValueError:
            self.logger.log('Please select at least one mzML file.')  
            return
        
        # Define output directory
        base_path = dirname(self.workflow_dir)

        # Set number of threads
        if 'threads' in self.executor.parameter_manager.get_parameters_from_json():
            threads = self.executor.parameter_manager.get_parameters_from_json()['threads']
        else:
            threads = DEFAULT_THREADS

        # Process files in sequence
        for in_mzml in in_mzmls:

            # Generate temporary output folder
            current_base = splitext(basename(in_mzml))[0]
            current_time = time.strftime("%Y%m%d-%H%M%S")
            dataset_id = '%s_%s'%(current_base, current_time)
            folder_path = join(base_path, 'FLASHDeconvOutput', '%s_%s'%(current_base, current_time))
            if exists(folder_path):
                rmtree(folder_path)
            makedirs(folder_path)

            self.logger.log(f'Processing {current_base}:')

            # Define output files
            out_tsv = join(folder_path, f'out.tsv')
            out_deconv = join(folder_path, f'out_deconv.mzML')
            out_anno = join(folder_path, f'anno_annotated.mzML')
            out_spec1 = join(folder_path, f'spec1.tsv')
            out_spec2 = join(folder_path, f'spec2.tsv')
            out_spec3 = join(folder_path, f'spec3.tsv')
            out_spec4 = join(folder_path, f'spec4.tsv')
            out_quant = join(folder_path, f'quant.tsv')
            out_msalign1 = join(folder_path, f'toppic_ms1.msalign')
            out_msalign2 = join(folder_path, f'toppic_ms2.msalign')
            out_feature1 = join(folder_path, f'toppic_ms1.feature')
            out_feature2 = join(folder_path, f'toppic_ms2.feature')

            self.logger.log(f"-> Running FLASHDeconv...")

            # Run FLASHDeconv
            self.executor.run_topp(
                'FLASHDeconv',
                input_output={
                    'in' : [in_mzml],
                    'out' : [out_tsv],
                    'out_spec1' : [out_spec1],
                    'out_spec2' : [out_spec2],
                    'out_spec3' : [out_spec3],
                    'out_spec4' : [out_spec4],
                    'out_mzml' : [out_deconv],
                    'out_quant' : [out_quant],
                    'out_annotated_mzml' : [out_anno],
                    'out_msalign1' : [out_msalign1],
                    'out_msalign2' : [out_msalign2],
                    'out_feature1' : [out_feature1],
                    'out_feature2' : [out_feature2],
                },
                custom_params = {
                    'threads' : threads
                }
            )

            self.logger.log(f"-> Processing Results...")

            # Store all files
            for file in listdir(folder_path):
                self.file_manager.store_file(
                    dataset_id, str(file).replace('.', '_'), 
                    Path(folder_path, file), file_name=file
                )
            results = self.file_manager.get_results(
                dataset_id, 
                ['out_deconv_mzML', 'anno_annotated_mzML']
            )
            out_tsv_ms1 = None
            if self.file_manager.result_exists(dataset_id, 'spec1_tsv'):
                out_tsv_ms1 = self.file_manager.get_results(
                    dataset_id, ['spec1_tsv']
                )['spec1_tsv']
            out_tsv_ms2 = None
            if self.file_manager.result_exists(dataset_id, 'spec2_tsv'):
                out_tsv_ms2 = self.file_manager.get_results(
                    dataset_id, ['spec2_tsv']
                )['spec2_tsv']
            parseDeconv(
                self.file_manager, dataset_id,
                results['out_deconv_mzML'], results['anno_annotated_mzML'], 
                out_tsv_ms1, out_tsv_ms2, logger=self.logger
            )
            
            FDsettings = self.executor.parameter_manager.get_parameters_from_json()['FLASHDeconv']
            self.file_manager.store_data(
                dataset_id, 'FD_parameters', FDsettings
            )
            json_file = Path(folder_path, 'FD_parameters.json')
            with open(json_file, 'w') as f:
                json.dump(FDsettings, f)
            self.file_manager.store_file(
                dataset_id, 'FD_parameters_json', json_file, 
                file_name='FD_parameters.json'
            )

            # Remove temporary folder
            rmtree(folder_path)



class IdaWorkflow(WorkflowManager):

    def __init__(self) -> None:
        # Initialize the parent class with the workflow name.
        super().__init__("FLASHIda", st.session_state["workspace"])
        self.script_path = join('src', 'FLASHIda', 'run.py')
        self.tool_name = 'FLASHIdaRunner'

    def configure(self) -> None:
        self.ui.input_widget(
            key="executable", name="Path to Flash.exe", default='',
            widget_type="text"
        )
        self.ui.input_widget(
            key="executable-secondary", name="Path to secondary Flash.exe", default='',
            widget_type="text"
        )
        self.ui.input_widget(
            key="raw-files", name="Path to raw files", default='',
            widget_type="text"
        )
        self.ui.input_widget(
            key="method-files", name="Path to method files", default='',
            widget_type="text"
        )

    def execution(self) -> None:
        params = self.parameter_manager.get_parameters_from_json()

        # Validate primary FLASHIda executable input
        flashida_path = Path(params['executable'])
        if flashida_path.suffix.lower() != '.exe':
            self.logger.log(
                f'FLASHIda executable was provided with extension '
                f'\'{flashida_path.suffix}\'. Expected \'.exe\''
            )
            return
        if flashida_path.is_file():
            self.logger.log(f'Found FLASHIda executable!')
        else:
            self.logger.log(f'{flashida_path} is not a file.')
            return
        
        # Validate secondary FLASHIda executable input
        flashida_secondary_path = Path(params['executable-secondary'])
        if flashida_secondary_path.suffix.lower() != '.exe':
            self.logger.log(
                f'Secondary FLASHIda executable was provided with extension '
                f'\'{flashida_secondary_path.suffix}\'. Expected \'.exe\''
            )
            return
        if flashida_secondary_path.is_file():
            self.logger.log(f'Found Secondary FLASHIda executable!')
        else:
            self.logger.log(f'{flashida_secondary_path} is not a file.')
            return
        
        # Validate method file input
        methods_folder_path = Path(params['method-files'])
        if methods_folder_path.is_dir() and (params['method-files'] != ''):
            self.logger.log(f'Found methods folder!')
        else:
            self.logger.log(
                f'Method folder \'{methods_folder_path}\' '
                f'is not a folder.'
            )
            return
        
        # Validate raw file input
        raw_folder_path = Path(params['raw-files'])
        if raw_folder_path.is_dir() and (params['raw-files'] != ''):
            self.logger.log(f'Found raw file folder!')
        else:
            self.logger.log(
                f'Raw folder \'{raw_folder_path}\' '
                f'is not a folder.'
            )
            return
        
        # Find existing raw files
        ign_raws, ign_methods, ign_secondary_flags = self._find_raws(
            raw_folder_path
        )
        if len(ign_raws) > 0:
            self.logger.log(
                'Found the following existing raw files that match the scheme:'
            )
            for i, (file, method, secondary) in enumerate(
                zip(ign_raws, ign_methods, ign_secondary_flags)
            ):
                self.logger.log(
                    f"{i+1}:\t{file}\t({method}.xml\t"
                    f"{'secondary' if secondary else 'primary'})"
                )
            self.logger.log('Ignoring these files!')

        self.logger.log('Listening for new raw files...')
        while(True):
            # Scan every 1s
            sleep(1)
            
            # Search for new raws
            new_raws, new_methods, new_secs = self._find_raws(raw_folder_path)
            for raw, method, secondary in zip(new_raws, new_methods, new_secs):
                if raw not in ign_raws:
                    break
            else:
                continue
            
            self.logger.log(f'Detected new raw \'{raw}\'')

            # Ignore raw in future cycles
            ign_raws.append(raw)
            ign_methods.append(method)
            ign_secondary_flags.append(secondary)
            
            # Validate method
            method_path = Path(methods_folder_path, f'{method}.xml')
            raw_path = Path(raw_folder_path, raw)
            self.logger.log(secondary)
            self.logger.log(flashida_secondary_path)
            exe_path = flashida_secondary_path if secondary else flashida_path
            if method_path.is_file():
                self.logger.log(f'Found method \'{method_path}\'!')
                self.logger.log(f'Starting FLASHDeconv...')
                self.executor.run_command(
                    [exe_path, '-m', method_path, '-r', raw_path],
                    cwd = exe_path.parent
                )
                self.logger.log('Listening for new raw files...')

            else:
                self.logger.log(
                    f'Method \'{method_path}\' is not valid. Ignoring...'
                )
                self.logger.log('Listening for new raw files...')


    def _find_raws(self, raw_path):
        # Find existing raw files
        raws = []
        methods = []
        secondary_flags = []
        method_pattern = r'.*FLASHIda_([^.]+)\.raw'
        for file in listdir(raw_path):
            if not Path(raw_path, file).is_file():
                continue
            match = re.search(method_pattern, file)
            if match:
                raws.append(str(file))
                full_suffix = match.group(1)
                parts = full_suffix.split('_', 1)
                method = parts[0]
                methods.append(method)
                suffix = parts[1] if len(parts) > 1 else None
                if suffix is not None:
                    suffix = suffix.split('_')[0] if '_' in suffix else suffix
                    suffix = suffix.split('.')[0] if '.' in suffix else suffix
                secondary = suffix == '2' if suffix else False
                secondary_flags.append(secondary)
        return raws, methods, secondary_flags



class IdaSimulatorWorkflow(WorkflowManager):

    def __init__(self) -> None:
        # Initialize the parent class with the workflow name.
        super().__init__("FLASHIdaSimulator", st.session_state["workspace"])
        self.script_path = join('scripts', 'flashida', 'write_method.py')
        self.tool_name = 'FLASHIdaRunner'
    
    def upload(self)-> None:
        self.ui.upload_widget(key="mzML-files", name="MS data", file_types="mzML")

    def configure(self) -> None:
        
        self.ui.select_input_file("mzML-files", name='Dataset')

        self.ui.input_widget(
            key="executable", name="Path to Flash.exe", default='',
            widget_type="text"
        )

        self.ui.input_python(self.script_path)

    def execution(self) -> None:
        params = self.parameter_manager.get_parameters_from_json()

        # Make sure output directory exists
        base_path = dirname(self.workflow_dir)
        
        # Get input files
        in_mzml = self.file_manager.get_files(self.params["mzML-files"])[0]

        # Generate output folder
        current_base = splitext(basename(in_mzml))[0]
        current_time = time.strftime("%Y%m%d-%H%M%S")
        dataset_id = '%s_%s'%(current_base, current_time)
        folder_path = join(base_path, 'FLASHIdaOutput', '%s_%s'%(current_base, current_time))
        makedirs(folder_path)

        # Generate temp paths for output files
        input_txt = join(folder_path, 'simulation_data.txt')
        input_xml = join(folder_path, 'method.xml')
        output_tsv = join(folder_path, 'simulation_results.tsv')


        # Convert input mzML to input format
        self.logger.log('Converting mzML to simulation input format...')
        exp = oms.MSExperiment()
        oms.MzMLFile().load(in_mzml, exp)
        output = []
        for s in exp.getSpectra():
            if s.getMSLevel() > 1:
                continue
            output.append(f'Spec\t{s.getRT()}\n')
            output += [f'{mz}\t{intensity}\n' for mz, intensity in zip(*s.get_peaks())]
        with open(input_txt, 'w') as f:
            f.writelines(output)

        # Write method.xml
        self.logger.log('Generating parameter file...')
        self.executor.run_python(
            self.script_path, {'input_xml' : join(folder_path, 'method.xml')}
        )

        # Run simulator
        self.logger.log('Running FLASHIda simulator...')
        self.executor.run_command(
            [params['executable'], abspath(input_txt), abspath(output_tsv), abspath(input_xml)],
            cwd = Path(params['executable']).parent
        )

        # Store all files
        for file in listdir(folder_path):
            self.file_manager.store_file(
                dataset_id, str(file).replace('.', '_'), 
                Path(folder_path, file), file_name=file
            )
    
        # Fetch results
        results = self.file_manager.get_results(
            dataset_id, ['simulation_results']
        )
        
        # Parse data
        parseIda(
            self.file_manager, dataset_id, results['simulation_results']
        )

        # Remove temporary folder
        rmtree(folder_path)


class QuantWorkflow(WorkflowManager):

    def __init__(self) -> None:
        # Initialize the parent class with the workflow name.
        super().__init__("FLASHQuant", st.session_state["workspace"])
        self.tool_name = 'FLASHQuantViewer'