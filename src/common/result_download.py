"""Downloading (online) or exporting (local app) the result files of a run.

Used by the Download page, which lists every run, and by the Viewer's
**Download** button, which offers the run on screen. The archive of a run is
kept in the workspace cache (``download_archive``) and built again when one
of its files is newer.
"""

from io import BytesIO
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

import streamlit as st

from src.common.common import export_files, export_folder_input

DOWNLOAD_PAGE = 'content/TopDown/TopDownDownload.py'

# The stored files a run's download holds; a run with any of them is listed
TARGETS = [
    'out_tsv', 'spec1_tsv', 'spec2_tsv', 'spec3_tsv', 'spec4_tsv', 'quant_tsv',
    'toppic_ms1_msalign', 'toppic_ms1_feature', 'toppic_ms2_msalign',
    'toppic_ms2_feature', 'out_deconv_mzML', 'anno_annotated_mzML',
    'FD_parameters_json', 'FTnT_parameters_json', 'tags_tsv', 'protein_tsv',
    'prsms_tsv', 'annotated_tsv', 'fragments_tsv',
]
# The archive the former Download page stored as a result; never part of a download
ARCHIVE = 'download_archive'
ARCHIVE_NAME = 'results.zip'


def is_local():
    """True in the local app, which exports to a folder instead of downloading."""
    return st.session_state.get('location') == 'local'


def verb():
    return 'Export' if is_local() else 'Download'


def downloadable_runs(file_manager):
    """Ids of the runs with result files, in the order the cache lists them."""
    return file_manager.get_results_list(TARGETS, partial=True)


def result_files(file_manager, dataset_id):
    """Paths of the run's result files, without its archive."""
    return [
        Path(p) for p in file_manager.get_all_files_except(dataset_id, [ARCHIVE]).values()
        if p is not None and Path(p).exists() and Path(p).name != ARCHIVE_NAME
    ]


def format_size(n_bytes):
    for unit in ('B', 'KB', 'MB', 'GB'):
        if n_bytes < 1024 or unit == 'GB':
            return f'{n_bytes:.0f} {unit}' if unit == 'B' else f'{n_bytes:.1f} {unit}'
        n_bytes /= 1024


def total_size(files):
    return sum(f.stat().st_size for f in files)


def archive_path(file_manager, dataset_id):
    """Where the run's ZIP archive is kept, next to its files."""
    return Path(file_manager.cache_path, 'files', dataset_id, ARCHIVE_NAME)


def build_archive(files, target):
    """The ZIP of `files` at `target`, built again when missing or older than a file.

    Streamlit runs a download's data callable in another thread, so this
    touches only the file system, never the FileManager's database.
    """
    target = Path(target)
    newest = max((f.stat().st_mtime for f in files), default=0)
    if not (target.exists() and target.stat().st_mtime >= newest):
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_suffix('.part')
        with ZipFile(partial, 'w', ZIP_DEFLATED) as zf:
            for f in files:
                zf.write(f, arcname=f.name)
        partial.replace(target)
    return target


def zip_folders(folders):
    """One ZIP with a folder per run: `folders` maps folder name to files."""
    buffer = BytesIO()
    with ZipFile(buffer, 'w', ZIP_DEFLATED) as zf:
        for folder, files in folders.items():
            for f in files:
                zf.write(f, arcname=f'{folder}/{f.name}')
    return buffer.getvalue()


def run_folders(file_manager, dataset_ids):
    """Folder name (the display name, made unique) to result files, per run."""
    folders = {}
    for dataset_id in dataset_ids:
        base = folder = safe_name(file_manager.get_display_name(dataset_id))
        n = 2
        while folder in folders:
            folder, n = f'{base} ({n})', n + 1
        folders[folder] = result_files(file_manager, dataset_id)
    return folders


def safe_name(name):
    """A display name usable as a file or folder name."""
    text = str(name).replace(': ', ' - ')
    cleaned = ''.join('_' if c in '/\\:*?"<>|' or ord(c) < 32 else c for c in text).strip(' .')
    return cleaned or 'results'


def run_button(file_manager, dataset_id, key, export_dir=None, label=None, primary=True,
               message_area=None):
    """One click downloads the run's ZIP (online) or copies its files to
    ``export_dir/<name>`` (local app), reporting the export in `message_area`."""
    name = file_manager.get_display_name(dataset_id)
    kind = 'primary' if primary else 'secondary'
    if is_local():
        if st.button(label or 'Export', key=f'export_{key}', icon=':material/drive_file_move:',
                     type=kind, width='stretch'):
            with message_area if message_area is not None else st.container():
                export_files(
                    [(f, f.name) for f in result_files(file_manager, dataset_id)],
                    Path(export_dir, safe_name(name)),
                )
        return
    st.download_button(
        label or 'Download ZIP',
        data=lambda files=result_files(file_manager, dataset_id), target=archive_path(file_manager, dataset_id):
            build_archive(files, target).read_bytes(),
        file_name=f'{safe_name(name)}.zip', mime='application/zip',
        key=f'download_{key}', icon=':material/download:', type=kind,
        on_click='ignore', width='stretch',
    )


def all_runs_button(file_manager, dataset_ids, export_dir=None):
    """Download (or export) every run at once."""
    label = f'{verb()} all {len(dataset_ids)} runs'
    if is_local():
        if st.button(label, key='export_all_runs', icon=':material/drive_file_move:', type='primary'):
            for dataset_id in dataset_ids:
                export_files(
                    [(f, f.name) for f in result_files(file_manager, dataset_id)],
                    Path(export_dir, safe_name(file_manager.get_display_name(dataset_id))),
                )
        return
    st.download_button(
        label, data=lambda folders=run_folders(file_manager, dataset_ids): zip_folders(folders),
        file_name='FLASHApp-results.zip', mime='application/zip', key='download_all_runs',
        icon=':material/download:', type='primary', on_click='ignore',
    )


def file_table(files):
    """The files a download holds, with their sizes."""
    st.dataframe(
        {'File': [f.name for f in files], 'Size': [format_size(f.stat().st_size) for f in files]},
        hide_index=True, width='stretch',
    )


def _run_dialog_body(file_manager, dataset_id):
    name = file_manager.get_display_name(dataset_id)
    files = result_files(file_manager, dataset_id)
    if not files:
        st.info(f'**{name}** has no result files to {verb().lower()}.')
        return
    where = 'copied to a folder on this computer' if is_local() else 'downloaded as one ZIP file'
    st.markdown(f'The {len(files)} result files of **{name}** '
                f'({format_size(total_size(files))}), {where}:')
    file_table(files)
    export_dir = export_folder_input('viewer-export') if is_local() else None
    run_button(file_manager, dataset_id, 'viewer_dialog', export_dir)
    try:
        st.page_link(DOWNLOAD_PAGE, label=f'All runs on the {verb()} page', icon=':material/folder_zip:')
    except KeyError:
        # A page loaded on its own (AppTest.from_file) has no navigation to link into
        pass


@st.dialog('Download results', width='large')
def _download_dialog(file_manager, dataset_id):
    _run_dialog_body(file_manager, dataset_id)


@st.dialog('Export results', width='large')
def _export_dialog(file_manager, dataset_id):
    _run_dialog_body(file_manager, dataset_id)


def _run_dialog(file_manager, dataset_id):
    (_export_dialog if is_local() else _download_dialog)(file_manager, dataset_id)


def viewer_button(file_manager, dataset_id=None):
    """The Viewer's button: the shown run for download, else the Download page."""
    if dataset_id is None:
        if st.button(verb(), icon=':material/download:', key='viewer_download', width='content',
                     type='primary', help=f'{verb()} the result files of your runs'):
            st.switch_page(DOWNLOAD_PAGE)
        return
    if st.button(verb(), icon=':material/download:', key='viewer_download', width='content',
                 type='primary', help=f'{verb()} the result files of the run shown below'):
        _run_dialog(file_manager, dataset_id)
