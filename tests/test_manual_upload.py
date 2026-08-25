"""
Tests for the manual result upload of FLASHDeconv results.

Uploading result files through the "Manual Result Upload" tab used to fail
with "AttributeError: 'UploadedFile' object has no attribute 'suffix'":
FileManager.store_file derived the stored file name from `file.suffix`, which
a Streamlit `UploadedFile` (a BytesIO subclass carrying the original name in
`name`) does not have.

The dataset id the page derives from the file name is covered as well.
FLASHDeconv names its output files `out_deconv.mzML` and
`anno_annotated.mzML`, so taking the part in front of the suffix filed the two
halves of one run under the datasets "out" and "anno", where neither could be
parsed. The mapping is a literal in the page (which cannot be imported without
a Streamlit runtime), so it is reproduced here.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from io import BytesIO
from pathlib import Path

import pytest

from src.workflow.FileManager import FileManager


# content/FLASHDeconv/FLASHDeconvWorkflow.py :: UPLOAD_FILE_TYPES
UPLOAD_FILE_TYPES = (
    ('deconv.mzML', 'out_deconv_mzML', 'out'),
    ('annotated.mzML', 'anno_annotated_mzML', 'anno'),
    ('spec1.tsv', 'spec1_tsv', ''),
    ('spec2.tsv', 'spec2_tsv', ''),
)


class FakeUploadedFile(BytesIO):
    """Stand-in for Streamlit's UploadedFile: a BytesIO with a file name."""

    def __init__(self, name, data=b'data'):
        super().__init__(data)
        self.name = name


@pytest.fixture
def file_manager(tmp_path):
    return FileManager(tmp_path, Path(tmp_path, 'cache'))


def store_upload(file_manager, file_names, default_dataset='uploaded'):
    """The storing loop of the manual result upload tab."""
    for file_name in file_names:
        for suffix, name_tag, tool_prefix in UPLOAD_FILE_TYPES:
            if not file_name.endswith(suffix):
                continue
            experiment = file_name[:-len(suffix)].rstrip('_')
            if experiment in ('', tool_prefix):
                experiment = default_dataset
            file_manager.store_file(
                experiment, name_tag, FakeUploadedFile(file_name)
            )
            break
        else:
            raise AssertionError(f'unmatched file name: {file_name}')


def test_store_uploaded_file_keeps_extension(file_manager):
    file_manager.store_file(
        'sample', 'out_deconv_mzML', FakeUploadedFile('out_deconv.mzML', b'mzML')
    )

    stored = file_manager.get_results('sample', ['out_deconv_mzML'])['out_deconv_mzML']
    assert stored.name == 'out_deconv_mzML.mzML'
    assert stored.read_bytes() == b'mzML'


def test_store_path_input_still_works(file_manager, tmp_path):
    source = Path(tmp_path, 'spec1.tsv')
    source.write_text('a\tb\n')

    file_manager.store_file('sample', 'spec1_tsv', source, remove=False)

    assert file_manager.get_results('sample', ['spec1_tsv'])['spec1_tsv'].name == 'spec1_tsv.tsv'


def test_unchanged_output_names_land_in_one_dataset(file_manager):
    store_upload(file_manager, ['out_deconv.mzML', 'anno_annotated.mzML', 'spec1.tsv'])

    assert file_manager.get_results_list(
        ['out_deconv_mzML', 'anno_annotated_mzML', 'spec1_tsv']
    ) == ['uploaded']


def test_renamed_experiments_stay_separate(file_manager):
    store_upload(
        file_manager,
        ['a_deconv.mzML', 'a_annotated.mzML', 'b_deconv.mzML', 'b_annotated.mzML']
    )

    assert sorted(file_manager.get_results_list(
        ['out_deconv_mzML', 'anno_annotated_mzML']
    )) == ['a', 'b']
