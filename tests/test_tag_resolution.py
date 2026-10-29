import numpy as np
import pandas as pd

from src.parse.tag_resolution import _split_ints, build_tagspace_to_proteoform_map


def test_split_ints_reads_float_formatted_indices():
    assert _split_ints(0.0) == [0]
    assert _split_ints("3;12") == [3, 12]
    assert _split_ints("2.0;5.0") == [2, 5]
    assert _split_ints(np.nan) == []
    assert _split_ints("") == []


def test_map_when_unmatched_tags_make_the_index_column_float():
    # Native-mode FLASHTnT: de novo tags carry no ProteoformIndex, so the
    # column is float and matched tags read as 0.0.
    tags = pd.DataFrame({
        "TagIndex": [0, 1, 2, 3],
        "ProteoformIndex": [np.nan, 0.0, 0.0, np.nan],
    })
    proteins = pd.DataFrame({"ProteoformIndex": [0], "TagIndices": ["1;2"]})
    assert build_tagspace_to_proteoform_map(tags, proteins) == {0: 0}
