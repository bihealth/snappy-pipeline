"""The PureCN panel converts the numbers that GATK writes as strings in a GenomicsDB workspace."""

import importlib.util
import json
from pathlib import Path

SCRIPT = (
    Path(__file__).parents[3]
    / "snappy_wrappers/wrappers/purecn/create_panel/genomicsdb_integers.py"
)

#: As GATK 4.6.2.0 GenomicsDBImport writes them
VIDMAP = {
    "fields": [
        {"name": "AD", "type": ["Integer"], "length": [{"variable_length_descriptor": "R"}]}
    ],
    "contigs": [{"name": "12", "length": "133851895", "tiledb_column_offset": "0"}],
}
CALLSET = {
    "callsets": [
        {"sample_name": "N1", "row_idx": "0", "idx_in_file": "0", "stream_name": "N1_stream"},
        {"sample_name": "N2", "row_idx": "1", "idx_in_file": "0", "stream_name": "N2_stream"},
    ]
}


def _script():
    spec = importlib.util.spec_from_file_location("genomicsdb_integers", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_genomicsdb_numbers_become_integers(tmp_path):
    (tmp_path / "vidmap.json").write_text(json.dumps(VIDMAP))
    (tmp_path / "callset.json").write_text(json.dumps(CALLSET))

    _script().write_integers(str(tmp_path))

    vidmap = json.loads((tmp_path / "vidmap.json").read_text())
    assert vidmap["contigs"] == [{"name": "12", "length": 133851895, "tiledb_column_offset": 0}]
    assert vidmap["fields"] == VIDMAP["fields"]
    callsets = json.loads((tmp_path / "callset.json").read_text())["callsets"]
    assert [(c["row_idx"], c["idx_in_file"]) for c in callsets] == [(0, 0), (1, 0)]
    assert [c["stream_name"] for c in callsets] == ["N1_stream", "N2_stream"]
