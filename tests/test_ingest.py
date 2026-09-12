from pathlib import Path

from mpatlas.census import run_census
from mpatlas.ingest import purificationdb, targettrack

FIXTURES = Path(__file__).parent / "fixtures"


def test_targettrack_fixture_keeps_nycomps_and_drops_soluble():
    rows = targettrack.load(FIXTURES / "targettrack_tiny.xml", tm_only=True)
    ids = {r.accession for r in rows}
    assert "NY001" in ids
    assert "NY002" in ids
    assert "SOL1" not in ids
    purified = [r for r in rows if r.status == "purified"]
    assert len(purified) == 1
    assert purified[0].label == 1.0


def test_purificationdb_flags_membrane_detergents():
    rows = purificationdb.load(FIXTURES / "purification_tiny.csv")
    assert len(rows) == 3
    dets = [r.extra.get("mp_detergent") for r in rows]
    assert "ddm" in dets
    assert "lmng" in dets
    assert rows[0].question == "B-conditions"


def test_census_stop_rule_on_tiny_bags(tmp_path: Path):
    tt = targettrack.load(FIXTURES / "targettrack_tiny.xml", tm_only=True)
    pdb = purificationdb.load(FIXTURES / "purification_tiny.csv")
    result = run_census(
        {
            "curnow": [],
            "curnow_unlabelled": [],
            "gfp": [],
            "uniprot": [],
            "mpstruc": [],
            "pdbtm": [],
            "topdb": [],
            "targettrack": tt,
            "purificationdb": pdb,
            "gpcrdb": [],
            "screens": [],
        },
        reports_dir=tmp_path,
        processed_dir=tmp_path,
        write_processed=False,
    )
    assert result["stop"]["B"] == "census_only"
    assert result["stop"]["B-conditions"] == "lookup_only"
    assert result["targettrack"]["purified"] == 1
