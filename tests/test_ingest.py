from pathlib import Path

from mpatlas.census import run_census
from mpatlas.ingest import gfp, purificationdb, targettrack, unitmp

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


def test_topdb_reads_numtm_and_seq():
    rows = unitmp.load_topdb(FIXTURES / "topdb_tiny.xml")
    assert len(rows) == 1
    assert rows[0].accession == "TEST_HUMAN"
    assert rows[0].n_tm == 2
    assert rows[0].pdb_id == "1ABC"
    assert rows[0].sequence.startswith("MIFLF")


def test_gfp_fixtures_have_paper_n():
    import pandas as pd
    from mpatlas.paths import FIXTURES

    daley = pd.read_csv(FIXTURES / "daley2005.csv")
    hammon = pd.read_csv(FIXTURES / "hammon2009.csv")
    assert len(daley) == 579
    assert set(daley["label"]) <= {0, 1}
    assert len(hammon) == 313
    assert int(hammon["label"].sum()) == 64
    if "sequence" in hammon.columns:
        assert int((hammon["sequence"].fillna("").str.len() > 10).sum()) == 313
