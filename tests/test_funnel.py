from pathlib import Path

from mpatlas.catalog import Record
from mpatlas.funnel import attrition
from mpatlas.ingest import targettrack


def test_attrition_is_cumulative():
    rows = [
        Record(source="t", question="B", sequence="A", status="cloned"),
        Record(source="t", question="B", sequence="A", status="purified"),
        Record(source="t", question="B", sequence="A", status="in pdb"),
    ]
    counts = attrition(rows)
    assert counts["cloned"] == 3
    assert counts["purified"] == 2
    assert counts["in pdb"] == 1


def test_targettrack_fasta_keeps_membrane_centers(tmp_path: Path):
    fasta = tmp_path / "tt.fasta"
    fasta.write_text(
        ">NYCOMPS-1 $ 1 $ NYCOMPS $ purified $ 1 $\n"
        "MIFLFVLLAIVVLITGLLMLNTSNSPYIRLIAFILLVITGLIMLKIFLFVLLAIVVLITGL\n"
        ">MCSG-1 $ 1 $ MCSG $ purified $ 1 $\n"
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKA\n"
    )
    rows = targettrack.load(fasta, tm_only=True)
    ids = {r.accession for r in rows}
    assert "NYCOMPS-1" in ids
    assert "MCSG-1" not in ids
    assert rows[0].label == 1.0
