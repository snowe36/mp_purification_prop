from pathlib import Path

from mpatlas.detergents import canonicalize, prior_for
from mpatlas.ingest import gpcrdb, screens


def test_canonicalize_ddm_aliases():
    assert canonicalize("DDM") == "DDM"
    assert canonicalize("n-dodecyl-β-D-maltoside") == "DDM"
    assert canonicalize("n-Dodecyl-b-D-maltopyranoside") == "DDM"
    assert canonicalize("LMNG") == "LMNG"
    assert canonicalize("Triton X-100") == "TX100"


def test_screens_hogbom_and_kotov_counts():
    rows = screens.load()
    counts = screens.screen_counts(rows)
    assert counts["hogbom2017"]["targets"] == 60
    assert counts["hogbom2017"]["detergents"] == 16
    assert counts["kotov2019"]["targets"] == 9
    assert counts["kotov2019"]["detergents"] >= 80
    assert counts["lantez2015"]["findings"] >= 6
    assert any((r.extra or {}).get("outcome") == "unfold" for r in rows if r.source == "kotov2019")


def test_gpcrdb_csv_solvent_type_is_the_chemical(tmp_path: Path):
    path = tmp_path / "gpcrdb.csv"
    path.write_text(
        "Receptor entry_name (uniprot),PDB,solvent type (detergent/additive),"
        "Solvent name,concentration,conc unit\n"
        "5ht1b_human,4IAQ,DDM,detergent,1,%w/v\n"
        "5ht1b_human,4IAQ,CHS,additive,0.2,%w/v\n"
    )
    rows = gpcrdb.load(path)
    assert len(rows) == 2
    names = {(r.extra or {}).get("detergent") for r in rows}
    assert names == {"DDM", "CHS"}
    assert all(r.question == "B-conditions" for r in rows)
    assert all(r.architecture == "gpcr" for r in rows)


def test_seven_tm_gets_chs():
    from mpatlas.features import AA, FeatureRow

    def feat(n_tm: int) -> FeatureRow:
        return FeatureRow(
            length=400,
            kd_mean=1.0,
            net_charge=0.0,
            n_tm=n_tm,
            tm_frac=0.3,
            longest_tm=22,
            longest_loop=40,
            extra_frac=0.4,
            n_sequons=0,
            n_cys=2,
            c_term_extra=1.0,
            n_term_extra=1.0,
            composition=tuple(0.05 for _ in AA),
        )

    assert "CHS" in prior_for(feat(7)).extract
    assert "CHS" not in prior_for(feat(12)).extract
