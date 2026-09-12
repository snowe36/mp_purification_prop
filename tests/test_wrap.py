from mpatlas.features import SEQUON, featurize_seq
from mpatlas.wrap import wrap_score


def test_sequons_skip_proline():
    assert SEQUON.search("NGT") 
    assert SEQUON.search("NAS")
    assert not SEQUON.search("NPT")


def test_wrap_amenable_needs_tm_belt():
    tm = "MIFLFVLLAIVVLITGLLMLNTSNSPYIRLIAFILLVITGLIMLK" * 2
    soluble = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIA"
    assert wrap_score(tm).amenable
    assert not wrap_score(soluble).amenable
    feat = featurize_seq(tm)
    assert feat.n_tm >= 1
