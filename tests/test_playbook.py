from mpatlas.features import featurize_seq
from mpatlas.playbook import decide
from mpatlas.wrap import wrap_score


def test_high_purify_score_means_standard_path():
    tm = "MIFLFVLLAIVVLITGLLMLNTSNSPYIRLIAFILLVITGLIMLK" * 2
    feat = featurize_seq(tm)
    wrap = wrap_score(tm, feat=feat)
    adv = decide(0.8, wrap, feat)
    assert adv.action == "standard"


def test_low_score_wrap_amenable_is_rescue():
    tm = "MIFLFVLLAIVVLITGLLMLNTSNSPYIRLIAFILLVITGLIMLK" * 2
    feat = featurize_seq(tm)
    wrap = wrap_score(tm, feat=feat)
    assert wrap.amenable
    adv = decide(0.2, wrap, feat)
    assert adv.action == "wrap_rescue"


def test_long_soluble_low_score_is_redesign():
    seq = (
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIA"
        * 14
    )
    feat = featurize_seq(seq)
    wrap = wrap_score(seq, feat=feat)
    adv = decide(0.2, wrap, feat)
    assert not wrap.amenable
    assert adv.action == "redesign"
    assert adv.levers
