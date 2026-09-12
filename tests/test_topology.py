from mpatlas.topology import predict_topology


def test_single_tm_n_in_is_bitopic():
    # Hydrophobic helix after a KR-rich N-terminus → N-in, C-out type II-ish.
    seq = "MKKKKK" + "LIVLFVLLAIVVLITGLLMLLIVLFVLLAI" + "NNNNNNNDDDDEEEE"
    topo = predict_topology(seq)
    assert topo.n_tm >= 1
    assert topo.architecture in {"bitopic_helical", "polytopic_helical"}
    assert topo.segments[0].kind in {"cyto", "extra", "unknown"}


def test_soluble_has_zero_tm():
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIA"
    topo = predict_topology(seq)
    assert topo.n_tm == 0
