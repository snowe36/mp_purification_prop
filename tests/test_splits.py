import numpy as np

from mpatlas.splits import cluster_split, hamming_clusters, random_split, split_ok


def test_random_split_is_disjoint_and_stratified():
    y = np.array([0] * 40 + [1] * 40)
    tr, te = random_split(len(y), y, test_size=0.3, seed=1)
    assert set(tr).isdisjoint(set(te))
    assert split_ok(y, tr, te, min_n=5)
    assert abs(y[tr].mean() - y[te].mean()) < 0.2


def test_hamming_clusters_group_near_sequences():
    a = "A" * 40
    b = "A" * 38 + "CC"
    c = "G" * 40
    labels = hamming_clusters([a, b, c], max_frac=0.12)
    assert labels[0] == labels[1]
    assert labels[2] != labels[0]


def test_cluster_split_holds_out_whole_clusters():
    clusters = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    tr, te = cluster_split(clusters, test_size=0.34, seed=0)
    te_c = set(clusters[te])
    tr_c = set(clusters[tr])
    assert te_c.isdisjoint(tr_c)
