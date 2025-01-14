import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib_logo import logo
import seaborn as sns

import pysam


RC = str.maketrans('ACGT', 'TGCA')


def rev_comp(seq):
    return seq.translate(RC)[::-1]


IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U',
    'S': 'GC', 'W': 'AU', 'R': 'AG', 'Y': 'CU',
    'B': 'CGU', 'H': 'ACU',
    'N': 'ACGU',
}

IUPAC_INV = {
    'A': 'B', 'C': 'D', 'G': 'H', 'U': 'V',
    'S': 'W', 'W': 'S', 'R': 'Y', 'Y': 'R',
    'N': 'N'
}


SE_CARTOON = mpimg.imread(os.path.join(os.path.split(__file__)[0], 'data/se_cartoon.png'))
A3_CARTOON = mpimg.imread(os.path.join(os.path.split(__file__)[0], 'data/a3_cartoon.png'))


def iupac_classify(seq, consensus):
    clss = []
    for o, e in zip(seq, consensus):
        if o in IUPAC[e]:
            clss.append(e)
        
        else:
            clss.append(IUPAC_INV[e])
    return ''.join(clss)


def u5_classify(donor_seq):
    assert not len(donor_seq) % 2
    ws = len(donor_seq) // 2
    return iupac_classify(donor_seq[ws - 2: ws], 'AG')


def u6_classify(donor_seq):
    assert not len(donor_seq) % 2
    ws = len(donor_seq) // 2
    return iupac_classify(donor_seq[ws + 2: ws + 5], 'RAG')


def acceptor_classify(acceptor_seq):
    assert not len(acceptor_seq) % 2
    ws = len(acceptor_seq) // 2
    return iupac_classify(acceptor_seq[ws - 5: ws - 2], 'UGC')


def edit_distance(seq1, seq2):
    ed = 0
    for i, j in zip(seq1, seq2):
        if i != j:
            ed += 1
    return ed


def _seq_test(seq, exp):
    exp = IUPAC[exp]
    return seq in exp


def perc_seq_pos(seqs, pos, nt):
    i = len(seqs[0]) // 2 + pos - 1
    res = [_seq_test(s[i], nt) for s in seqs]
    return res, np.mean(res) * 100


def perc_seq_switch_pos(seqs_1, seqs_2, pos, nt1, nt2):
    i = len(seqs_1[0]) // 2 + pos - 1
    res = [_seq_test(s1[i], nt1) & _seq_test(s2[i], nt2) for s1, s2 in zip(seqs_1, seqs_2)]
    return res, np.mean(res) * 100


def get_alt1_alt2_pos(df):
    alt1_alt2_donor_acceptors = {'alt1': {},
                                 'alt2': {}}
    for _, record in df.iterrows():
        chrom, strand = record.chrom, record.strand
        donor_idx = 0 if record.strand == '+' else 1
        acceptor_idx = 1 - donor_idx
        for alt_type in alt1_alt2_donor_acceptors:
            try:
                donor = (chrom, record[alt_type][donor_idx], strand)
                acceptor = (chrom, record[alt_type][acceptor_idx], strand)
            except TypeError:
                assert record[alt_type] is None
                continue
            if 'donor' not in alt1_alt2_donor_acceptors[alt_type]:
                alt1_alt2_donor_acceptors[alt_type]['donor'] = []
                alt1_alt2_donor_acceptors[alt_type]['acceptor'] = []
            alt1_alt2_donor_acceptors[alt_type]['donor'].append(donor)
            alt1_alt2_donor_acceptors[alt_type]['acceptor'].append(acceptor)
    return alt1_alt2_donor_acceptors


def get_donor_acceptor_seqs_from_df(df, fasta_fn, winsize=6):
    donor_acceptors = get_alt1_alt2_pos(df)
    donor_acceptor_seqs = {}
    with pysam.FastaFile(fasta_fn) as fasta:
        for alt_type in donor_acceptors:
            donor_acceptor_seqs[alt_type] = {}
            for ss_type in donor_acceptors[alt_type]:
                seqs = {}
                for chrom, pos, strand in donor_acceptors[alt_type][ss_type]:
                    seq = fasta.fetch(chrom, pos - winsize, pos + winsize)
                    if strand == '-':
                        seq = rev_comp(seq)
                    seqs[(chrom, pos, strand)] = seq.replace('T', 'U')
                donor_acceptor_seqs[alt_type][ss_type] = seqs
    return donor_acceptors, donor_acceptor_seqs


def get_pairwise_seqs(pos_dict, seq_dict, which='donor'):
    seqs = {'alt1': [], 'alt2': []}
    for alt_type in ['alt1', 'alt2']:
        for idx in pos_dict[alt_type][which]:
            seqs[alt_type].append(seq_dict[alt_type][which][idx])
    return seqs['alt1'], seqs['alt2']


def plot_donor_acceptor_logos(seqs, seq_type='donor', title=None, ax=None):
    w = len(seqs[0]) // 2
    if seq_type == 'donor':
        trim_method = lambda seqs: [s[w - 2: w + 5] for s in seqs]
        xticks = np.arange(7) + 0.5
        xticklabels = ['−2', '−1', '+1', '+2', '+3', '+4', '+5']
    elif seq_type == 'acceptor':
        trim_method = lambda seqs: [s[w - 5: w + 2] for s in seqs]
        xticks = np.arange(7) + 0.5
        xticklabels = ['−5', '−4', '−3', '−2', '−1', '+1', '+2']
    else:
        trim_method = lambda seqs: seqs
        xticks = []
        xticklabels = []
    ax = logo.draw_logo(
        trim_method(seqs),
        alphabet='rna',
        y_format='probability',
        ax=ax,
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks(xticks)
    ax.tick_params(bottom=False)
    ax.set_xticklabels(xticklabels)
    if title is not None:
        ax.set_title(f'{title} (n={len(seqs)})')
    return ax


U5_ORDER = ['AG', 'BG', 'AH', 'BH']
U6_ORDER = ['RAG', 'YAG', 'RBG', 'RAH', 'YAH', 'YBG', 'RBH', 'YBH']

def plot_u5_u6_heatmap(x_seqs, y_seqs=None,
                       x_method=u6_classify,
                       y_method=u5_classify,
                       x_label='U6 class',
                       y_label='U5 class',
                       x_order=U5_ORDER,
                       y_order=U6_ORDER,
                       normalize=False,
                       vmax=None,
                       vmin=0,
                       ax=None):
    if y_seqs is None:
        y_seqs = x_seqs
    x_grp = []
    y_grp = []
    for seq1, seq2 in zip(x_seqs, y_seqs):
        xc = x_method(seq1)
        yc = y_method(seq2)
        x_grp.append(xc)
        y_grp.append(yc)
    contingency_table = pd.crosstab(
        pd.Series(y_grp, name=y_label),
        pd.Series(x_grp, name=x_label),
        normalize=normalize
    )
    contingency_table = (contingency_table
        .reindex(index=x_order, columns=y_order)
        .fillna(0)
    )
    if not normalize:
        contingency_table = contingency_table.astype(int)
        
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(
        contingency_table,
        vmax=vmax,
        vmin=vmin,
        ax=ax,
        cmap='Blues',
        square=True,
        cbar=False,
        annot=True,
        fmt='d' if not normalize else '.2f'
    )
    plt.setp(ax.get_yticklabels(), va='center')
    return ax, y_grp, x_grp


def calculate_switch_nt(record, which='a5'):
    alt1, alt2 = record.alt1, record.alt2
    strand = record.strand
    if which  == 'a5':
        idx = 0 if strand == '+' else 1
    else:
        idx = 1 if strand == '+' else 0
    switch = alt2[idx] - alt1[idx]
    if strand == '-':
        switch = np.negative(switch)
    return switch


def g_test_seqs(seqs_1, seqs_2):
    freqs = np.round(
        [
            logo.calculate_normalised_counts(seqs_1, logo.ALPHABETS['rna']).ravel() * len(seqs_1),
            logo.calculate_normalised_counts(seqs_2, logo.ALPHABETS['rna']).ravel() * len(seqs_2),
        ]
    ).astype('int')
    freqs = freqs[:, freqs.sum(0) > 0]
    return stats.chi2_contingency(freqs, lambda_='log-likelihood')