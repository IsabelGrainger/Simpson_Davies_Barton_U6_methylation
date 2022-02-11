import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib_logo as logo
import seaborn as sns

import pysam


RC = str.maketrans('ACGT', 'TGCA')


def rev_comp(seq):
    return seq.translate(RC)[::-1]


IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U',
    'S': 'GC', 'W': 'AU', 'R': 'AG', 'Y': 'CU',
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


def get_donor_acceptor_seqs_from_df(df, fasta_fn, winsize=6):
    alt1_donor_seqs = []
    alt1_acceptor_seqs = []
    alt2_donor_seqs = []
    alt2_acceptor_seqs = []
    with pysam.FastaFile(fasta_fn) as fasta:
        for _, record in df.iterrows():
            donor_idx = 0 if record.strand == '+' else 1
            acceptor_idx = 1 - donor_idx
            if record.alt1 is not None:
                alt1_donor = record.alt1[donor_idx]
                alt1_acceptor = record.alt1[acceptor_idx]
                alt1_donor_seq = fasta.fetch(record.chrom, alt1_donor - winsize, alt1_donor + winsize)
                alt1_acceptor_seq = fasta.fetch(record.chrom, alt1_acceptor - winsize, alt1_acceptor + winsize)
                if record.strand == '-':
                    alt1_donor_seq = rev_comp(alt1_donor_seq)
                    alt1_acceptor_seq = rev_comp(alt1_acceptor_seq)
                alt1_donor_seqs.append(alt1_donor_seq.replace('T', 'U'))
                alt1_acceptor_seqs.append(alt1_acceptor_seq.replace('T', 'U'))
            else:
                alt1_donor_seqs.append(None)
                alt1_acceptor_seqs.append(None)
            if record.alt2 is not None:
                alt2_donor = record.alt2[donor_idx]
                alt2_acceptor = record.alt2[acceptor_idx]
                alt2_donor_seq = fasta.fetch(record.chrom, alt2_donor - winsize, alt2_donor + winsize)
                alt2_acceptor_seq = fasta.fetch(record.chrom, alt2_acceptor - winsize, alt2_acceptor + winsize)
                if record.strand == '-':
                    alt2_donor_seq = rev_comp(alt2_donor_seq)
                    alt2_acceptor_seq = rev_comp(alt2_acceptor_seq)
                alt2_donor_seqs.append(alt2_donor_seq.replace('T', 'U'))
                alt2_acceptor_seqs.append(alt2_acceptor_seq.replace('T', 'U'))
            else:
                alt2_donor_seqs.append(None)
                alt2_acceptor_seqs.append(None)
    return alt1_donor_seqs, alt1_acceptor_seqs, alt2_donor_seqs, alt2_acceptor_seqs


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