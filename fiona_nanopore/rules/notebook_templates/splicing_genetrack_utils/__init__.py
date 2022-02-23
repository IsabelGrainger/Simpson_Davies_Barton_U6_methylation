import random
from collections import defaultdict
import itertools as it

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import patches, gridspec
import seaborn as sns

import pysam
import pyBigWig as pybw

pal = sns.color_palette(['#0072b2', '#d55e00', '#009e73', '#f0e442', '#cc79a7', '#56b4e9', '#e69f00'])

RC = str.maketrans('ACGT', 'TGCA')


def bam_cigar_to_invs(aln):
    invs = []
    start = aln.reference_start
    end = aln.reference_end
    strand = '-' if aln.is_reverse else '+'
    left = start
    right = left
    has_ins = False
    for op, ln in aln.cigar:
        if op in (1, 4, 5):
            # does not consume reference
            continue
        elif op in (0, 2, 7, 8):
            # consume reference but do not add to invs yet
            right += ln
        elif op == 3:
            invs.append([left, right])
            left = right + ln
            right = left
    if right > left:
        invs.append([left, right])
    assert invs[0][0] == start
    assert invs[-1][1] == end
    return start, end, strand, np.array(invs)


def bam_query_iterator(bam, *args, **kwargs):
    strand = kwargs.pop('strand', None)
    if strand is None or strand == '.':
        for aln in bam.fetch(*args, **kwargs):
            yield bam_cigar_to_invs(aln)
    elif strand in '+-':
        is_reverse = strand == '-'
        for aln in bam.fetch(*args, **kwargs):
            if is_reverse == aln.is_reverse:
                yield bam_cigar_to_invs(aln)
    else:
        raise ValueError('strand is not one of +-.')


def random_undersample_bam_query(bam, chrom, start, end, strand=None,
                                 sample_size=100, random_state=None):
    invs = []
    for s, e, _, i in bam_query_iterator(bam, chrom, start, end, strand=strand):
        
        invs.append(i)
    if len(invs) <= sample_size:
        return invs
    else:
        idx = np.arange(len(invs))
        if random_state is None:
            idx = np.random.choice(idx, size=sample_size, replace=False)
        else:
            if isinstance(random_state, int):
                random_state = np.random.RandomState(random_state)
            idx = random_state.choice(idx, size=sample_size, replace=False)
        invs = [invs[i] for i in idx]
        return invs


def get_introns(invs):
    d = [i[0] for i in invs]
    a = [i[1] for i in invs]
    return list(zip(a[:-1], d[1:]))


def get_bam_coverage(bam_fn, chrom, start, end, strand,
                     donor_acceptor_pairs,
                     sample_size=100, random_state=1):
    with pysam.AlignmentFile(bam_fn) as bam:
        invs = random_undersample_bam_query(
            bam, chrom, start, end, strand,
            sample_size, random_state
        )
        n_invs = len(invs)
        introns = [set(get_introns(i)) for i in invs]
        intron_grouped = defaultdict(list)
        for i, v in zip(introns, invs):
            for p in donor_acceptor_pairs:
                if p in i:
                    intron_grouped[p].append(v)
                    break
            else:
                intron_grouped['ungrouped'].append(v)
        # within group sort by start
        for g in intron_grouped.values():
            g.sort(key=lambda invs: invs[0][0] if strand == '+' else invs[-1][1])
        return intron_grouped, n_invs


def rev_comp(seq):
    return seq.translate(RC)[::-1]


def get_range(pos, winsize):
    pos_range = (min(pos) - (winsize + 1), max(pos) + winsize)
    rel_pos = np.array(pos) - min(pos) + winsize - 1
    ncols = pos_range[1] - pos_range[0]
    return pos_range, rel_pos, ncols


def get_donor_acceptor_seq(fasta_fn, chrom, donor_range, acceptor_range, strand):
    with pysam.FastaFile(fasta_fn) as fasta:
        donor_seq = fasta.fetch(chrom, *donor_range)
        acceptor_seq = fasta.fetch(chrom, *acceptor_range)
    if strand == '-':
        donor_seq = rev_comp(donor_seq)
        acceptor_seq = rev_comp(acceptor_seq)
    return donor_seq.replace('T', 'U'), acceptor_seq.replace('T', 'U')


def plot_inv(i, ax, color, intron_color, y, h):
    if len(i) > 1:
        s, e = i[0, 1], i[-1, 0]
        ax.plot([s, e], [y + h / 2, y + h / 2], color=intron_color, zorder=0)
    for start, end in i:
        p = patches.Rectangle((start, y), width=end-start, height=h, color=color, zorder=1)
        ax.add_patch(p)



def plot_invs(grouped_invs, order, donor_range, acceptor_range, axes, color, y_offset=0, h=1):
    y = y_offset
    for p in order:
        g = grouped_invs[p]
        g = sorted(g, key=lambda invs: invs[0][0])
        for i in grouped_invs[p]:
            plot_inv(i, axes[0], color, '#eeeeee', y, h)
            plot_inv(i, axes[1], color, '#eeeeee', y, h)
            y = y + h + h / 10


def set_xticklabel_colors(ax):
    colors = {'A':  pal[2], 'C':  pal[0], 'G': pal[3], 'U': pal[1]}
    for lab  in ax.get_xticklabels():
        try:
            plt.setp(lab, color=colors[lab.get_text()])
        except KeyError:
            continue
            

def plot_donor_acceptor_nanopore(chrom, donors, acceptors, strand,
                                 bam_fns, labels, fasta_fn,
                                 winsize, title, height=1, ypad=10,):
    donor_range, d_rel_pos, donor_ncols = get_range(donors, winsize)
    acceptor_range, a_rel_pos, acceptor_ncols = get_range(acceptors, winsize)

    donor_seq, acceptor_seq = get_donor_acceptor_seq(
        fasta_fn,
        chrom, donor_range, acceptor_range, strand
    )
    n_conds = len(bam_fns)
    
    fig = plt.figure(figsize=((donor_ncols + acceptor_ncols) // 4, 3 * n_conds))
    axes = []
    axes.append(
        plt.subplot2grid((1, donor_ncols + acceptor_ncols), loc=(0, 0), colspan=donor_ncols)
    )
    axes.append(
        plt.subplot2grid((1, donor_ncols + acceptor_ncols),
                         loc=(0, donor_ncols),
                         colspan=acceptor_ncols,
                         sharey=axes[0])
    )

    if strand == '+':
        donor_acceptor_pairs = [(d - 1, a - 1) for d, a in it.product(donors, acceptors)]
    else:
        donor_acceptor_pairs = list(it.product(acceptors, donors))
    order = donor_acceptor_pairs + ['ungrouped']
    y_offset = 0
    handles = []
    for bam_fn, label, colour in zip(bam_fns, labels, pal[:2]):

        invs, n_invs = get_bam_coverage(
            bam_fn, chrom,
            min(*donor_range, *acceptor_range),
            max(*donor_range, *acceptor_range),
            strand, donor_acceptor_pairs,
            sample_size=100, random_state=101
        )
        plot_invs(invs, order, donor_range, acceptor_range, axes, colour, y_offset)
        h = axes[1].fill_between([], [], [], color=colour, label=label)
        handles.append(h)
        y_offset += n_invs * height + ypad

    axes[0].set_ylim(y_offset + ypad, -ypad)

    if strand == '-':
        donor_range = donor_range[::-1]
        acceptor_range = acceptor_range[::-1]

    for ax, range_, seq, pos, hide_side, ss_type in zip(axes,
                                                        [donor_range, acceptor_range],
                                                        [donor_seq, acceptor_seq],
                                                        [donors, acceptors],
                                                        ['left', 'right'],
                                                        ['5\'SS', '3\'SS']):
        ax.set_xlim(range_[0], range_[1] + (1 if strand == '-' else -1))
        ax.tick_params(bottom=False, left=False, labelleft=False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines[hide_side].set_visible(False)
        ax.set_xticks(np.arange(*range_, 1 if strand == '+' else -1)[:-1] + (0.5 if strand == '+' else -0.5))
        ax.set_xticklabels(list(seq)[:-1], fontfamily='monospace', fontsize=25, weight='bold')
        set_xticklabel_colors(ax)
        for i, d in enumerate(pos, 1):
            if len(pos) > 1:
                _label = f'Alt. {ss_type} {i}'
            else:
                _label = ss_type
            d = d - 1 if strand == '+' else d
            ax.plot([d, d], [0, y_offset], ls='--', color='#252525', zorder=4)
            ax.annotate(s=_label, xy=(d, 0), va='bottom', ha='center')
        
    axes[1].legend(handles, labels, title='Genotype', loc=4)
    fig.suptitle(title, y=1.01)
    return fig, axes


def plot_gene_track(record, bam_fns, labels, fasta_fn, winsize=6, title=None):
    if record.event_class == 'A5':
        if record.strand == '+':
            donors = [record.alt1[0] + 1, record.alt2[0] + 1]
            acceptors = [record.alt1[1] + 1]
        else:
            donors = [record.alt1[1], record.alt2[1]]
            acceptors = [record.alt1[0]]
    elif record.event_class == 'RI':
        if record.strand == '+':
            donors = [record.alt2[0] + 1]
            acceptors = [record.alt2[1] + 1]
        else:
            donors = [record.alt2[1]]
            acceptors = [record.alt2[0]]
    elif record.event_class == 'A3':
        if record.strand == '+':
            donors = [record.alt1[0] + 1]
            acceptors = [record.alt1[1] + 1, record.alt2[1] + 1]
        else:
            donors = [record.alt1[1]]
            acceptors = [record.alt1[0], record.alt2[0]]
    elif record.event_class == 'SE':
        raise NotImplementedError()
    else:
        raise NotImplementedError()
    # skip if the distance between donors and acceptors is too big as figures become huge
    if (np.diff(donors) > 120).any() or (np.diff(acceptors) > 120).any():
        raise NotImplementedError()
    else:
        return plot_donor_acceptor_nanopore(
            record.chrom, donors, acceptors, record.strand,
            bam_fns, labels, fasta_fn,
            winsize=winsize, title=title
        )

