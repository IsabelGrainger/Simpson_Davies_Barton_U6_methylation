import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import patches, gridspec
import seaborn as sns

import pysam
import pyBigWig as pybw


pal = sns.color_palette(['#0072b2', '#d55e00', '#009e73', '#f0e442', '#cc79a7', '#56b4e9', '#e69f00'])


def get_coverage(bw, chrom, start, end, strand):
    vals = bw.values(chrom, start, end, numpy=True)
    if strand == '-':
        vals = vals[::-1]
    vals[np.isnan(vals)] = 0
    return vals


def get_donor_acceptor_coverage(bw_fn, chrom, donor_range, acceptor_range, strand):
    bw = pybw.open(bw_fn)
    donor_cov = get_coverage(bw, chrom, *donor_range, strand)
    acceptor_cov = get_coverage(bw, chrom, *acceptor_range, strand)
    bw.close()
    return donor_cov, acceptor_cov


RC = str.maketrans('ACGT', 'TGCA')


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


def set_xticklabel_colors(ax):
    colors = {'A':  pal[2], 'C':  pal[0], 'G': pal[3], 'U': pal[1]}
    for lab  in ax.get_xticklabels():
        try:
            plt.setp(lab, color=colors[lab.get_text()])
        except KeyError:
            continue


def plot_psi(psi_melt, event_id_col='index', ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    ax = sns.boxplot(
        x='temp',
        y='psi',
        hue='geno',
        boxprops=dict(alpha=0.2),
        whiskerprops=dict(alpha=0.2),
        medianprops=dict(alpha=0.2),
        capprops=dict(alpha=0.2),
        showfliers=False,
        order=[4, 12, 20, 28],
        dodge=False,
        data=psi_melt,
        ax=ax
    )
    ax = sns.stripplot(
        x='temp',
        y='psi',
        hue='geno',
        order=[4, 12, 20, 28],
        hue_order=['col0', 'fio1'],
        data=psi_melt,
        jitter=True,
        ax=ax
    )
    ax.set_ylim(0, 1)
    ax.set_ylabel('Fraction spliced in')
    ax.set_xlabel('Temperature (°C)')
    ax.legend_.remove()
    h1 = ax.scatter([], [], color=pal[0])
    h2 = ax.scatter([], [], color=pal[1])
    ax.legend([h1, h2], ['Col-0', 'fio1-2'])
    return ax


def plot_donor_acceptor_gene_track_and_boxplot(chrom, donors, acceptors, strand,
                                               bw_fns, labels, fasta_fn, psi,
                                               winsize=6, title=None,):

    donor_range, d_rel_pos, donor_ncols = get_range(donors, winsize)
    acceptor_range, a_rel_pos, acceptor_ncols = get_range(acceptors, winsize)

    donor_seq, acceptor_seq = get_donor_acceptor_seq(
        fasta_fn,
        chrom, donor_range, acceptor_range, strand
    )

    fig = plt.figure(figsize=((donor_ncols + acceptor_ncols) // 4 + 6, 5))
    axes = []
    axes.append(
        plt.subplot2grid((1, donor_ncols + acceptor_ncols + 24), loc=(0, 0), colspan=donor_ncols)
    )
    axes.append(
        plt.subplot2grid((1, donor_ncols + acceptor_ncols + 24),
                         loc=(0, donor_ncols),
                         colspan=acceptor_ncols,
                         sharey=axes[0])
    )
    axes.append(
        plt.subplot2grid((1, donor_ncols + acceptor_ncols + 24),
                         loc=(0, donor_ncols + acceptor_ncols + 4),
                         colspan=20)
    )
    for i, (bw_fn, lab, c) in enumerate(zip(bw_fns, labels, pal)):
        dcov, acov = get_donor_acceptor_coverage(bw_fn, chrom, donor_range, acceptor_range, strand)
        norm = max(dcov[d_rel_pos].max(), acov[a_rel_pos + 1].max())
        for ax, cov, seq, pos, splice_type in zip(axes[:2],
                                                  [dcov, acov],
                                                  [donor_seq, acceptor_seq],
                                                  [d_rel_pos, a_rel_pos],
                                                  ['5\'SS', '3\'SS']):
            ax.step(
                np.arange(len(cov)), cov / norm,
                where='post', color=c,
                lw=3, alpha=0.85
            )
            ax.fill_between(
                np.arange(len(cov)), cov / norm,
                step='post', color=c, alpha=0.25,
                label=lab,
            )
            if i == 0:
                ax.set_xticks(np.arange(len(cov) - 1) + 0.5)
                ax.set_xticklabels(list(seq)[:-1], fontfamily='monospace', fontsize=25, weight='bold')
                set_xticklabel_colors(ax)
                for j, d in enumerate(pos, 1):
                    splice_label = f'Alt. {splice_type} {j}' if len(pos) > 1 else splice_type
                    ax.plot([d + 1, d + 1], [0, 1 + 0.1 * j], ls='--', color='#252525', zorder=-1)
                    ax.annotate(text=splice_label, xy=(d + 1, 1.05 + 0.1 * j), va='bottom', ha='center')
            
            
            
    axes[0].tick_params(bottom=False)
    axes[1].tick_params(left=False, bottom=False, labelleft=False)
    fig.suptitle(title, y=1.01)
    axes[0].set_ylim(0, 1.35)
    axes[1].legend(title='Genotype', loc=4)
    axes[0].set_ylabel('Normalised coverage')

    plot_psi(psi, ax=axes[2])
    
    return axes


def plot_skipped_exon_and_boxplot(chrom, upstream_donor, downstream_donor,
                                  upstream_acceptor, downstream_acceptor,
                                  strand, bw_fns, labels, fasta_fn, psi,
                                  winsize=6, title=None,):

    udonor_range, ud_rel_pos, udonor_ncols = get_range([upstream_donor,], winsize)
    ddonor_range, dd_rel_pos, ddonor_ncols = get_range([downstream_donor,], winsize)
    uacceptor_range, ua_rel_pos, uacceptor_ncols = get_range([upstream_acceptor], winsize)
    dacceptor_range, da_rel_pos, dacceptor_ncols = get_range([downstream_acceptor], winsize)

    udonor_seq, uacceptor_seq = get_donor_acceptor_seq(
        fasta_fn,
        chrom, udonor_range, uacceptor_range, strand
    )
    ddonor_seq, dacceptor_seq = get_donor_acceptor_seq(
        fasta_fn,
        chrom, ddonor_range, dacceptor_range, strand
    )

    total_genetrack_cols = udonor_ncols + uacceptor_ncols + ddonor_ncols + dacceptor_ncols
    fig = plt.figure(figsize=(total_genetrack_cols // 4 + 6, 5))
    axes = []
    axes.append(
        plt.subplot2grid((1, total_genetrack_cols + 24), loc=(0, 0), colspan=udonor_ncols)
    )
    axes.append(
        plt.subplot2grid((1, total_genetrack_cols + 24),
                         loc=(0, udonor_ncols),
                         colspan=uacceptor_ncols,
                         sharey=axes[0])
    )
    axes.append(
        plt.subplot2grid((1, total_genetrack_cols + 24),
                         loc=(0, udonor_ncols + uacceptor_ncols),
                         colspan=ddonor_ncols,
                         sharey=axes[0])
    )
    axes.append(
        plt.subplot2grid((1, total_genetrack_cols + 24),
                         loc=(0, udonor_ncols + uacceptor_ncols + ddonor_ncols),
                         colspan=dacceptor_ncols,
                         sharey=axes[0])
    )
    axes.append(
        plt.subplot2grid((1, total_genetrack_cols + 24),
                         loc=(0, total_genetrack_cols + 4),
                         colspan=20)
    )
    for i, (bw_fn, lab, c) in enumerate(zip(bw_fns, labels, pal)):
        udcov, uacov = get_donor_acceptor_coverage(bw_fn, chrom, udonor_range, uacceptor_range, strand)
        ddcov, dacov = get_donor_acceptor_coverage(bw_fn, chrom, ddonor_range, dacceptor_range, strand)
        norm = max(udcov[ud_rel_pos].max(), dacov[da_rel_pos + 1].max())
        for ax, cov, seq, pos, splice_type in zip(axes[:4],
                                                  [udcov, uacov, ddcov, dacov],
                                                  [udonor_seq, uacceptor_seq, ddonor_seq, dacceptor_seq],
                                                  [ud_rel_pos, ua_rel_pos, dd_rel_pos, da_rel_pos],
                                                  ['5\'SS', '3\'SS', '5\'SS', '3\'SS']):
            ax.step(
                np.arange(len(cov)), cov / norm,
                where='post', color=c,
                lw=3, alpha=0.85
            )
            ax.fill_between(
                np.arange(len(cov)), cov / norm,
                step='post', color=c, alpha=0.25,
                label=lab,
            )
            if i == 0:
                ax.set_xticks(np.arange(len(cov) - 1) + 0.5)
                ax.set_xticklabels(list(seq)[:-1], fontfamily='monospace', fontsize=25, weight='bold')
                set_xticklabel_colors(ax)
                for j, d in enumerate(pos, 1):
                    splice_label = f'Alt. {splice_type} {j}' if len(pos) > 1 else splice_type
                    ax.plot([d + 1, d + 1], [0, 1 + 0.1 * j], ls='--', color='#252525', zorder=-1)
                    ax.annotate(text=splice_label, xy=(d + 1, 1.05 + 0.1 * j), va='bottom', ha='center')
            
            
            
    axes[0].tick_params(bottom=False)
    axes[1].tick_params(left=False, bottom=False, labelleft=False)
    axes[2].tick_params(left=False, bottom=False, labelleft=False)
    axes[3].tick_params(left=False, bottom=False, labelleft=False)
    fig.suptitle(title, y=1.01)
    axes[0].set_ylim(0, 1.35)
    axes[3].legend(title='Genotype', loc=4)
    axes[0].set_ylabel('Normalised coverage')

    plot_psi(psi, ax=axes[4])
    
    return axes


def plot_gene_track(record, event_psi, bw_fns, labels, fasta_fn, winsize=6, title=None):
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
        if record.strand == '+':
            udonor = record.alt2[0] + 1
            ddonor = record.alt1[0] + 1
            uacceptor = record.alt1[1] + 1
            dacceptor = record.alt2[1] + 1
        if record.strand == '-':
            udonor = record.alt2[1]
            ddonor = record.alt1[1]
            uacceptor = record.alt1[0]
            dacceptor = record.alt2[0]
        # special case for skipped exons requiring different plot type
        return plot_skipped_exon_and_boxplot(
            record.chrom, udonor, ddonor, uacceptor, dacceptor,
            record.strand,
            bw_fns, labels, fasta_fn, event_psi,
            winsize=winsize, title=title
        )
    else:
        raise NotImplementedError()
    # skip if the distance between donors and acceptors is too big as figures become huge
    if (np.diff(donors) > 120).any() or (np.diff(acceptors) > 120).any():
        raise NotImplementedError()
    else:
        return plot_donor_acceptor_gene_track_and_boxplot(
            record.chrom, donors, acceptors, record.strand,
            bw_fns, labels, fasta_fn, event_psi,
            winsize=winsize, title=title
        )

