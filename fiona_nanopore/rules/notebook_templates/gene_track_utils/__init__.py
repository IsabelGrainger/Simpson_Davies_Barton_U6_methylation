import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import patches, gridspec
import seaborn as sns

import pysam
import pyBigWig as pybw

pal = sns.color_palette(['#0072b2', '#d55e00', '#009e73', '#f0e442', '#cc79a7', '#56b4e9', '#e69f00'])


def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search(f'{attribute} "(.+?)";', gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError(
            f'Could not parse attribute {attribute} '
            f'from GTF with feature type {record[2]}'
        )
    return attr


def flatten_invs(invs):
    flattened = []
    all_invs = iter(sorted(invs))
    inv_start, inv_end = next(all_invs)
    for start, end in all_invs:
        if start <= inv_end:
            inv_end = max(inv_end, end)
        else:
            flattened.append([inv_start, inv_end])
            inv_start, inv_end = start, end
    if not flattened or flattened[-1] != [inv_start, inv_end]:
        flattened.append([inv_start, inv_end])
    return flattened


def get_exons_introns_flank(invs, flanksize=200):
    start = invs[0][0]
    end = invs[-1][1]
    left_flank = [[max(0, start - flanksize), start]]
    right_flank = [[end, end + flanksize]]
    if len(invs) > 1:
        introns = list(zip([i[1] for i in invs[:-1]],
                           [i[0] for i in invs[1:]]))
    else:
        introns = []
    introns = [list(i) for i in introns]
    return invs, introns, left_flank, right_flank


def split_intervals(invs, pos, side='left'):
    invs = np.array(invs)
    idx = np.searchsorted(invs.ravel(), pos)
    split = np.insert(invs.ravel(), idx, [pos, pos]).reshape(-1, 2)
    split_idx = (idx + 1) // 2
    return split[:split_idx].tolist(), split[split_idx:].tolist()


def get_cds_utr_introns_flank(invs, cds_pos, flanksize=200):
    exons, introns, left_flank, right_flank = get_exons_introns_flank(invs, flanksize)
    utr1, cds = split_intervals(exons, cds_pos[0])
    cds, utr2 = split_intervals(cds, cds_pos[1])
    return utr1, cds, utr2, introns, left_flank, right_flank


def gtf_parser(gtf_fn, id_attrb='transcript_id'):
    gtf_invs = {}
    gtf_cds_range = {}
    with open(gtf_fn) as gtf:
        for i, record in enumerate(gtf):
            if record.startswith('#'):
                continue
            record = record.split('\t')
            chrom, _, feat_type, start, end, _, strand = record[:7]
            start = int(start)
            end = int(end)
            transcript_id = get_gtf_attribute(record, id_attrb)
            idx = (chrom, transcript_id, strand)
            if feat_type == 'exon':
                if idx not in gtf_invs:
                    gtf_invs[idx] = []
                gtf_invs[idx].append((start, end))
            elif feat_type == 'CDS':
                if idx not in gtf_cds_range:
                    gtf_cds_range[idx] = [np.inf, -np.inf]
                gtf_cds_range[idx] = [min(start, gtf_cds_range[idx][0]),
                                      max(end, gtf_cds_range[idx][1])]
    for chrom, transcript_id, strand in gtf_invs.keys():
        invs = flatten_invs(gtf_invs[(chrom, transcript_id, strand)])
        try:
            cds_range = gtf_cds_range[(chrom, transcript_id, strand)]
        except KeyError:
            # non-coding, ignore
            continue
        utr5, cds, utr3, introns, upstream, downstream = get_cds_utr_introns_flank(
            invs, cds_range
        )
        if strand == '-':
            utr5, utr3 = utr3, utr5
            upstream, downstream = downstream, upstream
        features = {
            'Upstream': upstream,
            '5\'UTR': utr5,
            'CDS': cds,
            'Introns': introns,
            '3\'UTR': utr3,
            'Downstream': downstream
        }
        yield chrom, transcript_id, strand, features


def count_sites_in_genic_features(bed_fn, gtf_fn):
    feature_counts = {}
    n_records = 0
    with pysam.TabixFile(bed_fn) as tabix:
        for chrom, transcript_id, strand, features in gtf_parser(gtf_fn):
            if chrom in ['Mt', 'Pt']:
                continue
            for feat_type, invs in features.items():
                feat_inv_len = 0
                ov_count = 0
                for inv in invs:
                    if not inv[0] == inv[1]:
                        feat_inv_len += inv[1] - inv[0]
                        for ov in tabix.fetch(chrom, *inv):
                            ov = ov.split('\t')
                            if ov[5] == strand:
                                ov_count += 1
                if feat_inv_len:
                    feature_counts[(feat_type, transcript_id)] = (ov_count / feat_inv_len) * 1000
                else:
                    feature_counts[(feat_type, transcript_id)] = 0
    feature_counts = pd.Series(feature_counts, name='sites_per_kb')
    feature_counts.index = feature_counts.index.rename(names=['feat_type', 'transcript_id'])
    return feature_counts.reset_index()


def get_gene_features_for_coverage_track(gene_id, gtf_fn):
    for chrom, g_id, strand, features in gtf_parser(gtf_fn, 'gene_id'):
        if gene_id == g_id:
            break
    if strand == '+':
        gene_range = [features['Upstream'][0][1], features['Downstream'][0][0]]
    else:
        gene_range = [features['Downstream'][0][1], features['Upstream'][0][0]]
    gene_range_pad = int((gene_range[1] - gene_range[0]) * 0.05)
    features['Range'] = [gene_range[0] - gene_range_pad, gene_range[1] + gene_range_pad]
    return chrom, strand, features


def get_binned(bw, chrom, start, end, strand, binsize):
    vals = bw.values(chrom, start, end, numpy=True)
    vals[np.isnan(vals)] = 0
    if binsize != 1:
        assert binsize > 0
        vals = np.array([sum(vals[i: i + binsize]) for i in np.arange(0, len(vals), binsize)])
    return vals

def plot_bw(bw_fn, colour,
            chrom, start, end, strand,
            ax=None, binsize=10,
            label=None, font_size=10):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    bw = pybw.open(bw_fn)
    vals = get_binned(bw, chrom, start, end, strand, binsize)
    ax.bar(
        np.arange(start, end, binsize),
        vals,
        width=binsize,
        color=colour,
        edgecolor='none',
    )
    if label is not None:
        ax.annotate(
            s=label,
            xy=(0.01, 0.95),
            ha='left',
            va='top',
            xycoords='axes fraction',
            fontsize=font_size
        )
    if strand == '+':
        ax.set_xlim(start, end)
    else:
        ax.set_xlim(end, start)
    ax.set_ylabel(f'miCLIP 5\' coverage ({binsize:d} nt bins)')
    return ax


def plot_inv(i, ax, color, intron_color, y, h):
    if len(i) > 1:
        s, e = i[0][1], i[-1][0]
        ax.plot([s, e], [y, y], color=intron_color, zorder=0)
    for start, end in i:
        p = patches.Rectangle((start, y - h/2), width=end-start, height=h, color=color, zorder=1)
        ax.add_patch(p)


def plot_annot(chrom, features, strand,
               ax=None, utr=False):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 1))
    if not utr:
        for invs, thickness in zip([features['5\'UTR'], features['CDS'], features['3\'UTR']], [1, 1.5, 1]):
            plot_inv(invs, ax, '#252525', '#252525', 1, thickness)
        range_ = features['Range']
    else:
        plot_inv(features['3\'UTR'], ax, '#252525', '#252525', 1, 1)
        stop_codon = features['3\'UTR'][0][0] if strand == '+' else features['3\'UTR'][-1][1]
        stop_codon_inv = [[stop_codon - 20, stop_codon]] if strand == '+' else [[stop_codon, stop_codon + 20]]
        plot_inv(stop_codon_inv, ax, '#252525', '#252525', 1, 1.5)
        range_ = [features['3\'UTR'][0][0] - 20, features['3\'UTR'][-1][1] + 20]
    ax.set_ylim(-1.5, 2)
    if strand == '+':
        ax.set_xlim(*range_)
    else:
        ax.set_xlim(*range_[::-1])
    ax.set_axis_off()
    ax.annotate(s='Araport11 annotation', xy=(0.99, 0.05), ha='right', va='bottom', xycoords='axes fraction')
    return ax


def get_peaks(tabix_fn, chrom, start, end, strand):
    peaks = []
    try:
        with pysam.TabixFile(tabix_fn) as f:
            pass
        tmp_bed = tabix_fn
    except OSError:
        tmp_bed = os.path.join(os.environ['TMPDIR'], 'tmp.bed.gz')
        pysam.tabix_compress(tabix_fn, tmp_bed, force=True)
        pysam.tabix_index(tmp_bed, force=True, preset='bed')

    with pysam.TabixFile(tmp_bed) as f:
        for record in f.fetch(chrom, start, end):
            record = record.split()
            if record[5] == strand:
                range_ = [int(record[1]), int(record[2])]
                peaks.append(range_)
    return peaks


def plot_peaks(tabix_fn, chrom, plot_start, plot_end, strand, color, ax=None, label=None):
    peaks = get_peaks(tabix_fn, chrom, plot_start, plot_end, strand)
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 1))
    for start, end in peaks:
        p = patches.Rectangle((start, 0.5), width=end-start, height=1, color=color, zorder=1)
        ax.add_patch(p)
    ax.set_ylim(-1, 2)
    if strand == '+':
        ax.set_xlim(plot_start, plot_end)
    else:
        ax.set_xlim(plot_end, plot_start)
    if label is not None:
        ax.annotate(s=label, xy=(0.99, 0.05), ha='right', va='bottom', xycoords='axes fraction')
    ax.set_axis_off()
    return ax


def plot_yanoncomp_miclip_gene_or_utr(gene_id, chrom, strand, features, grid, cloc, cspan, utr,
                                      miclip_cov_fns, miclip_peaks_fn, der_sites_fn, yanocomp_fns):
    miclip_cov_ax = plt.subplot2grid(grid, loc=(0, cloc), rowspan=5, colspan=cspan,)
    miclip_peak_ax = plt.subplot2grid(grid, loc=(5, cloc), rowspan=1, colspan=cspan, sharex=miclip_cov_ax)
    annot_ax = plt.subplot2grid(grid, loc=(6, cloc), rowspan=1, colspan=cspan, sharex=miclip_cov_ax)
    der_site_ax = plt.subplot2grid(grid, loc=(7, cloc), rowspan=1, colspan=cspan, sharex=miclip_cov_ax)
    yanocomp_axes = [plt.subplot2grid(grid, loc=(8 + i, cloc), rowspan=1, colspan=cspan, sharex=miclip_cov_ax)
                     for i in range(len(yanocomp_fns))]
    if utr:
        range_ = [features['3\'UTR'][0][0] - 20, features['3\'UTR'][-1][1] + 20]
    else:
        range_ =  features['Range']
    plot_bw(
        miclip_cov_fns[strand == '-'],
        pal[0], chrom, *range_, strand,
        ax=miclip_cov_ax, binsize=1 if utr else 10,
    )
    plot_peaks(
        miclip_peaks_fn,
        chrom, *range_, strand, pal[0], ax=miclip_peak_ax, label='miCLIP peaks'
    )
    plot_annot(
        chrom, features, strand, ax=annot_ax, utr=utr
    )
    plot_peaks(
        der_sites_fn,
        chrom, *range_, strand, pal[1], ax=der_site_ax, label='VIR-dependent differential error sites'
    )
    for (label, fn), ax, c in zip(yanocomp_fns.items(), yanocomp_axes, pal[2:]):
        plot_peaks(
            fn, chrom, *range_, strand, c, ax=ax, label=label
        )
    # convert miclip xaxis to gene track bar
    if not utr:
        miclip_cov_ax.xaxis.tick_top()
        miclip_cov_ax.spines['bottom'].set_visible(False)
        miclip_cov_ax.spines['right'].set_visible(False)
        miclip_cov_ax.set_xlabel(f'{gene_id} (Chr{chrom}, {strand}ve strand)')
        miclip_cov_ax.xaxis.set_label_position('top')
        miclip_cov_ax.set_xticks([
            x for x in miclip_cov_ax.get_xticks() if not x % 1000
        ])
        miclip_cov_ax.set_xticklabels([
            f'{x / 1000:.0f} kb' for x in miclip_cov_ax.get_xticks()
        ])
        if strand == '+':
            miclip_cov_ax.set_xlim(*range_)
        else:
            miclip_cov_ax.set_xlim(*range_[::-1])
    else:
        miclip_cov_ax.set_xticks([])
        miclip_cov_ax.spines['top'].set_visible(False)
        miclip_cov_ax.spines['bottom'].set_visible(False)
        miclip_cov_ax.spines['right'].set_visible(False)

def plot_yanocomp_miclip(gene_id, miclip_cov_fns, miclip_peaks_fn, der_sites_fn, yanocomp_fns, gtf_fn):
    # create figure
    fig = plt.figure(figsize=(16, 10))
    grid = (8 + len(yanocomp_fns), 3)
    chrom, strand, features = get_gene_features_for_coverage_track(gene_id, gtf_fn)
    # whole gene
    plot_yanoncomp_miclip_gene_or_utr(
        gene_id, chrom, strand, features, grid, 0, 2, False,
        miclip_cov_fns, miclip_peaks_fn, der_sites_fn, yanocomp_fns
    )
    # 3utr only
    plot_yanoncomp_miclip_gene_or_utr(
        gene_id, chrom, strand, features, grid, 2, 1, True,
        miclip_cov_fns, miclip_peaks_fn, der_sites_fn, yanocomp_fns
    )
    return fig