import pandas as pd
from vcfmotifs.vcfsearch.score import gscore


def goodness_of_fit(ref_seq, alt_seq, motif, threshold):
    gscr = gscore(ref_seq, alt_seq, motif)
    return gscr > threshold


def max_score(ref_seq, alt_seq, threshold):
    """Specify threshold for motif match"""
    ind, max_score = max(enumerate((ref_seq, alt_seq)))
    return max_score > threshold


def min_score(ref_scr, alt_scr, threshold):
    """ How small score need to be to consider as disturbing linkage change"""
    minimal = min(ref_scr, alt_scr)
    return minimal < threshold


def refmotif_differ(ref_seq, alt_seq, threshold):
    """How relevant change may be ref -> alt,
    if motif occurs in referenced seq"""
    diff = ref_seq - alt_seq
    return diff > threshold


def altmotif_differ(ref_seq, alt_seq, threshold):
    """How relevant change may be ref -> alt,
    if motif occurs in altered seq"""
    diff = alt_seq - ref_seq
    return diff > threshold


def abs_difference(ref_scr, alt_scr, threshold=4):
    """Abs difference - change ref <-> alt,
    returns diff and which variant has positive score"""
    if ref_scr >= 0 or alt_scr >= 0:
        diff = abs(ref_scr - alt_scr)
        # ind, val = max(enumerate((ref_scr, alt_scr)))
        return diff > threshold


def difference(ref_scr, alt_scr):
    diff = abs(ref_scr - alt_scr)
    return diff


def apply_filters(scores, motif, **defined_thresholds):  # defined thresholds dictionary
    """apply given filers on results dataset
    Given filters: min_score, max_score, difference, ref_diff, alt_diff"""
    ref_scr, alt_scr = scores
    filters_results = {}

    for filters, threshold in defined_thresholds.items():
        if filters == "min_score":
            filter = min_score(ref_scr, alt_scr, threshold)
            filters_results[filters] = filter
        if filters == "max_score":
            filter = max_score(ref_scr, alt_scr, threshold)
            filters_results[filters] = filter
        if filters == "difference":
            filter = abs_difference(ref_scr, alt_scr, threshold)
            filters_results[filters] = filter
        if filters == "ref_diff":
            filter = refmotif_differ(ref_scr, alt_scr, threshold)
            filters_results[filters] = filter
        if filters == "alt_diff":
            filter = altmotif_differ(ref_scr, alt_scr, threshold)
            filters_results[filters] = filter
        if filters == "g_score":
            filter = goodness_of_fit(ref_scr, alt_scr, motif, threshold)
    if all(filters_results.values()):
        diff = difference(ref_scr, alt_scr)
        return diff


def make_filter(min_score=None, max_score=None, difference=4, ref_differ=None, alt_differ=None,
                results_key='difference', g_score=None):
    """create filter -  dictionary"""
    filters = {'results_key': results_key, 'difference': difference, 'alt_differ': alt_differ,
               'ref_differ': ref_differ, 'min_score': min_score, 'max_score': max_score, 'g_score': g_score}
    return {k: v for k, v in filters.items() if v is not None}


def filter_data(res_frame, **defined_filters):
    """ separate function to apply filters on results"""
    for filt, threshold in defined_filters.items():
        if filt == "difference":
            condition = abs(res_frame["ref_score"] - res_frame["alt_score"]) > threshold
            res_frame = res_frame[condition]
        if filt == "min_score":
            condition1 = res_frame["ref_score"] < threshold
            condition2 = res_frame["alt_score"] < threshold
            res_frame = res_frame[condition1 | condition2]
        if filt == "max_score":
            condition1 = res_frame["ref_score"] > threshold
            condition2 = res_frame["alt_score"] > threshold
            res_frame = res_frame[condition1 | condition2]
        if filt == "ref_diff":
            condition = (res_frame["ref_score"] - res_frame["alt_score"]) > threshold
            res_frame = res_frame[condition]
        if filt == "alt_diff":
            condition = (res_frame["alt_score"] - res_frame["ref_score"]) > threshold
            res_frame = res_frame[condition]
        if filt == "chrom": # get scores only from chosen chromosome
            res_frame = res_frame.loc[res_frame['chrom'] == threshold]
    return res_frame


def filter_results(mtf_csvfile, outfile=None, **filters):
    """readin as pandas DataFrame in batches filter by columns operations on alt and ref seq"""
    mtfdata = pd.read_csv(mtf_csvfile, skiprows=3)
    if len(mtfdata.columns) == 9:
        mtfdata.columns = ['ID', 'chrom', 'positions', 'polymorphism', 'sequences', 'ref_score', 'alt_score', 'motif_loc', 'clinvar']
    else:
        mtfdata.columns = ['ID', 'chrom', 'positions', 'polymorphism', 'sequences', 'ref_score', 'alt_score',
                           'motif_loc']
    mtfdata = filter_data(mtfdata, **filters)
    mtfdata.to_csv(outfile)
    return mtfdata


