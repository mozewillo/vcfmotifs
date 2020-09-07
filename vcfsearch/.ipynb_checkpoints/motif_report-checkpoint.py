### for motif --how many detected changes, --where the changes,
### --max score, --avrg score,--extra presence from altered 
### --what TF resp
import pandas as pd

from mtfromvcf.vcfmotifs.vcf_search import find_motif


def motif_frame(*args):
    """initialize data frame - specify what to include in report"""
    motif_frame = pd.DataFrame(columns=['SNP_POS', 'POS_MTF', 'REF -> ALT',
                                        'Ref_score', 'Alt_score', 'Difference'])
    for scr in args:
        if scr == "G_score":
            motif_frame["G_score"] = None
        elif scr == 'ref_diff':
            motif_frame["Ref_diff"] = None
        elif scr == 'alt_diff':
            motif_frame["Alt_diff"] = None
    return motif_frame


def motif_search(record, motif, motif_data=None):
    """single motif for one SNP"""
    if motif_data is None:
        motif_data = motif_frame()
    for snp in record:
        results = find_motif(record, motif, search_threshold=-10, difference_threshold=4, g_threshold=None)
        motif_data = motif_data.append()
    return motif_data


def motif_analyses(record, motif):
    """single motif vcf"""

    pass


def motifs_search(record, motifs):
    """group of motifs -- for vcf """
    pass
