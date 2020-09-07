""" single change research"""

from vcfmotifs.vcfsearch.read_in import name_csv
from vcfmotifs.vcfsearch.score import estimate_threshold
from vcfmotifs.vcfsearch.vcf_search import find_motif


class simpleSNP:
    """ simulate vcf record """
    def __init__(self, CHROM, POS, REF, ALT):
        self.CHROM = CHROM
        self.POS = POS
        self.REF = REF
        self.ALT = ALT


def single_change(recordLike, motifs_readin, resultsfile, fastafile=None, search_ref_threshold='estimate_threshold',
                  custom_function=None, search_alt_threshold=-20, pseudocounts=0.2,
                  background={'A': 0.3, 'C': 0.2, 'G': 0.2, 'T': 0.3}, **filters):
    """a single variation or list of neighboring variations,  """
    if type(recordLike) == list:
        position = []
        ref = ""
        alt = ""
        for record in recordLike:
            position.append(record.POS)
            ref.append(record.REF)
            alt.append(record.ALT)
        polymorphism = (recordLike[0].CHROM, position, str(ref) + " -> " + str(alt))
    else:
        polymorphism = (recordLike.CHROM, recordLike.POS, str(recordLike.REF) + " -> " + str(recordLike.ALT))
    start_info = "# polymorphism: " + polymorphism + "\n" + "# thresholds: ref: " + str(search_ref_threshold) + \
                 " alt: " + str(search_alt_threshold) + "\n" + "# filters: " + str(filters) + "\n"

    resultsfile = name_csv(resultsfile)
    with open(resultsfile, 'a') as resfile:
        resfile.write(start_info)
        resfile.write("matrix_id,'refseq -> altseq,ref_score,alt_score,motif_loc \n")

    for mtf in motifs_readin:
        if search_ref_threshold == 'estimate_threshold':
            search_ref_threshold = estimate_threshold(mtf)
        elif search_ref_threshold == 'custom_threshold_function':
            try:
                custom_function(mtf)
            except ValueError:
                print("Define custom function that takes only motif as argument")

        results = find_motif(recordLike, mtf, search_ref_threshold, search_alt_threshold, pseudocounts=pseudocounts,
                             background=background, fastafile=fastafile, **filters)
        with open(resultsfile, 'a') as resfile:
            for res in results:
                result, diff = res
                minfo = mtf.matrix_id + "," + str(result["sequences"]) + "," + str(result["ref_score"]) \
                        + "," + str(result["alt_score"]) + "," + str(result["motif_loc"]) + "\n"
                resfile.write(minfo)


def single_search(snp_chrom, snp_position, snp_ref, snp_alt, motifs_readin, resultsfile, fastafile=None,
                  search_ref_threshold='estimate_threshold', custom_function=None, search_alt_threshold=-20,
                  background={'A': 0.3, 'C': 0.2, 'G': 0.2, 'T': 0.3}, **filters):
    """single SNP analysis, no vcf file required """
    record = simpleSNP(snp_chrom, snp_position, snp_ref, snp_alt)
    return single_change(record, motifs_readin, resultsfile, fastafile=fastafile,
                         search_ref_threshold=search_ref_threshold,
                         custom_function=custom_function, search_alt_threshold=search_alt_threshold,
                         background=background, **filters)