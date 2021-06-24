"""Create CSV report for searches"

from vcfmotifs.vcfsearch import read_in
import pandas as pd


def start_report(filename, mtf, search_ref_threshold, clinvar=None, **defined_filters):
    """ create file with meta data """
    filename = read_in.name_csv(filename)
    info = "#" + "motif name : " + mtf.name + "motif_id : " + mtf.matrix_id + "\t" + "# motif_consensus_sequence : " + \
           str(mtf.consensus) + "\n" + "# search_ref_threshold : " + str(search_ref_threshold) + "\n" \
           + "# applied filters : " + str(defined_filters) + "\n"
    with open(filename, 'a') as r:
        r.write(info)
        r.write("variant_ID,chrom,var_pos,polymorphism,refseq -> altseq,ref_score,alt_score,motif_loc ")
        if clinvar is not None:
            r.write(',ClinVar')
        r.write('\n')
    return filename


def add_to_report(filename, result, clinvar):
    """ add one record to report"""
    with open(filename, 'a') as mr:
        minfo = str(result['variant_ID']) + "," + str(result['chrom']) + "," + str(result['var_pos']) + "," + \
                str(result["polymorphism"]) + "," + str(result["sequences"]) + "," + str(result["ref_score"]) \
                + "," + str(result["alt_score"]) + "," + str(result["motif_loc"])
        if clinvar is not None:
            minfo += ',' + clinvar[0]
        minfo += '\n'
        mr.write(minfo)


def vcf_frame(*args):
    """ create vcf pandas DataFrame"""
    motif_frame = pd.DataFrame(
        columns=['variant_ID','chrom', 'var_pos', 'matrix_id', 'polymorphism', 'altered_seq', 'ref_score', 'alt_score', 'motif_loc'])
    motif_frame.set_index('matrix_id', inplace=True)
    return motif_frame
