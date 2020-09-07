import os

import vcf
from seqseek import Chromosome, BUILD38
import pandas as pd

from mtfromvcf.vcfmotifs.filters import G_score, abs_difference, apply_filters, make_filter
from mtfromvcf.vcfmotifs.read_in import read_fasta, multipleSNP, motif_handle, vcf_handle
from mtfromvcf.vcfmotifs.score import ref_score, alt_score, alt_scope, pssm


def snps_position(records):
    """get all positions if more snps in motif scope"""
    if type(records) == list:
        # for snps in the same motif
        position = []
        for snps in records:
            position.append(snps.POS, snps.CHROM[3:])
    else:
        position = (records.POS, records.CHROM[3:])

    return position


def motif_scope(record, motif_len, hg38=False, fastafile=None):
    """return part of reference genome sequence of analysed scope"""
    if type(record) == list:
        locus, chrom = record[0].POS
        diff = record[-1].POS - record[0].POS
        stop = locus + motif_len + diff
    else:
        locus, chrom = record.POS, record.CHROM
        stop = locus + motif_len
    start = locus - motif_len
    if hg38:
        scope = Chromosome(chrom, assembly=BUILD38).sequence(start, stop)
    elif fastafile is not None:
        scope = read_fasta(fastafile, start, stop)
    else:
        raise ValueError("No reference sequence given. Choose hs38 assembly or upload fasta file.")
    return scope


def find_motif(record, motif, search_ref_threshold=6, search_alt_threshold=-20, fastafile=None, **defined_filters):
    """find motif in sequence from given vcf record"""
    if type(record) == list:
        # for snps in the same motif
        position = []
        for snps in record:
            position.append(record.POS, record.CHROM[3:])
    else:
        position = (record.POS, record.CHROM[3:])

    motif_len = motif.length
    ref_seq = motif_scope(record, motif_len, fastafile=fastafile)
    print(ref_seq)
    motif_pssm = pssm(motif)
    # get scores
    ref = ref_score(motif_pssm, ref_seq, search_ref_threshold)
    alt_seq = alt_scope(record, ref_seq, motif_len)
    alt = alt_score(alt_seq, motif_pssm, search_alt_threshold)
    # basic difference
    print(ref_seq)
    # differ = abs_difference(ref, alt, **defined_filters)
    filtered = apply_filters((ref, alt), **defined_filters)

    results = {'matrix_id': motif.matrix_id, 'altered_seq': alt_seq, 'ref_score': ref, 'alt_score': alt}
    print("rezultat", results)
    return results


def vcf_frame(*args):
    motif_frame = pd.DataFrame(columns=['matrix_id'])
    motif_frame.set_index('matrix_id', inplace=True)
    return motif_frame


def research_motifs(motifs_readin, vcf_readin, multiple=False, search_ref_threshold=6, search_alt_threshold=-20,
                    fastafile=None, results_key='difference', **defined_filters):
    """ For every motif iterate thought whole vcf -- enable multiple function"""
    vcfFrame = vcf_frame()
    res = []
    for mtf in motifs_readin:
        mtf_len = len(mtf)
        if multiple:
            vcf_readin = multipleSNP(vcf_readin, mtf_len)
        for snp in vcf_readin:
            print("snp", snp)
            make_filter(**defined_filters)
            rc = find_motif(snp, mtf, search_ref_threshold, search_alt_threshold, fastafile=fastafile)
            res.append(rc)
    return res


#
# print(os.path.isdir('/home/aleksandra/Desktop/lic_py/dmel_data/moti'))
mtfs = motif_handle('/home/aleksandra/Desktop/lic_py/dmel_data/moti')
vcfs = vcf_handle('/home/aleksandra/Desktop/lic_py/dmel_data/chr_X.vcf.gz')

#
# for m in mtfs:
#     print(m.consensus)


# def opened():
#     with open('/home/aleksandra/Desktop/lic_py/dmel_data/chr_X.vcf.gz', 'rb') as f:
#         vcfr = vcf.Reader(f, compressed=True)
#         try:
#             for v in vcfr:
#                 yield v
#         except IndexError:
#             pass

# vcfr = opened()
l = 0
# for f in vcfs:
#     l+=1
#     print(f)
# print(l)

res = research_motifs(mtfs, vcfs, fastafile='/home/aleksandra/Desktop/lic_py/dmel_data/dmel_chrx.fasta', difference=4)

print(res)
