import os

from vcfmotifs.vcfsearch.csv_report import start_report, add_to_report
from vcfmotifs.vcfsearch.filters import apply_filters, make_filter
from vcfmotifs.vcfsearch import scopes, html_reports
from vcfmotifs.vcfsearch import score
from vcfmotifs.vcfsearch import read_in
from vcfmotifs.vcfsearch.read_in import clinvar_url
from vcfmotifs.vcfsearch.score import estimate_backgroud
import numpy as np


def find_motif(record, motif, refseq,  search_ref_threshold=6, search_alt_threshold=-20, pseudocounts=0.2,
               background={'A': 0.3, 'C': 0.2, 'G': 0.2, 'T': 0.3},  **defined_filters):
    """PWM motif match in sequence scope from given vcf record"""
    ref_seq = scopes.motif_scope(record, motif.length, refseq=refseq)
    if ref_seq == None:
        return None
    motif_pssm = score.pssm(motif, pseudocounts=pseudocounts, background=background)
    # get scores
    ref = score.ref_score(motif_pssm, ref_seq, search_ref_threshold)
    alt_seq = scopes.alt_scope(record, ref_seq, motif.length)
    alt = score.alt_score(alt_seq, motif_pssm, search_alt_threshold)

    # pairing scores in ref and alt scopes, applying further filters
    pairs = score.pair_scores(ref, alt)
    for pos, scr in pairs:
        scr = np.round(scr, 2)
        # applying defined filters
        filtered = apply_filters((scr[0], scr[1]), motif, **defined_filters)
        if filtered is not None:
            result = create_result(record=record, motif=motif, position=pos, score=scr, ref_seq=ref_seq, alt_seq=alt_seq)
            yield result, filtered


def create_result(record, motif, position, score, ref_seq, alt_seq):
    end = position + motif.length
    if type(record) == list:
        end = position + motif.length
        mloc = int(record[0].POS - motif.length + 1 + position)
        positions = []
        ref = []
        alt = []
        varID = []
        for r in record:
            posin = str(r.POS)
            if mloc <= int(posin) <= mloc + motif.length:
                positions.append(posin)
                ref.append(r.REF)
                alt.append(str(r.ALT))
                varID.append(r.ID)
        positions = '/ '.join(positions)
        alt = '/ '.join(alt)
        ref = '/ '.join(ref)
        varID = '/ '.join(varID)
        polymorphism =(str(ref) + " -> " + str(alt))
        chrom =  record[0].CHROM
    else:
        polymorphism = str(str(record.REF) + " -> " + str(record.ALT))
        mloc = int(record.POS - motif.length + 1 + position)
        varID = record.ID
        positions = record.POS
        chrom = record.CHROM
    result = {'motif_loc': mloc, 'sequences': str(ref_seq)[position:end] + " -> " + str(alt_seq)[position:end],
              'ref_score': score[0], 'alt_score': score[1], 'polymorphism': polymorphism, 'variant_ID': varID,
              'chrom' : chrom, 'var_pos': positions}
    return result


def vcf_search(vcf_path, vcf_splice=None, motifs_path=None, motifs_readin=None, refseq='hg38', multiple=False,
               search_ref_threshold='estimate_threshold',g_score=0.8, pseudocounts="calculate_pseudocounts",
               background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}, custom_function=None,
               search_alt_threshold=-20, mtf_report=None, html_report=None, clinvar=None, **defined_filters):
    """Search for given motifs in reference and altered sequence
       iterating through whole vcf, recommended to enable multiple function """
    if motifs_path is not None:
        motifs_readin = read_in.motif_handle(motif_path=motifs_path)
    filters = make_filter(**defined_filters)

    for mtf in motifs_readin:
        # estimated threshold for motif discovery if not given
        if search_ref_threshold == 'estimate_threshold':
            search_ref_threshold = score.estimate_threshold(mtf)
        elif search_ref_threshold == 'custom_threshold_function':
            try:
                custom_function(mtf)
            except ValueError:
                print("Define custom function that takes only motif as argument")
        motif_search(vcf_path=vcf_path, vcf_splice=vcf_splice, motif=mtf, refseq=refseq, multiple=multiple,
                     search_ref_threshold=search_ref_threshold, search_alt_threshold=search_alt_threshold, g_score=g_score,
                     pseudocounts=pseudocounts, background=background, mtf_report=mtf_report, html_report=html_report,
                     clinvar=clinvar, **filters)


def motif_search(vcf_path, motif, vcf_splice=None, refseq='hg38', multiple=False, search_ref_threshold='estimate_threshold',
                 search_alt_threshold=-100, g_score=0.8, pseudocounts=0.8, background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
                 mtf_report=None, html_report=None, clinvar=None, **defined_filters):
    """search for single motif in vcf file"""
    # start reports
    if search_ref_threshold == 'estimate_threshold':
        search_ref_threshold = score.estimate_threshold(motif, g_score=g_score)
    if mtf_report is not None:
        mtf_report = start_report(mtf_report, motif, search_ref_threshold, clinvar=clinvar, **defined_filters)
    if html_report is not None:
        if not os.path.isfile(html_report):
            html_reports.start_HTML(html_report)
            html_reports.next_motif(motif, html_report, search_ref_threshold, clinvar=clinvar, first=True)
        else:
            html_reports.next_motif(motif, html_report, search_ref_threshold, clinvar=clinvar)
    # vcf readin
    if vcf_splice is not None:
        vcf_readin = vcf_splice(vcf_path, vcf_splice[0], vcf_splice[1])
    else:
        vcf_readin = read_in.vcf_handle(vcf_path=vcf_path)

    if multiple:
        vcf_readin = read_in.multipleSNP(vcf_readin, motif.length)
    if background == 'estimate_function':
        background = estimate_backgroud(refseq, motif)
    for var in vcf_readin:
        results = find_motif(record=var, motif=motif, refseq=refseq, search_ref_threshold=search_ref_threshold,
                             search_alt_threshold=search_alt_threshold, pseudocounts=pseudocounts,
                             background=background, **defined_filters)
        for res in results:
            result, diff = res
            if clinvar is not None: # if in vcf clinvar IDs
                clinvar = clinvar_url(record=var)
            if mtf_report:
                add_to_report(mtf_report, result, clinvar=clinvar)
            if html_report:
                html_reports.add_to_mtfHTML(result, html_report, clinvar=clinvar)
    if html_report is not None:
        html_reports.end_mtfHTML(html_report)
