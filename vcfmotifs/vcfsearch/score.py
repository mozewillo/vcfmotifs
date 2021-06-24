from Bio import motifs, SeqIO
from seqseek import Chromosome, BUILD38, BUILD37

from vcfmotifs.vcfsearch.read_in import read_fasta


def estimate_backgroud(refseq, start=0, stop=False, chrom=None):
    """estimate the background frequency for nucleotides in referenced sequence"""
    if refseq == 'hg38':
        sequence = Chromosome(chrom, assembly=BUILD38).sequence(start, stop)
    elif refseq == 'hg37':
        sequence = Chromosome(chrom, assembly=BUILD37).sequence(start, stop)
    else:
        sequence = read_fasta(refseq, start=start, stop=stop)
    c = 0
    a = 0
    l = 0
    for i, n in enumerate(sequence):
        if i < start:
            continue
        l += 1
        if n == "C":
            c += 1
        if n == "A":
            a += 1
        if stop and i == stop:
            break
    a = a / l
    c = c / l
    background = {'A': a, 'C': c, 'G': c, 'T': a}
    return background


def pssm(motif, pseudocounts="calculate_pseudocounts", background={'A': 0.3, 'C': 0.2, 'G': 0.2, 'T': 0.3}):
    """return position-specific scoring matrix from motif counts"""
    if pseudocounts == "calculate_pseudocounts":
        pseudocounts = motifs.jaspar.calculate_pseudocounts(motif)
    pwm = motif.counts.normalize(pseudocounts=pseudocounts)
    pssm = pwm.log_odds(background=background)
    return pssm


def estimate_threshold(motif, g_score=0.8, pssmtx=None):
    """ set the threshold """
    if pssmtx is None:
        pssmtx = pssm(motif)
    max_scr = pssmtx.max
    thr = g_score * max_scr
    return thr


def ref_score(pssm, ref_scope, search_threshold=10):
    """ referenced sequence score """
    ref = list(pssm.search(ref_scope, search_threshold, both=False))
    return ref


def alt_score(altered_seq, motif_pssm, search_threshold=-10):
    alt = list(motif_pssm.search(altered_seq, search_threshold, both=False))
    return alt


def pair_scores(ref_scores, alt_scores):
    for r in ref_scores:
        for a in alt_scores:
            if r[0] == a[0]:
                pair = (r[0], (r[1], a[1]))
                yield pair


def gscore(ref, alt, motif):
    """Estimate non-length depended motif search value"""
    max_score = max(ref, alt)
    motif_pssm = pssm(motif)
    m_score = motif_pssm.search(motif.consesnus)
    G = max_score / m_score
    return G

