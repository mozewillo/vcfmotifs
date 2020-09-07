
from Bio.Alphabet import generic_dna
from Bio.Seq import MutableSeq
from seqseek import Chromosome, BUILD38, BUILD37

from vcfmotifs.vcfsearch.read_in import read_fasta


def motif_scope(record, motif_len, refseq='hg38'):
    """return part of reference genome sequence of analysed scope"""
    if type(record) == list:
        locus, chrom = record[0].POS, record[0].CHROM
        diff = record[-1].POS - record[0].POS
        stop = locus + motif_len + diff
    else:
        locus, chrom = record.POS, record.CHROM
        stop = locus + motif_len
    start = locus - motif_len + 1
    if refseq == 'hg38':
        try:
            scope = Chromosome(chrom, assembly=BUILD38).sequence(start - 1, stop - 1)
        except ValueError:
            return None
    elif refseq == 'hg37':
        scope = Chromosome(chrom, assembly=BUILD37).sequence(start - 1, stop - 1)
    elif refseq is not None:
        scope = read_fasta(refseq, start, stop)
    else:
        raise ValueError("No reference sequence given. Choose 'hs38' or 'hg37' assembly, or upload fasta file.")
    return scope



def alt_scope(records, ref_scope, motif_len):
    """altered by VCF record sequence """
    ref_scope = str(ref_scope)
    altered_seq = MutableSeq(ref_scope, generic_dna)

    def alter(record, altered_seq, motif_len):
        if record.ALT[0] == None:
            return altered_seq
        # deletion
        if len(record.ALT) == 0:
            if len(record.REF) > 1:
                refs = list(record.REF)
                altered_seq[motif_len - 1] = str(record.ALT[0])
                for r in refs:
                    if motif_len < len(altered_seq):
                        del altered_seq[motif_len]
            else:
                del altered_seq[motif_len - 1]
        # substitution
        elif len(record.ALT[0]) == 1:
            altered_seq[motif_len - 1] = str(record.ALT[0])
        # insertion
        else:
            del altered_seq[motif_len - 1]
            alts = list(str(record.ALT[0]))
            alts.reverse()
            for n in alts:
                altered_seq.insert(motif_len, n)
        return altered_seq

    # multiple snps in motif range
    if type(records) == list:
        inital_len = len(ref_scope)
        inital=records[0].POS
        for i, record in enumerate(records):
            print(record)
            if i == 0:
                initial = record.POS
            diff = record.POS - initial
            shift = len(altered_seq) - inital_len
            altered_seq = alter(record, altered_seq, motif_len + diff)
    else:
        altered_seq = alter(records, altered_seq, motif_len)

    altered_seq = altered_seq.toseq()
    return altered_seq


def snps_position(records):
    """get all positions if more snps in motif scope"""
    if type(records) == list:
        # for snps in the same motif
        position = []
        for snps in records:
            position.append((snps.POS, snps.CHROM))
    else:
        position = (records.POS, records.CHROM)
    return position

