import os
import zipfile
import vcf
import re
import pandas as pd
from Bio import motifs
from Bio import SeqIO



def vcf_handle(vcf_path):
    """search for a single vcf file or directory with multiple vcf files"""
    def get_records(vcf_path):
        with open(vcf_path, 'r') as v:
            yield from vcf.Reader(v)
    # gzip
    if vcf_path.endswith('.vcf.gz'):
        f = open(vcf_path, 'rb')
        # with open(vcf_path, 'rb') as f:
        vcf_read = vcf.Reader(f, compressed=True)
        try:
            for v in vcf_read:
                yield v
        except StopIteration:
            f.close()
        except IndexError:
            f.close()
            return None
    # zip
    if zipfile.is_zipfile(vcf_path):
        with zipfile.ZipFile(vcf_path) as f:
            for name in f.namelist():
                if name.endswith('.vcf'):
                    get_records(vcf_path)
                # change vcf
    # directory
    elif os.path.isdir(vcf_path):
        files = os.listdir(vcf_path)
        vcfs = [file for file in files if file.endswith('.vcf')]
        for file in vcfs:
            single_path = vcf_path + file
            get_records(single_path)
    # vcf file
    elif os.path.isfile(vcf_path) and vcf_path.endswith('.vcf'):
        print("read")
        with open(vcf_path, 'r') as v:
            yield from vcf.Reader(v)
    else:
        raise ValueError("Given path %s does not contain vcf file" % vcf_path)


def motif_handle(motif_path):
    """return Motif object from given path"""
    # zip
    if zipfile.is_zipfile(motif_path):
        with zipfile.ZipFile(motif_path) as f:
            for name in f.namelist():
                if name.endswith(('.jaspar', '.transfac', '.meme', 'pfm')):
                    fext = name.split(".")[-1]
                    with f.open(name) as mtf:
                        motif = motifs.read(mtf, fext)
                        yield motif
    # motif files
    elif os.path.isfile(motif_path) and motif_path.endswith(('.jaspar', '.transfac', '.meme', 'pfm')):
        fext = motif_path.split(".")[-1]
        with open(motif_path) as handle:
            motif = motifs.read(handle, fext)
            return motif

    # directory
    elif os.path.isdir(motif_path):
        files = os.listdir(motif_path)
        mtfs = [file for file in files if file.endswith(('.jaspar', '.transfac', '.meme', 'pfm'))]
        for file in mtfs:
            single_path = motif_path + '/' + file
            with open(single_path) as handle:
                fname, fext = os.path.splitext(single_path)
                fext = fext[1:]  #
                motif = motifs.read(handle, fext)
                yield motif


def multipleSNP(recordReader, motif_len):
    """more than one SNP in scope of one motif"""
    new_scope = recordReader.__next__()
    chrom = new_scope.CHROM
    in_scope = [new_scope]
    new_scope = None
    for record in recordReader:
        if new_scope is not None:
            chrom = new_scope.CHROM
            in_scope = [new_scope]
            new_scope = None
        elif record.CHROM == chrom and record.POS - in_scope[-1].POS < motif_len - 1:
            in_scope.append(record)
        else:
            new_scope = record
            yield in_scope
    # last yield
    if new_scope is not None:
        ns = [new_scope]
        yield ns
    else:
        yield in_scope


def read_fasta(fasta_file, start, stop ):
    # readin reference sequence from fasta file
    with open(fasta_file, 'r') as f:
        sequence = SeqIO.read(f, "fasta")
        return sequence.seq[start - 1: stop - 1]


def name_csv(filename):
    if filename is not None and not filename.endswith('.csv'):
        filename = filename + '.csv'
    return filename


def vcf_splice(vcffile, start, stop):
    with open (vcffile, 'r') as vcfs:
        records = vcf.Reader(vcfs)
        for record in records:
            if start <= record.POS <= stop:
                yield record


def clinvar_url(record):
    if type(record) == list:
        info = []
        urls = []
        for r in record:
            url = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/' + r.ID
            urls.append(url)
            if 'CLNSIG' in r.INFO.keys():
                info.append(r.INFO['CLNSIG'])
            else:
                info = r.ID
                break
        htmllink = ""
        csvlink = ""
        for i in range(len(info)):
            htmllink  += '<a href=\"'+ str(urls[i]) + '">' + str(info[i]) + '</a>'
            csvlink = '=HYPERLINK("' + str(urls[i]) + '";"' + str(info[i]) + '")'
    else:
        if 'CLNSIG' in record.INFO.keys():
            info = record.INFO['CLNSIG']
        else:
            info = "ClinVar : " + record.ID
        var_url = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/' + record.ID
        htmllink  = '<a href=\"'+ str(var_url) + '">' + str(info) + '</a>'
        csvlink = '=HYPERLINK("' + str(var_url) + '";"' + str(info) + '")'
    return csvlink, htmllink


def process_sites(sitesfile):
    """ process sites file from JASPAR sites fasta file"""
    with open('/home/aleksandra/Desktop/lic_py/human/ETS/site/MA0080.5.sites', 'r') as sites:
        for site in sites:
            if site.startswith(">"):
                site = site[site.rfind('chr'):]
                site = re.findall('\d+', site)
                chrom = site[0]
                start = site[1]
                end = site[2]
                yield chrom, start, end

def dataframe_csv(mtf_csvfile):
    """readin as pandas DataFrame in batches filter by columns operations on alt and ref seq"""
    mtfdata = pd.read_csv(mtf_csvfile, skiprows=3)
    if len(mtfdata.columns) == 9:
        mtfdata.columns = ['ID', 'chrom', 'positions', 'polymorphism', 'sequences', 'ref_score', 'alt_score', 'motif_loc', 'clinvar']
    else:
        mtfdata.columns = ['ID', 'chrom', 'positions', 'polymorphism', 'sequences', 'ref_score', 'alt_score',
                           'motif_loc']
    return mtfdata
