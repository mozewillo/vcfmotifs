import os
import requests
from Bio import motifs

from vcfmotifs.jasparapi.check import format_check


def save_motif(motif_id, dir_path, mtf_format="jaspar", user=True):
    if user:
        format_check(mtf_format)
    motif_url = "http://jaspar.genereg.net/api/v1/matrix/" + motif_id  # + "/?format=" + mtf_format
    from_url(motif_url, dir_path=dir_path, mtf_format=mtf_format, user=False)


def from_url(motif_url, dir_path, mtf_format="jaspar", filename=None, user=True, timeout=10):
    """Save motif in given format from motif_url"""
    if user:
        format_check(mtf_format)
    mtf = requests.get(motif_url, timeout=timeout)
    # file name
    if filename is None:
        mtf_name = mtf.json()
        mtf_name = mtf_name['name']
    file_name = mtf_name + "." + mtf_format
    if not dir_path.endswith('/'):
        dir_path = dir_path + '/'
    mtf_path = dir_path + file_name

    # save motif in given format
    record = requests.get(motif_url + "/?format=" + mtf_format, timeout)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    with open(mtf_path, 'w+') as f:
        f.write(record.text)
        print("Motif %s saved at %s" % (mtf_name, mtf_path))


def get_Motif(mtf_path, mtf_format):
    """Get Bio.motifs.Motif object from file"""
    format_check(mtf_format)
    with open(mtf_path, 'r') as f:
        mtf = motifs.read(f, mtf_format)
    return mtf
