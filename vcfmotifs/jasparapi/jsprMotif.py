import io
import requests
from Bio import motifs

"""Create Bio.Motifs using matrix_id or matrix_url"""


def idMotif(matrix_id, user=True, timeout=10):
    """get Motif object by matrix id"""
    matrix_url = "http://jaspar.genereg.net/api/v1/matrix/" + matrix_id + "/?format=jaspar"
    response = requests.get(matrix_url, timeout=timeout)
    mf_data = io.StringIO()
    mf_data.write(response.text)
    mf_data.seek(0)
    motif = motifs.read(mf_data, "jaspar")
    mf_data.close()
    return motif


def urlMotif(motif_url, timeout=10):
    """get Motif object by motif_url"""
    response = requests.get(motif_url, timeout=timeout)
    mf_data = io.StringIO()
    mf_data.write(response.text)
    mf_data.seek(0)
    motif = motifs.read(mf_data, "jaspar")
    mf_data.close()
    return motif


def genMotifs(motif_ids):
    """generate Motifs from list of matrix ids"""
    for id in motif_ids:
        motif = idMotif(matrix_id=id)
        yield motif
