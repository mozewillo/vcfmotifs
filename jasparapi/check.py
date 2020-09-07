import requests

from vcfmotifs.jasparapi.exceptions import MatrixFormatError, JasparCollectionError


def format_check(mtf_format):
    formats = ["jaspar", "pfm", "transfac", "meme"]
    if mtf_format not in formats:
        exception = MatrixFormatError(mtf_format)
        raise exception


def collection_check(coll_name):
    # jaspar_coll = ["CNE", "CORE", "FAM", "PBM", "PBM_HLH", "PBM_HOMEO", "PHYLOFACTS", "POLII", "SPLICE",
    # "UNVALIDATED"]
    colls = requests.get("http://jaspar.genereg.net/api/v1/collections/")
    colls = colls.json()['results']
    for coll in colls:
        if coll_name == coll['name']:
            break
    else:
        exception = JasparCollectionError(coll_name)
        raise exception


def taxon_check(tax_name):
    pass
