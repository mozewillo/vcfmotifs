"""
Collect motifs for one or more species
Jaspar has 63 species collections
"""

import requests

from vcfmotifs.jasparapi.check import taxon_check
from vcfmotifs.jasparapi.collection import get_collection


def tax_id(name, timeout=10):
    """Get tax_id from species name"""
    species = requests.get('http://jaspar.genereg.net/api/v1/species/?page=1&page_size=100', timeout=timeout)
    species = species.json()['results']
    for spc in species:
        if name == spc['species']:
            spc_id = spc['tax_id']
            return spc_id
    else:
        return "Species not found"


def get_ids(collection, pagination, timeout=10):
    """get only matrixes ids, no saving"""
    motifs_ids = []
    next_url = "http://jaspar.genereg.net/api/v1/collections/" + collection + '/?pagesize' + str(pagination)
    while next_url is not None:
        coll_page = requests.get(next_url, timeout=timeout).json()
        next_url = coll_page['next']
        records = coll_page['results']
        for record in records:
            record_id = record['matrix_id']
            motifs_ids.append(record_id)
    return motifs_ids


def species_collection(spc_id, dir_path=None, mtf_format="jaspar", pagination=25, save=False, timeout=10):
    """get species collection by species id"""
    collection = "species/" + spc_id
    return get_collection(collection, dir_path, mtf_format, pagination, save, timeout)


def species_collect(name, dir_path=None, mtf_format="jaspar", pagination=25, save=False, timeout=100):
    """get (and save) species collection by species name"""
    spc_id = tax_id(name, timeout)
    return species_collection(spc_id, dir_path, mtf_format, pagination, save, timeout)


def taxon_collect(tax_name, dir_path=None, mtf_format="jaspar", pagination=25, save=False, timeout=10):
    taxon_check(tax_name)
    collection = 'taxon/' + tax_name
    return get_collection(collection, dir_path, mtf_format, pagination, save, timeout)
