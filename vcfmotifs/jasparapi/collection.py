"""
Collect motifs from collections
There are 10 collections avaliable
Collection names are: 'CORE', 'CNE', 'PHYLOFACTS', 'SPLICE', 'POLII', 'FAM','PBM','PBM_HOMEO','PBM_HLH', 'UNVALIDATED
"""
import requests

from vcfmotifs.jasparapi import check
from vcfmotifs.jasparapi.jsprMotif import urlMotif
from vcfmotifs.jasparapi.mfile import from_url


def jaspar_collection(coll_name, dir_path=None, mtf_format="jaspar", pagination=25, save=False, timeout=10):
    coll_name = coll_name.upper()  # all collections have uppercase name
    check.collection_check(coll_name)
    check.format_check(mtf_format)
    coll_name = 'collections/' + coll_name
    return get_collection(coll_name, dir_path, mtf_format, pagination, save, timeout)


def get_collection(coll_name, dir_path=None, mtf_format="jaspar", pagination=25, save=False, timeout=10):
    if save:
        return save_collection(coll_name, dir_path, mtf_format, pagination, timeout)
    else:
        return generate_collection(coll_name, pagination, timeout)


def save_collection(collection, dir_path=None, mtf_format="jaspar", pagination=25, timeout=10):
    next_url = "http://jaspar.genereg.net/api/v1/" + collection + '/?pagesize' + str(pagination)
    while next_url is not None:
        coll_page = requests.get(next_url, timeout=timeout).json()
        next_url = coll_page['next']
        records = coll_page['results']
        for record in records:
            mtx_url = record['url']
            from_url(mtx_url, dir_path=dir_path, mtf_format=mtf_format, user=False)


def generate_collection(collection, pagination=25, timeout=10):
    next_url = "http://jaspar.genereg.net/api/v1/" + collection + '/?pagesize' + str(pagination)
    while next_url is not None:
        coll_page = requests.get(next_url, timeout=timeout)
        if coll_page.status_code == 200:
            coll_page = coll_page.json()
            next_url = coll_page['next']
            records = coll_page['results']
            for record in records:
                mtx_url = record['url'] + '?format=jaspar'
                mtf = urlMotif(mtx_url)
                yield mtf
        else: raise requests.exceptions.HTTPError("HTTP response %s" %coll_page.status_code)


