import requests

def get_logo(matrix_id,path, timeout=10):
    matrix_url = "http://jaspar.genereg.net/api/v1/matrix/" + matrix_id
    response = requests.get(matrix_url, timeout=timeout).json()
    logo_url = response['sequence_logo']
    r = requests.get(logo_url, stream=True)
    if r.status_code == 200:
        with open(path, 'wb') as f:
            for chunk in r.iter_content(1024):
                f.write(chunk)
    return logo_url
