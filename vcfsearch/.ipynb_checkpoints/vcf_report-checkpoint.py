""" frequency for all vcf return "modified vcf" with changes score information about motif"""
""" report for all vcf  --frequency through sequence --most frequent motif changes
detected changes -- same motif changes(different SNP same motif) """
import csv
import pandas as pd

### stworz tabele
# dla pojedynczego snp
# jakie motywy zaburza
# jak czesto zaburza
"""DataFrame for separate snp"""

def search_SNP(snp_position, snp_ref, snp_alt):


"""DataFrame for whole vcf"""
def create_note(record):
    """create dataframe from first record analyses(stored in dict)"""
    df = pd.DataFrame(data=record)


def note_record(record_data, record):
    """ -altered<->ref -motif columns"""
    record_data.append(record, ignore_index=True)


def write_record(csv_file, record_data):
    """use dataFrame to store results then put them into csv"""
    with open(csv_file, 'w', newLine='') as f:
        f_writer = csv.writer(f)
