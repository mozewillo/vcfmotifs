"""CREATE HTML REPORTS"""


def add_to_mtfHTML(result, fileout, r_type="mtf", clinvar=None, motif=None):
    """ add record to html """
    if r_type == "mtf":
        info = [str(result["variant_ID"]), str(result["chrom"]), str(result["var_pos"]),  str(result["polymorphism"]), str(result["sequences"]), str(result["ref_score"]),
                str(result["alt_score"]),
                str(result["motif_loc"])]
    elif r_type == "snp":
        info = [str(motif.matrix_id), str(result["sequences"]), str(result["ref_score"]), str(result["alt_score"]),
                str(result["motif_loc"])]
    report = "<tr>"
    for i in info:
        report += "<td>" + i + "</td>\n"
    if clinvar is not None:
        report += "<td>" + clinvar[1] + "</td>\n"
    report += "</tr>\n"
    with open(fileout, 'a') as out:
        out.writelines(report)


def start_HTML(fileout):
    """star html"""
    report = '<!DOCTYPE html> \n<html lang="en">\n' + '<head>\n' + '<meta charset="utf-8">' + '<link rel="stylesheet" ' + \
             'href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" ' + \
             'integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" ' \
             'crossorigin="anonymous">\n' + '</head>' + '<body style="font-family:' + "'Monospace'" + ';size:8px;color:#3e3e5b;">'
    with open(fileout, 'a') as out:
        out.writelines(report)


def start_single(record, fileout, **filters):
    """Start html report for sinle change"""
    start_HTML(fileout)
    report = ''
    if type(record) == list:
        positions = []
        ref = []
        alt = []
        varID = []
        for r in record:
            positions.append(str(r.POS))
            ref.append(r.REF)
            alt.append(str(r.ALT))
            varID.append(r.ID)
        header = "<h1>" + " Neighboring variation analysis VarID: " + str(record.ID) + "</h1>"
        info = ["polymorphism : " + str(record.ID), "record position : " + str(record.POS),
                "\tREF :" + str(record.ALT)
                + str(record.ALT) + "applied filters: " + str(filters)]
    else:
        header = "<h1>" + " Single variation analysis VarID: " + str(record.ID) + "</h1>"
        info = ["record position : " + str(record.POS), "\tREF :" + str(record.ALT)
                + str(record.ALT) + "applied filters: " + str(filters)]
    columns = ['matrix_id', 'sequences', 'ref_score', 'alt_score', 'motif_loc']
    report += header
    for i in info:
        report += "<p>" + i + "</p>\n"
    report += '<table class="table table-striped">\n<thead class="thead-dark">\n'
    # add labels
    report += "<tr>\n"
    for l in columns:
        report += "<th>" + l + "</th>\n"
    report += "</tr>\n"
    with open(fileout, 'a') as out:
        out.writelines(report)


def next_motif(motif, fileout, ref_threshold, clinvar=None, first=False):
    "start table for next motif"
    report = ""
    if not first:
        report = '</table>\n' + "</div>\n"
    info = ["motif_id : " + motif.matrix_id, " motif_consensus_sequence : " + str(motif.consensus),
            "search_ref_threshold : " + str(ref_threshold)]
    labels = ['ID', 'chrom', 'position', 'polymorphism', 'altered_seq', 'ref_score', 'alt_score', 'motif_loc']
    report += "<h1>" + motif.name + "</h1>\n"
    for i in info:
        report += "<p>" + i + "</p>\n"
    report += '<div class="table-responsive">\n' + '<table class="table table-striped">\n<thead class="thead-dark">\n'
    # add labels
    report += "<tr>\n"
    for l in labels:
        report += "<th>" + l + "</th>\n"
    if clinvar is not None:
        report += "<th>" + "ClinVar" + "</th>\n"
    report += "</tr>\n"
    with open(fileout, 'a') as out:
        out.writelines(report)


def end_mtfHTML(fileout):
    """close html tags"""
    closing = "</table>\n</div>\n</html>"
    with open(fileout, 'a') as out:
        out.writelines(closing)
