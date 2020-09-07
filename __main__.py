import argparse

from vcfmotifs.jasparapi.collection import generate_collection
from vcfmotifs.jasparapi.jsprMotif import idMotif
from vcfmotifs.jasparapi.species import species_collect
from vcfmotifs.vcfsearch.read_in import motif_handle
from vcfmotifs.vcfsearch.single_search import single_change, simpleSNP
from vcfmotifs.vcfsearch.vcf_search import motif_search, vcf_search


def main():
    parser = argparse.ArgumentParser(description="Specify parameters to run analyses.")

    subparsers = parser.add_subparsers(dest="command")
    parser.add_argument('-alt', type=float, required=False, help="search threshold for alternate sequence", default=-50)
    parser.add_argument('-ref', type=float, required=False, help="search threshold for referenced sequence", default='estimate_thresholds')
    parser.add_argument('-gscr', type=float,required=False, help="googness of fit value (0-1)")
    parser.add_argument('-pc', '--pseudocounts', required=False, help="pseudocounts for creating PWM matrix", default=0.8)
    parser.add_argument('-b','--background', required=False, help="background for creating PWM matrix",
                        default={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
    parser.add_argument('-dif', '--difference', type=float, required=False, help="minimal absolute difference between ref and alt score", default=5)
    parser.add_argument('-min', type=float, required=False, help="minimal search score is smaller then min threshold ")
    parser.add_argument('-max',type=float, required=False, help="maximal search score is bigger then max threshold ")
    parser.add_argument('-rdif',type=float, required=False, help="threshold for ref - alt scores")
    parser.add_argument('-adif', type=float, required=False, help="threshold for alt - ref scores")
    parser.add_argument('clinvar', type=bool, required=False, help="if clinvar info occurs in VCF", default=False)
    parser.add_argument('-refseq', required=True, help="reference sequence for vcf, for human genome type 'hg38' or 'hg37'")
    parser.add_argument('-multi', required=False, help="search for neighboring variants", default=True)

    # subparser for single variation search
    single = subparsers.add_parser('single', help='search for single variant')
    single.add_argument('-chrom', type=str, required=True, help="Chromosome")
    single.add_argument('-pos', type=int, required=True, help="Variation position")
    single.add_argument('-refn', type=str, required=True, help="REF nucleotide")
    single.add_argument('-altn', type=str, required=True, help="ALT nucleotide")
    # subparser for VCF search
    vcf = subparsers.add_parser('vcf', help='search in VCF file(s)')
    vcf.add_argument('-path', required=True, help="VCF file path" )
    vcf.add_argument('-start', help="start of searched region")
    vcf.add_argument('-end', help="end of searched region")
    # subparser for uploading motifs
    mtf = subparsers.add_parser('motif', help='motif JASPAR ID/collection/species')
    mtf.add_argument('-col',  help="type name of collection")
    mtf.add_argument('-spc', help="type name of the species")
    mtf.add_argument('-id', help="JASPAR ")
    mtf.add_argument('-mpath', help='path to directory with motif files')

    parser.add_argument('-html', required=False, help="output filename for html report")
    parser.add_argument('-csv', required=False, help="output filename for csv report")



    args = parser.parse_args()


    if args.command = 'jaspar':
        if args.id is not None:
            motif = idMotif(args.id)
        if args.col is not None:
            motif = generate_collection(args.col)
        if args.col is not None:
            motif = species_collect(args.spc)
        if args.mpath is not None:
            motif = motif_handle(args.path)

    if args.clinvar:
        cv = args.clinvar
        cv = cv.upper()
        if cv == "TRUE":
            args.clinvar = True

    defined_filters = {}

    if args.dif is not None:
        defined_filters['difference'] = args.dif
    if args.max is not None:
        defined_filters['max_score'] = args.max
    if args.min is not None:
        defined_filters['min_score'] = args.min
    if args.adif is not None:
        defined_filters['alt_diff'] = args.adif
    if args.rdfif is not None:
        defined_filters['ref_diff'] = args.rdif

    if args.command = 'single':
        snp = simpleSNP( CHROM=args.chrom, POS=args.pos, REF=args.refn , ALT=args.altn)
        single_change(snp, motif_readin=motif, resultsfile=args.csv, fastafile=args.refseq,
                      search_ref_threshold=args.ref, search_alt_threshold=args.alt, pseudocounts=args.pc,
                      background=args.b, **defined_filters)

    if args.command = 'vcf':
        if args.start is not None and args.end is not None:
            splice = (args.start, args.end)
        else: splice=None
        vcf_search(vcf_path=args.path, vcf_splice=splice,motifs_readin=motif, refseq=args.refseq,
                   multiple=args.multi, search_ref_threshold=args.ref, search_alt_threshold=args.alt, g_score=args.gscr,
                   pseudocounts=args.pc, background=args.b, mtf_report=args.csv, html_report=args.html,
                   clinvar=args.clinvar, **defined_filters)


if __name__ == '__main__':
    main()