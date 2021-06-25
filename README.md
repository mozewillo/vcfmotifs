![alt text](https://github.com/mozewillo/vcfmotifs/blob/master/additional_files/logo.png) \
*Tool for Detecting DNA polymorphisms influencing transcription factorbinding*

**Background** \
Human genome sequence contains millions of DNA polymorphisms, that don’t affect pheno-type in an obvious way. Most of them are located in the non-coding part of the genome,making it more difficult to understand their function. The aim of this work is to develop atool to detect the subset of polymorphisms most likely to disrupt transcription. This tool isimplemented in Python programming language and can be used by programmers as a library

**Method** \
Vcfmotifs is searching through a sequence by matching motif position weight matrix that can be obtained by conversions (mentioned in the Data section) from the position frequencymatrix(PFM). 
We can use motif data from different databases mostly they provide at least one of the given formats: ’.clusterbuster’, ’.pfm’, ’.jaspar’, ’.transfac’ -then we use the Bio.motifspackage to read it in. 
We also included in Vcfmotifs subpackage called *jasparapi* that givesus access to JASPAR RESTful API and giving us ready generators yielding Bio.Motif objects.
So we can easily search for motifs by Jaspar IDs or whole collections.

The score is counted based on log-odds calculated from the frequencies of nucleotides observed on certain position in the motif.
The Position Frequency Matrix is corrected by the pseudocounts to avoid probabilties equal zero (pseudocounts can be given by the user).
To get the proper odds user can set the background frequencies of the nucleotides elsewise the backgroung probabilities will be estimated from the sequence given.
The scores of the motif is calculated for the consensus motif sequence (*ref_score*) and the one altered by the polymorphism (*alt_score*).
Main obstacle in analysis using this tool is setting the threshold for the change that can be significant in motif functionality.

**Parameters of search** \
The proposed search requires thresholds for scores of interest of the user. 
There are few scores that the filtering can be applied on, multiple filters can be applied on one search.\
*Scores for each sequence (thresholds can be set):*
<ul>
<li>ref_score - score that is obtained by the consensus sequence of the motif </li>
<li> alt_score - score that is obtained by the sequence containing the polymorphism </li>
</ul>
List of filters:
<ul>
<li>abs_difference - the absolute difference between the scores of motif searches for referenced and altered sequences </li>
<li>max_score - specify threshold for motif match that one of two scores needs to exceed </li>
<li> min_score - defining how small score need to be to consider as disturbing linkage change </li>
<li> altmotif_differ - how relevant change may be ref -> alt, if motif occurs in altered seq </li>
<li> refmotif_differ - how relevant change may be ref -> alt, if motif occurs in referenced seq </li>
</ul>

The tool proposes the adjusting the threshold for each motif, since the scores may vary for long vs. short motifs, and those more or less conserved. 
The adjustment is based on the upper threshold proportion of change (alt_score divided by ref_score) that user finds interesing for reporting.

```
usage: pars.py [-h] [-alt ALT] [-ref REF] [-gscr GSCR] [-pc PSEUDOCOUNTS] [-b BACKGROUND] [-diff DIFFERENCE]
               [-min MIN] [-max MAX] [-rdif RDIF] [-adif ADIF] [-clinvar CLINVAR] -refseq REFSEQ [-multi MULTI]
               [-html HTML] [-csv CSV]
               {single,vcf,motif} ...

Specify parameters to run analyses.

positional arguments:
  {single,vcf,motif}
    single              search for single variant
    vcf                 search in VCF file(s)
    motif               motif JASPAR ID/collection/species

optional arguments:
  -h, --help            show this help message and exit
  -alt ALT              search threshold for alternate sequence
  -ref REF              search threshold for referenced sequence
  -gscr GSCR            googness of fit value (0-1)
  -pc PSEUDOCOUNTS, --pseudocounts PSEUDOCOUNTS
                        pseudocounts for creating PWM matrix
  -b BACKGROUND, --background BACKGROUND
                        background for creating PWM matrix
  -diff DIFFERENCE, --difference DIFFERENCE
                        minimal absolute difference between ref and alt score
  -min MIN              minimal search score is smaller then min threshold
  -max MAX              maximal search score is bigger then max threshold
  -rdif RDIF            threshold for ref - alt scores
  -adif ADIF            threshold for alt - ref scores
  -clinvar CLINVAR      if clinvar info occurs in VCF
  -refseq REFSEQ        reference sequence for vcf, for human genome type 'hg38' or 'hg37'
  -multi MULTI          search for neighboring variants
  -html HTML            output filename for html report
  -csv CSV              output filename for csv report

```



**Example of html report** \
Example of result of search for transcription factor ETV6 binding side variations. (This factor plays a role in hematopoiesis and malignant transformation. See motif: [jaspar motif](http://jaspar.genereg.net/matrix/MA0645.1/)
![alt text](https://github.com/mozewillo/vcfmotifs/blob/master/additional_files/example_html.png)
