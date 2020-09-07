import seaborn
from Bio import motifs
from matplotlib import transforms

"""https://omarwagih.github.io/ggseqlogo/ <- provides custom color schemes"""



import numpy as np
import pandas as pd
import seqlogo

pseudocounts = {"A":0.8, "C": 0.3, "G": 0.3, "T": 0.8}
random_ppm = np.random.dirichlet(np.ones(4), size=6)


#print(random_ppm)
# ppm = seqlogo.Ppm(random_ppm)
# pwn = seqlogo.ppm2pwm(ppm)
# # print(ppm)
#print(pwn)

with open ("/home/aleksandra/lic_py/test_data/motif/MA0506.1.jaspar") as handle:
    motif = motifs.read(handle, "jaspar")
    print("PWM", motif.counts)
    counts =(motif.counts)
    #pwn = np.ndarray(pwnn)
    motif = pd.DataFrame(motif.counts)
    print(motif)
    #ppm = seqlogo.pwm2ppm(pwn)
    pwm = seqlogo.Pfm(motif)
    # ppm = seqlogo.pfm2ppm(pwm)
    # ppm.alphabet = "DNA"
    # print(ppm)
    seqlogo.seqlogo(pwm, ic_scale=False, format='png',color_scheme='monochrome', size='medium', filename="/home/aleksandra/lic_py/lic/logo.png")

"""no seqlogo try --> easier to modify"""
