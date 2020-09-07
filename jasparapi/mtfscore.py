""" plot maximal threshold to length for all motifs for choosen motifs"""

import matplotlib.pyplot as plt

from vcfmotifs.vcfsearch import score

def max_thresholds_plot(motifs, save=None):
    """ plot maximal thresholds for all motifs in collection"""
    lengths = []
    refscr = []

    for mtf in motifs:
        lengths.append(mtf.length)
        pssm = score.pssm(mtf)
        scr = pssm.max
        refscr.append(scr)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set(title='Maximal PWM scores',
           xlabel='motif length',
           ylabel='max_score')
    ax.scatter(lengths, refscr, color='mediumblue')
    if save is not None:
        plt.savefig('max_scores.png')
    plt.show()
