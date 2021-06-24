
## Copyright © 2017 Saket Choudhary<saketkc__AT__gmail>
##https://github.com/saketkc/pyseqlogo/blob/master/pyseqlogo/pyseqlogo.py
## Modfied 2020 Aleksandra Mozwillo

import seaborn
import matplotlib.pyplot as plt

from vcfmotifs import vcfsearch

plt.style.use('seaborn-ticks')
from matplotlib import transforms
import matplotlib.patheffects
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
import matplotlib.font_manager

import numpy as np

COLOR_SCHEME = {'G': 'dimgray',
                'A': 'lightgray',
                'C': 'gray',
                'T': 'darkgrey'}

BASES = list(COLOR_SCHEME.keys())


class Scale(matplotlib.patheffects.RendererBase):

    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy) + affine
        renderer.draw_path(gc, tpath, affine, rgbFace)


def draw_logo(motif, ref, alt, position, pseudocounts=0.8, background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
              fontfamily='Arial', size=80, alt_color='gold', ref_color='red'):
    """Takes Motif object, draws motif_logo monochrome with color-marked SNP ref->alt"""
    if fontfamily == 'xkcd':
        plt.xkcd()
    else:
        rcParams['font.family'] = fontfamily
    SNP_COLORS = {alt: alt_color, ref: ref_color}
    counts = vcfsearch.score.pssm(motif, pseudocounts=pseudocounts, background=background)
    for n in counts:
        print(n)
        for i in range(motif.length):
            print(i)
            if counts[n][i] < 0:
                counts[n][i] = 0
            else:
                counts[n][i] = counts[n][i] * 0.65
            if i == 2:
                print('jej')
                if n == 'G':
                    counts[n][i] = counts[n][i] - 0.12
                if n == 'A':
                    counts[n][i] = counts[n][i] + 0.12
    scores_len = motif.length
    
    # set proper size in the end!
    fig, ax = plt.subplots(figsize=(scores_len, 2.5))
    # font
    font = FontProperties()
    font.set_size(size)
    font.set_weight('bold')
    # set axis
    ax.set_xticks(range(1, scores_len + 1))
    ax.set_yticks(range(0, 3))
    ax.set_xticklabels(range(1, scores_len + 1), rotation=90)
    ax.set_yticklabels(np.arange(0, 3, 1))
    seaborn.despine(ax=ax, trim=True)

    trans_offset = transforms.offset_copy(ax.transData,
                                          fig=fig,
                                          x=1,
                                          y=0,
                                          units='dots')

    for pos in range(scores_len):  # position,if i == snp.pos -> colored
        yshift = 0
        i = 0
        for base in counts:
            score = counts[base][pos]
            i+=1
            print(i, base, score)
            if pos == position and (base == alt or base == ref):
                txt = ax.text(pos + 1,
                              0,
                              base,
                              transform=trans_offset,
                              color=SNP_COLORS[base],
                              ha='center',
                              fontproperties=font)
            else:
                # if score[pos] < 0.03:
                #     continue
                txt = ax.text(pos + 1,
                              0,
                              base,
                              transform=trans_offset,
                              color=COLOR_SCHEME[base],
                              ha='center',
                              fontproperties=font)
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height * score
            trans_offset = transforms.offset_copy(txt._transform,
                                                  fig=fig,
                                                  y=yshift,
                                                  units='points')
        trans_offset = transforms.offset_copy(ax.transData,
                                              fig=fig,
                                              x=1,
                                              y=0,
                                              units='points')
    plt.show()
