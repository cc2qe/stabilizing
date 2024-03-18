#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from itertools import chain
import pandas as pd
from sys import argv

pos_l = []
S_l = []
genotypes_l = []
freq_l = []

# "low": minor allele frequency 0.005 - 0.05
def is_low(x):
    return ((x >= 0.005 and x < 0.05) or (x <= 0.995 and x > 0.95))

# "common": minor allele frequency >= 0.05
def is_common(x):
    return (x >= 0.05 and x <= 0.95)

with open(argv[1]) as f:
    for line in f:
        if '#' in line:
            continue
        l = line.split('\t')
        info = l[7]
        pos = int(l[1])
        S = float(info.split(";")[1].split('=')[1])
        genos = l[9:]
        genotypes = [[int(j) for j in i.split('|')] for i in genos]
        genotypes_flattened = list(chain.from_iterable(genotypes))
        pos_l.append(pos)
        S_l.append(S)
        freq_l.append(np.sum(genotypes_flattened)/len(genotypes_flattened))
        genotypes_l.append(genotypes_flattened)

S_mean = np.mean(S_l)
S_sd = np.std(S_l)

S_norm = [(S_i - S_mean)/S_sd for S_i in S_l]

dist_l = []
prod_l = []
class_l = []
corr_l = []
corr_class_l = []
cov_l = []
for idx in range(0, len(pos_l)):
    # print(idx)
    for jdx in range(idx+1, len(pos_l)):
        if freq_l[jdx] == 1.0 or freq_l[idx] == 1.0:
            continue
        dist = np.abs(pos_l[jdx] - pos_l[idx])
        if dist > 100000:
            continue
        Bi = S_norm[idx]
        Bj = S_norm[jdx]
        corr = np.corrcoef(genotypes_l[jdx], genotypes_l[idx])[0][1]
        cov = np.cov(genotypes_l[jdx], genotypes_l[idx])[0][1]
        prod_l.append(Bi*Bj)
        dist_l.append(dist)
        corr_l.append(corr)
        cov_l.append(cov)

        if is_low(freq_l[idx]) and is_low(freq_l[jdx]):
            class_l.append("low-low")
        elif is_common(freq_l[idx]) and is_common(freq_l[jdx]):
            class_l.append("common-common")
        elif (is_common(freq_l[idx]) and is_low(freq_l[jdx])) or (is_low(freq_l[idx]) and is_common(freq_l[jdx])):
            class_l.append("common-low")
        else:
            class_l.append("NA")
            
        if corr < 0:
            corr_class_l.append('neg')
        elif corr > 0:
            corr_class_l.append('pos')
        else:
            corr_class_l.append('NA')


df = pd.DataFrame({
    'dist':dist_l,
    'prod':prod_l,
    'corr':corr_l,
    'class':class_l,
    'corr_class':corr_class_l,
    'cov':cov_l
})

bins = [0, 100, 1000, 10000, 100000]
labels = [100, 1000, 10000, 100000]
df['bin'] = pd.cut(x = df['dist'], bins = bins, labels = labels, include_lowest = True)
df.groupby(['bin']).mean().to_csv(argv[2] + ".covar.csv")
df.groupby(['bin', 'class', 'corr_class']).mean().to_csv(argv[2] + ".bin_covar.csv")
