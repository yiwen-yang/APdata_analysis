# l-DNB(set mean score as every sample's score)

import os
import scipy.stats as stat
import math
import numpy as np
import random
import time


begin=time.strftime("%Y-%m-%d:%H-%M-%S", time.localtime(time.time()))


normal = {}
f = open(r"/Pancreatitis/Huaxi_Pancreatitis_202206/raw data/reference_sample.csv")
flag = 0
for p in f:
    flag += 1
    if flag == 1:
        continue
    t = p.split(sep=",")
    normal[t[0]] = [float(t[i]) for i in range(1, len(t))]
f.close()

network = {}
keys = list(normal.keys())
n = len(keys)

# load the PPI network
background = []
f = open(r"/Pancreatitis/Huaxi_Pancreatitis_202206/raw data/string_interactions0.7AB.tsv", "r")
line = f.readlines()
for l in line:
    t = l.split()
    background.append(t)


f.close()


fw = open(r"/Pancreatitis/Huaxi_Pancreatitis/raw_data/reference_network_ppI_0.7AB.txt", "w")
pwd_n = 0
for bg in background:
    pwd_n += 1
    # print(pwd_n)
    if (bg[0] not in normal.keys()) | (bg[1] not in normal.keys()):
        continue
    r = stat.pearsonr(normal[bg[0]], normal[bg[1]])
    print(r)
    # The threshold of P-value need be set in here for Pearson Correlation Coefficient
    if r[1] < 0.05:
        fw.write(bg[0] + "\t" + bg[1] + "\t" + str(r[0]) + "\n")

fw.close()


print("Begin time is "+begin)
print("End time is "+time.strftime("%Y-%m-%d:%H-%M-%S",time.localtime(time.time())))

