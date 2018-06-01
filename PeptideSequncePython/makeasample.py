# -*- coding: utf-8 -*-
"""
Created on Thu May 31 10:48:08 2018

@author: jtfl2
"""

import itertools as it
import pc_self as pcs
import numpy as np
import pickle
import progressbar as pb
import pandas as pd
import time

aromatic = {}
hydrophobic = {}
hydrophilic = {}
per = {}
ranges = [[],[],[],[],[]]
structure =  np.zeros((3,3))
branch = pcs.branch(structure)
tree = pcs.tree()
tree.addtotree(branch)

hydrophobic['3'] = it.product('IMVAL', repeat=3)
hydrophilic['5'] = it.product('HRKDETSGNQ', repeat=5)
aromatic['4'] = it.product('YWF','IMVALHRKDETSGNQYWF','IMVALHRKDETSGNQYWF','IMVALHRKDETSGNQYWF')

structure[0,0] = structure[1,1] = structure[2,2] = 1

phob = []
for y in hydrophobic['3']:
    phob.append("".join(y[0]+y[1]+y[2]))           
r = pd.DataFrame(phob)
rand= r.sample(frac=0.1,replace=False)
phobrand = []
for i in rand[0]:
    phobrand.append(i)

phil = []
for y in hydrophilic['5']:
    phil.append("".join(y[0]+y[1]+y[2]+y[3]+y[4]))
r = pd.DataFrame(phil)
rand= r.sample(frac=0.1,replace=False)
philrand = []
for i in rand[0]:
    philrand.append(i)
        
aro = []
for y in aromatic['4']:
    aro.append("".join(y[0]+y[1]+y[2]+y[3]))
r = pd.DataFrame(aro)
rand= r.sample(frac=0.1,replace=False)
arorand = []
for i in rand[0]:
    arorand.append(i)
        
sample = {}
sample['12'] = it.product(phobrand,philrand,arorand)


f = open('sampleof12lengthpeptides.csv', 'w')
f.write("Tree Branch, Branch Part 1, Branch Part 2, Branch Part 3,Mass,Isoelectric Point (pI),Net Charge,Hydrophobicity,Extinction coefficient 1,Extinction coefficient2" + "\n")
p = ""
total = len(phobrand)*len(philrand)*len(arorand)
i = 0
pb.printProgressBar(i, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
for p in sample['12']:
    seq = "".join(p[0] + p[1] + p[2])
    Info = pcs.calcProp(seq,ranges)
    ranges = Info[7]
    fingerprint = pcs.build_fingerprint(Info, structure)
    [first,second,third] = [branch.add_to_branch(fingerprint[0], 0), branch.add_to_branch(fingerprint[1],1), branch.add_to_branch(fingerprint[2],2)]
    f.write(str(tree.getbranchindex(branch)) + ',' + str(first) + ',' + str(second) + ',' + str(third) + ',' + str(Info[1])+ ','+ str(Info[2])+ ','+ str(Info[3])+ ','+ str(Info[4])+ ','+ str(Info[5])+ ','+ str(Info[6])   + "\n")
    i = i + 1
    pb.printProgressBar(i, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
f.close()

tree.addtotree(branch)

with open("peptidetree.txt", "wb") as fp:
   pickle.dump(tree, fp)