# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 13:38:09 2018

@author: jtfl2
"""
import itertools as it
import pc_self as pcs
import numpy as np
import pickle


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
per['12'] = it.product(hydrophobic['3'],hydrophilic['5'],aromatic['4'], repeat = 1)

f = open('12lengthpeptides.csv', 'w')
f.write("Tree Branch, Branch Part 1, Branch Part 2, Branch Part 3,Mass,Isoelectric Point (pI),Net Charge,Hydrophobicity,Extinction coefficient 1,Extinction coefficient2" + "\n")
p = ""

i = 0
for p in per['12']:
    seq = "".join(p[0] + p[1] + p[2])
    Info = pcs.calcProp(seq,ranges)
    ranges = Info[7]
    fingerprint = pcs.build_fingerprint(Info, structure)
    [first,second,third] = [branch.add_to_branch(fingerprint[0], 0), branch.add_to_branch(fingerprint[1],1), branch.add_to_branch(fingerprint[2],2)]
    f.write(str(tree.getbranchindex(branch)) + ',' + str(first) + ',' + str(second) + ',' + str(third) + ',' + str(Info[1])+ ','+ str(Info[2])+ ','+ str(Info[3])+ ','+ str(Info[4])+ ','+ str(Info[5])+ ','+ str(Info[6])   + "\n")
    print(i)
    i = i + 1
f.close()

tree.addtotree(branch)

with open("peptidetree.txt", "wb") as fp:
   pickle.dump(tree, fp)
