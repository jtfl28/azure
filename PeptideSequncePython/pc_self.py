# -*- coding: utf-8 -*-
"""
Created on Mon May  7 16:16:42 2018

@author: jtfl2
"""
import numpy as np

masteraminos = {
    'A':{   '3letter':  'Ala',
	    'sc_mass':  15.0234,
	    'pk1':      2.35,
	    'pk2':      9.87,
	    'sc_hphob': 0.5},

    'R':{   '3letter':  'Arg',
	    'sc_mass':  100.0873,
	    'pk1':      1.82,
	    'pk2':      8.99,
	    'pk3':      12.48,
	    'sc_hphob': 1.81},

    'N':{   '3letter':  'Asn',
	    'sc_mass':  58.0292,
	    'pk1':      2.14,
	    'pk2':      8.72,
	    'sc_hphob': 0.85},

    'D':{   '3letter':  'Asp',
	    'sc_mass':  59.0132,
	    'pk1':      1.99,
	    'pk2':      9.9,
	    'pk3':      3.9,
	    'sc_hphob': 3.64},

    'C':{   '3letter':  'Cys',
	    'sc_mass':  46.9955,
	    'pk1':      1.92,
	    'pk2':      10.7,
	    'pk3':      8.3,
	    'sc_hphob': -0.02,
	    'extco':    125},

    'Q':{   '3letter':  'Gln',
	    'sc_mass':  72.0448,
	    'pk1':      2.17,
	    'pk2':      9.13,
	    'sc_hphob': 0.77},

    'E':{   '3letter':  'Glu',
	    'sc_mass':  73.0288,
	    'pk1':      2.1,
	    'pk2':      9.47,
	    'pk3':      4.07,
	    'sc_hphob': 3.63},

    'G':{   '3letter':  'Gly',
	    'sc_mass':  1.0078,
	    'pk1':      2.35,
	    'pk2':      9.78,
	    'sc_hphob': 1.15},

    'H':{   '3letter':  'His',
	    'sc_mass':  81.0452,
	    'pk1':      1.8,
	    'pk2':      9.33,
	    'pk3':      6.04,
	    'sc_hphob': 2.33},

    'I':{   '3letter':  'Ile',
	    'sc_mass':  57.0702,
	    'pk1':      2.32,
	    'pk2':      9.76,
	    'sc_hphob': -1.12},

    'L':{   '3letter':  'Leu',
	    'sc_mass':  57.0702,
	    'pk1':      2.33,
	    'pk2':      9.74,
	    'sc_hphob': -1.25},

    'K':{   '3letter':  'Lys',
	    'sc_mass':  72.0811,
	    'pk1':      2.16,
	    'pk2':      9.06,
	    'pk3':      10.54,
	    'sc_hphob': 2.8},

    'M':{   '3letter':  'Met',
	    'sc_mass':  75.0267,
	    'pk1':      2.13,
	    'pk2':      9.28,
	    'sc_hphob': -0.67},

    'F':{   '3letter':  'Phe',
	    'sc_mass':  91.0546,
	    'pk1':      2.2,
	    'pk2':      9.31,
	    'sc_hphob': -1.71},

    'P':{   '3letter':  'Pro',
	    'sc_mass':  41.039,
	    'pk1':      1.95,
	    'pk2':      10.64,
	    'sc_hphob': 0.14},

    'S':{   '3letter':  'Ser',
	    'sc_mass':  31.0183,
	    'pk1':      2.19,
	    'pk2':      9.21,
	    'sc_hphob': 0.46},

    'T':{   '3letter':  'Thr',
	    'sc_mass':  45.0339,
	    'pk1':      2.09,
	    'pk2':      9.1,
	    'sc_hphob': 0.25},

    'W':{	'3letter':  'Trp',
	    'sc_mass':  130.0655,
	    'pk1':      2.46,
	    'pk2':      9.41,
	    'sc_hphob': -2.09,
	    'extco':    5500},

    'Y':{   '3letter':  'Tyr',
	    'sc_mass':  107.0495,
	    'pk1':      2.2,
	    'pk2':      9.21,
	    'pk3':      10.07,
	    'sc_hphob': -0.71,
	    'extco':    1490},

    'V':{   '3letter':  'Val',
	    'sc_mass':  43.0546,
	    'pk1':      2.39,
	    'pk2':      9.74,
	    'sc_hphob': -0.46}}


def CharCount(theString,theChar): 
    result = 0
    for i in range(0, len(theString)):
        if theString[i] == theChar:
            result = result + 1
    return result


def count_aminos(sequence):
    rescounts = {}
    rescounts['A'] = CharCount(sequence, "A")
    rescounts['R'] = CharCount(sequence, "R")
    rescounts['N'] = CharCount(sequence, "N")
    rescounts['D'] = CharCount(sequence, "D")
    rescounts['C'] = CharCount(sequence, "C")
    rescounts['Q'] = CharCount(sequence, "Q")
    rescounts['E'] = CharCount(sequence, "E")
    rescounts['G'] = CharCount(sequence, "G")
    rescounts['H'] = CharCount(sequence, "H")
    rescounts['I'] = CharCount(sequence, "I")
    rescounts['L'] = CharCount(sequence, "L")
    rescounts['K'] = CharCount(sequence, "K")
    rescounts['M'] = CharCount(sequence, "M")
    rescounts['F'] = CharCount(sequence, "F")
    rescounts['P'] = CharCount(sequence, "P")
    rescounts['S'] = CharCount(sequence, "S")
    rescounts['T'] = CharCount(sequence, "T")
    rescounts['W'] = CharCount(sequence, "W")
    rescounts['Y'] = CharCount(sequence, "Y")
    rescounts['V'] = CharCount(sequence, "V")
    return rescounts

def calcmass(counts,seq):
    alphamass = 56.0136
    h2o_mass = 18.0105
    mass = alphamass*len(seq) + h2o_mass
    for key in counts:
        mass = mass + counts[key]*masteraminos[key]['sc_mass']
    mass = round(mass,4)
    
    return mass 

   
def cystine_count(cysteines):
    return (cysteines-(cysteines%2))/2

    
    
def calcec(counts,ranges):
    ec2 = counts['W']*masteraminos['W']['extco'] + counts['Y']*masteraminos['Y']['extco']
    ec1 = ec2 + cystine_count(counts['C'])*masteraminos['C']['extco']
        
    if not ranges[3]:
        ranges[3] = [ec1,ec1]
    elif ranges[3][0] > ec1:
        ranges[3][0] = ec1
    elif ranges[3][1] < ec1:
        ranges[3][1] = ec1
        
    if not ranges[4]:
        ranges[4] = [ec2,ec2]
    elif ranges[4][0] > ec2:
        ranges[4][0] = ec2
    elif ranges[4][1] < ec2:
        ranges[4][1] = ec2
        
    return [ec1, ec2,ranges]


def find_charge(a,b,pH):
    c = 0
    for key in a:
        if a[key]['count'] > 0:
            c = c + -a[key]['count']/(1 + pow(10, (a[key]['pk']-pH)))
	
    for key in b:
        if b[key]['count'] > 0:
            c += b[key]['count']/(1 + pow(10, (pH - b[key]['pk'])))
    c = round(c,3)
    return c

def calcpi(counts,sequence,ranges):
    #alert("pI being calculated for " + sequence)
    first_res = sequence[0]
    #alert("first residue is "+first_res)
    last_res = sequence[len(sequence)-1]
    #alert("last residue is "+last_res)
    acids = {   'C-term': {'count': 1,           'pk':masteraminos[first_res]['pk1']},
		'D':      {'count': counts['D'], 'pk':masteraminos['D']['pk3']},
		'E':      {'count': counts['E'], 'pk':masteraminos['E']['pk3']},
		'C':      {'count': counts['C'], 'pk':masteraminos['C']['pk3']},
		'Y':      {'count': counts['Y'], 'pk':masteraminos['Y']['pk3']}}
		
    bases = {   'N-term': {'count': 1,           'pk':masteraminos[last_res]['pk2']},
		'K':      {'count': counts['K'], 'pk':masteraminos['K']['pk3']},
		'R':      {'count': counts['R'], 'pk':masteraminos['R']['pk3']},
		'H':      {'count': counts['H'], 'pk':masteraminos['H']['pk3']}}    
    

    for pH in np.arange(0,14,0.01):
        if find_charge(acids,bases,pH) <= 0:
            break
        
    pI = round(pH,2)
    net_charge = round(find_charge(acids, bases, 7))
    if not ranges[0]:
        ranges[0] = [pI,pI]
    elif ranges[0][0] > pI:
        ranges[0][0] = pI
    elif ranges[0][1] < pI:
        ranges[0][1] = pI
        
    if not ranges[1]:
        ranges[1] = [net_charge,net_charge]
    elif ranges[1][0] > net_charge:
        ranges[1][0] = net_charge
    elif ranges[1][1] < net_charge:
        ranges[1][1] = net_charge
        
    return [pI,net_charge,ranges]

    
def calchphob(counts,ranges):
    hydrophobicity = 7.9
    for key in counts:
        hydrophobicity = hydrophobicity + counts[key] * masteraminos[key]['sc_hphob']
    hydrophobicity = round(hydrophobicity,2)  
    
    if not ranges[2]:
        ranges[2] = [hydrophobicity,hydrophobicity]
    elif ranges[2][0] > hydrophobicity:
        ranges[2][0] = hydrophobicity
    elif ranges[2][1] < hydrophobicity:
        ranges[2][1] = hydrophobicity
        
    return [round(hydrophobicity,2),ranges]

    
    
def calcProp(seq,ranges):
    P_count = count_aminos(seq)
    P_mass = calcmass(P_count,seq)
    [P_pI,P_charge,ranges] = calcpi(P_count,seq,ranges)
    [P_hphob,ranges] = calchphob(P_count,ranges)
    [P_ec1,P_ec2,ranges] = calcec(P_count,ranges)
    return [seq,float(P_mass),float(P_pI),P_charge,float(P_hphob),P_ec1,P_ec2,ranges]
    

def make_barcodes(properties):
    ranges = properties[7]
    binary = list(''.join(format(ord(x), 'b') for x in properties[0]))
    binary = [255*int(i) for i in binary]
    prop = [((float(properties[2])-ranges[0][0])/(ranges[0][1]-ranges[0][0]))*255,(((float(properties[3])-ranges[1][0])/(ranges[1][1]-ranges[1][0]))*255),(((float(properties[4])-ranges[2][0])/(ranges[2][1]-ranges[2][0]))*255),(((float(properties[5])-ranges[3][0])/(ranges[3][1]-ranges[3][0]))*255),(((float(properties[6])-ranges[4][0])/(ranges[4][1]-ranges[4][0]))*255)]
    return [binary,prop]

def make_matrix(aminos,n,current):
    matrix = np.zeros((n,n))
    column = 0
    if aminos[current] in 'IMVAL':
        order = list('IMVAL')
    elif aminos[current] in 'RHKDESTNQG':
        order = list('RHKDESTNQG')
    elif aminos[current] in 'YWF':
        order = list('IMVALRHKDESTNQGYWF')
    while current < len(aminos) and aminos[current] in order:
        for i in range(n):
            if aminos[current] == order[i]:
                matrix[i,column] = 1
                current = current + 1  
                column = column + 1
                break
    return [matrix, current]
        
        

def build_fingerprint(info, s):
    current = 0
    aminos = list(info[0])
    if s[0,0] == 1:
        n = 5
    elif s[1,0] == 1:
        n = 10
    elif s[2,0] == 1:
        n = 18
    [part1, current] = make_matrix(aminos,n,current)
        
    if s[0,1] == 1:
        n = 5
    elif s[1,1] == 1:
        n = 10
    elif s[2,1] == 1:
        n = 18
    [part2, current] = make_matrix(aminos,n,current)
    
    if s[0,2] == 1:
        n = 5
    elif s[1,2] == 1:
        n = 10
    elif s[2,2] == 1:
        n = 18
    [part3, current] = make_matrix(aminos,n,current)
    
    return [part1,part2,part3]


class branch:
    def __init__(self,s):
        self.domain = s
        empty = [[]]
        self.parts = empty*len(s)
        
    def take_from_old(self,s,oldpart):
        self.domain = s
        self.parts = oldpart

    def add_to_branch(self,new,part):
        i = 0
        if len(self.parts[part]) == 0:
            self.parts[part] = [new]
            return i
        else:
            for prints in self.parts[part]:
                if (prints == new).all():
                    return i
                    break
                else: 
                    i = i + 1
            self.parts[part].append(new)
            return i
        
    def mergebranch(self,new):
        for i in range(3):
            for newparts in new.parts[i]:
                isin = 0
                for oldparts in self.parts[i]:
                    if (newparts == oldparts).all():
                        isin = 1
                if isin == 0:
                    self.parts[i].append(newparts)
                    
    def getindex(self, fingerprint):
        index = [-1,-1,-1]
        for i in range(3):
            j = 0
            for part in self.parts[i]:
                if (fingerprint[i] == part).all():
                    index[i] = j
                    break
                j = j + 1
        return index


class tree:
    def __init__(self):
        self.branches = []
        
    def addtotree(self, newbranch):
        if len(self.branches) == 0:
            self.branches.append(newbranch)
        else:
            i = 0
            for bran in self.branches:
                if (bran.domain == newbranch.domain).all():
                    self.branches[i].mergebranch(newbranch)
                    return
                else:
                    i = i + 1
            self.branches.append(newbranch)
            return
            
    def getbranchindex(self, branch):
        i = 0
        for b in self.branches:
            if (branch.domain == b.domain).all():
                return i
            i = i + 1

                
        
 
    
    