#!/usr/bin/python
import numpy as np
import os

def get_ame_masses(path):
    masses = np.zeros((400,400))
    for line in open(path):
        if line[0]=='#':
            continue
        line = line.split()
        #print line
        if '-' in line[5]:
            continue
        if line[5][-1]=='#':
            line[5] = line[5][:-1]
        A = int(line[3])
        Z = int(line[2])
        mass = 0.931494028*1.0e-6*float(line[5])
        masses[A,Z]=mass
    return masses

def get_ame_q(A,Zi,Zf=None):
    if Zf == None:
        Zf = Zi-1
    if Zi == Zf:
        return 0.0
    m1 = mass_table[A,Zi]
    m2 = mass_table[A,Zf]
    if m1 == 0.0 or m2 == 0.0:
        return False
    else:
        #print m1*1000
         return (m1 - m2 - 0.000511)*1000

mass_table = get_ame_masses(os.path.dirname(__file__)+'/mass.mas12')
