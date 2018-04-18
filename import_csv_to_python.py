#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 14:05:44 2018

@author: jake
"""
datadir = './Data/'
#infilename = 'CXA1820_spectrum.csv'
#outfilename = 'CXA1820_spectrum.json'

#infilename = 'linCIE2008v2e_1.csv'
#outfilename = 'luminosity_spectrum.json'

#infilename = 'epc660_sensitivity.csv'
#outfilename = 'epc660_sensitivity.json'

infilename = 'XB-D_Blue_spectrum.csv'
outfilename = 'XB-D_Blue_spectrum.json'

import csv
import json
import matplotlib.pyplot as plt

x = []
y = []
with open(datadir+infilename, newline='') as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
    for row in reader:
        x.append(row[0])
        y.append(row[1])
    csvfile.close()
        
#Round first and last elements
x[0] = round(x[0])
x[-1] = round(x[-1])
#Linear interpolate between end elements
beginjump = x[1]-x[0]
if beginjump > 1:
    delta = (y[1]-y[0])/beginjump
    while (x[1]-x[0]) > 1:
        x.insert(1,x[1]-1)
        y.insert(1,y[1]-delta)
        
endjump = x[-1]-x[-2]
if endjump >1:
    delta = (y[-1]-y[-2])/endjump
    while (x[-1]-x[-2]>1):
        x.insert(-1,x[-2]+1)
        y.insert(-1,y[-2]+delta)
        
for i in range(0,len(y)):
    y[i] = round(y[i],4)                
    

with open(datadir+outfilename,'w') as outfile:
    json.dump({'x':x,'y':y},outfile)
   
fig = plt.figure()    
plt.plot(x,y)