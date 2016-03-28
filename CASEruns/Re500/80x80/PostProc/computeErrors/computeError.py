#!/usr/bin/env python
import matplotlib.pyplot as plt
import csv

# compute error by comparing the numerical solution of v-velocity to ghia's data.
ghiaFile = '../plot_Vvel_along_Horizontal_Line/ghia_Re100.csv'
dataFile = '../../v-velocity_in_x.csv'

xList1 = []
xList2 = []
vVelList1 = []
vVelList2 = []

# read Ghia's data
with open(ghiaFile, 'rb') as csvfile:
   reader = csv.DictReader(csvfile)
   for row in reader:
      xList1.append(float(row['x']))
      vVelList1.append(float(row['v-velocity']))


# read my simulation data
with open(dataFile, 'rb') as csvfile:
   reader = csv.DictReader(csvfile)
   for row in reader:
      xList2.append(float(row['x']))
      vVelList2.append(float(row['v-velocity']))

for x in xList1:
   for n in range(len(xList2)-1):
      
