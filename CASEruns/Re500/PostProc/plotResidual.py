#!/usr/bin/env python
import matplotlib.pyplot as plt
import csv

dataFile1 = '../20x20/u-residualLog.csv'
dataFile2 = '../40x40/u-residualLog.csv'
dataFile3 = '../80x80/u-residualLog.csv'

iterList1 = []
iterList2 = []
iterList3 = []
uResList1 = []
uResList2 = []
uResList3 = []

# read 20x20 case
with open(dataFile1, 'rb') as csvfile:
   reader = csv.DictReader(csvfile)
   for row in reader:
      iterList1.append(float(row['nIter']))
      uResList1.append(float(row['residual']))

# read 40x40 case
with open(dataFile2, 'rb') as csvfile:
   reader = csv.DictReader(csvfile)
   for row in reader:
      iterList2.append(float(row['nIter']))
      uResList2.append(float(row['residual']))

# read 80x80 case
with open(dataFile3, 'rb') as csvfile:
   reader = csv.DictReader(csvfile)
   for row in reader:
      iterList3.append(float(row['nIter']))
      uResList3.append(float(row['residual']))

iterList1 = [i/1000 for i in iterList1]
iterList2 = [i/1000 for i in iterList2]
iterList3 = [i/1000 for i in iterList3]

MinX = 0
MaxX = max(iterList3)
MinY = 0.00005
MaxY = 1.0

p = plt.plot(iterList1,uResList1, 'k-', label="20x20")
plt.setp(p, linewidth='3.0')
p = plt.plot(iterList2,uResList2, 'r-', label="40x40")
plt.setp(p, linewidth='3.0')
p = plt.plot(iterList3,uResList3, 'b-', label="80x80")
plt.setp(p, linewidth='3.0')

plt.axis([MinX,MaxX, MinY, MaxY])
plt.xscale('linear')
plt.yscale('log')
plt.xlabel('# iterations [x1000]', fontsize=22)
plt.ylabel('u-residual', fontsize=22)

plt.grid(True)
ax = plt.gca()
xlabels = plt.getp(ax, 'xticklabels')
ylabels = plt.getp(ax, 'yticklabels')
plt.setp(xlabels, fontsize=18)
plt.setp(ylabels, fontsize=18)
plt.legend(
          loc='upper right',
          borderpad=0.25,
          handletextpad=0.25,
          borderaxespad=0.25,
          labelspacing=0.0,
          handlelength=2.0,
          numpoints=1)
legendText = plt.gca().get_legend().get_texts()
plt.setp(legendText, fontsize=18)
legend = plt.gca().get_legend()
legend.draw_frame(False)
fig = plt.gcf()
fig.set_size_inches(7,5)
plt.tight_layout()
pltFile = 'u-residual.png'
plt.savefig(pltFile, format='png')
plt.close()

print "%s DONE!!" % (pltFile)
