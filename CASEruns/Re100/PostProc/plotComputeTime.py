#!/usr/bin/env python
import matplotlib.pyplot as plt
import csv

errFile = 'computeTime.csv'

nList = []
errList = []
with open(errFile, 'rb') as csvfile:
   reader = csv.DictReader(csvfile)
   for row in reader:
      nList.append(float(row['n']))
      errList.append(float(row['computeTime']))

MinX = 15
MaxX = 85
MinY = 50.0
MaxY = 4000.0

p = plt.plot(nList,errList, 'k-o', markersize=15)
plt.setp(p, linewidth='3.0')

plt.axis([MinX,MaxX, MinY, MaxY])
plt.xscale('linear')
plt.yscale('log')
plt.xlabel('grid resolution, N', fontsize=22)
plt.ylabel('Compute time [sec]', fontsize=22)

plt.grid(True)
ax = plt.gca()
xlabels = plt.getp(ax, 'xticklabels')
ylabels = plt.getp(ax, 'yticklabels')
plt.setp(xlabels, fontsize=18)
plt.setp(ylabels, fontsize=18)
#plt.legend(
#          loc='upper right',
#          borderpad=0.25,
#          handletextpad=0.25,
#          borderaxespad=0.25,
#          labelspacing=0.0,
#          handlelength=2.0,
#          numpoints=1)
#legendText = plt.gca().get_legend().get_texts()
#plt.setp(legendText, fontsize=18)
#legend = plt.gca().get_legend()
#legend.draw_frame(False)
fig = plt.gcf()
fig.set_size_inches(7,5)
plt.tight_layout()
pltFile = 'computeTime.png'
plt.savefig(pltFile, format='png')
plt.close()

print "%s DONE!!" % (pltFile)
