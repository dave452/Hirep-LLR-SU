import matplotlib.pyplot as plt
import sys
import numpy as np

faddr=sys.argv[1]

rawdata = np.genfromtxt( faddr, usecols=(0,1,2,3,4), names=['beta','t','st','tsym','stsym'])

X = rawdata['beta']
Y = rawdata['t']
sY = rawdata['st']

Y2 = rawdata['tsym']
sY2 = rawdata['stsym']

plt.figure()
plt.xlim(7.5,8.5)
plt.ylabel(r'$\sqrt{8 t_c}a$')
plt.xlabel(r'$\beta$')
#plt.title(r'$\sqrt{t_c}(\beta)$')
plt.errorbar(X,Y,yerr=sY,fmt='o',linestyle='None', label='Plaquette')
plt.errorbar(X,Y2,yerr=sY2,fmt='o',linestyle='None', label='Clover')
plt.legend()
plt.grid(True)
plt.show()
