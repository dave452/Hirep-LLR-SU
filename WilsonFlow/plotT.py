import matplotlib.pyplot as plt
import sys
import numpy as np

faddr=sys.argv[1]

rawdata = np.genfromtxt( faddr, usecols=(0,1,2), names=['beta','t','st'])

X = rawdata['beta']
Y = np.sqrt(rawdata['t'])
sY = 0.5/np.sqrt(rawdata['t'])*rawdata['st']


plt.figure()
plt.xlim(7.5,8.5)
plt.ylabel(r'$\sqrt{t_c}a$')
plt.xlabel(r'$\beta$')
plt.title(r'$\sqrt{t_c}(\beta)$')
plt.errorbar(X,Y,yerr=sY,fmt='o',linestyle='None')
plt.grid(True)
plt.show()
