from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import sys

faddr=sys.argv[1]
ymin=0.25
ymax=0.32
beta = np.float(sys.argv[2])

t2Escale = 0.3
avr=0.
avr2=0.
rawdata = np.genfromtxt(faddr, usecols =(0,1,2,3,4,5), names=['Nconf', 't', 'E', 't2E', 'Esym', 't2Esym', 'TC'])
nconf = np.int(np.max(rawdata['Nconf'])-1)

tc = np.empty(0,dtype='f8');

for nc in range(1,nconf):
 data = np.compress( rawdata['Nconf']==nc , rawdata)
 X = data['t']
 Y = data['t2Esym']
 if( np.max(Y) > 0.3 ):
   f2 = interp1d(Y,X,kind='cubic')
   f2scale = f2(0.3)
   tc = np.append(tc,f2scale)

#print len(tc)
#print tc

BinSize = np.int(sys.argv[3])
#print len(tc)/BinSize
bindata = tc[:(len(tc) / BinSize) * BinSize].reshape(-1, BinSize).mean(axis=1)
#
NBin = len(bindata)
Nboot = np.int(sys.argv[4])
bootvar = np.empty(Nboot,dtype='f8')
#binarray = np.empty(NBin,dtype='f8')
#
#binavr=0
#bincounter=0
#for i in range(0,len(tc)):
# #print "tc", i, tc[i]
# binavr += tc[i]
# bincounter +=1
# #print i,
# if( bincounter  % BinSize == 0 ):
#  #print "binned value", binavr/BinSize
#  binarray[i/BinSize] = binavr/BinSize
#  binavr=0
#print binarray 
#
for i in range(0,len(bootvar)):
 #print i, binarray[np.random.randint(0,NBin)]
 bootvar[i] = np.mean(np.sqrt(8.*np.random.choice(bindata, len(bindata))),dtype=np.float32)
 
# 
# 
print np.mean(bootvar,dtype=np.float32),
print np.sqrt(np.var(bootvar,ddof=1,dtype=np.float32))
##print binarray
#print np.average(tc), np.std(tc)
#print np.average(binarray), np.std(binarray)
##print np.average(bootvar), np.std(bootvar)
#  
#
#	
#
#
##avr /= nconf
##avr2 /= nconf
##avr2 -= avr*avr 
##err = np.sqrt(avr2/(nconf-1))
#
##plt.plot(X,Y,'o')
##plt.plot(Xnew, f1(Xnew), '^')
##plt.plot(Xnew, f2(Xnew))
##plt.grid(True)
##plt.show()
