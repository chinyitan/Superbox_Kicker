import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import numpy as np
import struct
import time
from scipy.optimize import curve_fit


runit=3.5
tunit=262.
Munit=5.6e10

radmax=100/(runit*1000)
lradmin=-4
lradmax=-1.5

bins=500
bins2=50
	

def binreader(filename1):
	r=[]    
	data = open(filename1, 'rb').read(4*dlen)
	values = struct.unpack('f'*dlen, data)
	del data	
	
	for i in range (0,(int(len(values)/Numvar))):
		#Position velcoity 6-vectors
		r.append([])
		r[i].append(values[Numvar*i])
		r[i].append(values[Numvar*i+1])
		r[i].append(values[Numvar*i+2])
	return r
        
		
def binreaderCONT(filename1):
	floatdata=[]
	r=[]
	Dat_file=open (filename1,"rb") 
	Datdata = list(Dat_file.read())
	for i in range (0,(int(len(Datdata)/4))):
		hexval=((Datdata[4*i+3])*256**3)+((Datdata[4*i+2])*256**2)+((Datdata[4*i+1])*256) +((Datdata[4*i]))
		binval=f"{hexval:0>32b}"
		f = int(str(binval), 2)
		floatdata.append(struct.unpack('f', struct.pack('I', f))[0])	
	for i in range (0,(int((len(floatdata)-120)/6))):
		#Position velcoity 6-vectors
		r.append([])
		r[i].append(floatdata[6*i+120])
		r[i].append(floatdata[6*i+1+120])
		r[i].append(floatdata[6*i+2+120])

	return r



Dr={}                         #Dictionarites of Values
Drsat={}
Numvar = 17
dlen=50000*17   #Numvar







Dr[0]=binreader("m1e3r05_3-g01.0072000")   # Read output file , assign key (tp=0)

#Sorting Starts Here 
Radius=[]
RadiusSI=[]
logRadius=[]
logRadiusSI=[]


# Radius points (Linear)
dr=radmax/bins	           
for i in range (bins):
	Radius.append(dr*i+dr/2)
	RadiusSI.append((dr*i+dr/2)*3500) 

# Radius points (log)
logbins=np.logspace(lradmin,lradmax, bins2+1)
for i in range (bins2):
	logRadius.append((logbins[i]+logbins[i+1])/2)
	logRadiusSI.append(math.log10(((logbins[i]+logbins[i+1])/2)*runit*1000))
    



def RadSort(tp):       # Obtain 3D-Density per bin (Linear)
	RList=[]
	Binlist=[]
	Density=[]
	for i in range(len(Dr[tp])):
		RList.append(math.sqrt((Dr[tp][i][0])**2+(Dr[tp][i][1])**2+(Dr[tp][i][2])**2))
	for i in range(bins):
		Binlist.append(len(list(x for x in RList if dr*i <= x < dr*(i+1)))) 
	for i in range(len(Binlist)):
		rho = (Binlist[i])/(4*math.pi*Radius[i]**2*dr)
		rhoSI = rho/(3500**3)      #Convert to SI
		Density.append(rhoSI)
	return Density



def SigmaSort(tp):    # Obtain Surface Density per bin (Linear)
	RList=[]
	Binlist=[]
	Density=[]
	for i in range(len(Dr[tp])):
		RList.append(math.sqrt((Dr[tp][i][0])**2+(Dr[tp][i][1])**2))
	for i in range(bins):
		Binlist.append(len(list(x for x in RList if dr*i <= x < dr*(i+1)))) 
	for i in range(len(Binlist)):
		rho = (Binlist[i])/(2*math.pi*Radius[i]*dr)
		rhoSI = rho/(3500**2)      #Convert to SI
		Density.append(rhoSI)
	return Density



def logRadSort(tp):    # Obtain 3D-Density per bin (binning in Log)
	RList=[]
	Binlist=[]
	lnDensity=[]
	for i in range(len(Dr[tp])):
		RList.append(math.sqrt((Dr[tp][i][0])**2+(Dr[tp][i][1])**2+(Dr[tp][i][2])**2))
	for i in range(bins2):
		Binlist.append(len(list(x for x in RList if logbins[i]<= x < logbins[i+1]))) 
	for i in range(len(Binlist)):
		rho = (Binlist[i])/(4*math.pi*logRadius[i]**2*(logbins[i+1]-logbins[i]))
		rhoSI = rho/((runit*1000)**3)      #Convert to pc
		if rhoSI <= 0:
			lnDensity.append(float("-inf"))
		else:
			lnDensity.append(math.log10(rhoSI))
	return lnDensity



def logSigmaSort(tp):    # Obtain Surface Density per bin (Log)
	RList=[]
	Binlist=[]
	lnDensity=[]
	for i in range(len(Dr[tp])):
		RList.append(math.sqrt((Dr[tp][i][0])**2+(Dr[tp][i][1])**2))#+(Dr[tp][i][2])**2))
	for i in range(bins2):
		Binlist.append(len(list(x for x in RList if logbins[i]<= x < logbins[i+1]))) 
	for i in range(len(Binlist)):
		rho = (Binlist[i])/(2*math.pi*logRadius[i]*(logbins[i+1]-logbins[i]))
		rhoSI = rho/((runit*1000)**2)      #Convert to pc
		if rhoSI <= 0:
			lnDensity.append(float("-inf"))
		else:
			lnDensity.append(math.log10(rhoSI))
	return lnDensity



Sigmabins=SigmaSort(0)         # Call function  
lnSigmabins=logSigmaSort(0) 




#Fitting with Gaussian 

expval=math.exp(1)

def gaus(x,a,sigma,y0):
    return a*expval**(-(x/sigma)**2)+y0

popt,pcov = curve_fit(gaus,RadiusSI,Sigmabins)
print("Fitted Parameters (a,sigma,yo)  :  " +str(popt))

#Coamparing with Theoretical Initial Profile

#Initial Stellar Profile
M = 1350   #Msol
a = 0.78   #pc
NumStar = 50000

Mstar=M/NumStar   #Mass of star for Normalization to correct number of stars

def Dehnen(x):      # Density Profile for Cored Dehnen model
	rho=3*M/(4*math.pi)*a/(x+a)**4
	return rho


def DehnenSurface(x):   # Surface Density Profile for Cored Dehnen model
	s=x/a
	if s>1:
		y=(s**2-1)**-0.5*np.arccos(s**(-1))
	elif s==1:
		y=1.0
	else:	
		y=(1-s**2)**-0.5*np.arccosh(s**(-1))
	sigma = M/(4*math.pi*a**2)/(s**2-1)**3 *(-2-13*s**2+3*s**2*(4+s**2)*y)
	return sigma

ThSigma=[DehnenSurface(x)/Mstar for x in RadiusSI]             # Adding Dehnen models for linespace
logThSigma=[math.log10(DehnenSurface(10**x)/Mstar) for x in logRadiusSI]   #Note: items in logRadiusSI have been converted into log



#Plotting



plt.plot(RadiusSI, Sigmabins,label="Superbox profile")
plt.plot(RadiusSI, gaus(RadiusSI, *popt),color="k",label="Fitted Gaussian profile")
plt.xlabel("Radius,r (pc)", fontsize=13)
plt.ylabel("Surface Brightness Density, $\\Sigma $ (pc$^{-2}$)", fontsize=13)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.legend()
plt.show()





plt.plot(logRadiusSI, lnSigmabins,label="Current Profile")
plt.plot(logRadiusSI, logThSigma,ls="dashed",label="Initial Profile")
plt.xlabel("log$_{10}$ Radius,r (pc)", fontsize=13)
plt.ylabel("log$_{10}$  Surface Brightness Density, $\\Sigma $ (pc$^{-2}$)", fontsize=13)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.legend()
plt.show()

