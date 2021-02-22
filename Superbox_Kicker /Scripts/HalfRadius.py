import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import math



#################  # Settings  

"""
ODE Units:
Length pc
Mass   Msol
Time   Gyrs
"""


#Age
tfinal= 1 #Gyrs

# DM Halo Paramter
Mhalo= 1.4e9 
rhalo=  1.6    #kpc

#Cluster paramter
Ms = 1350.     #Solar Mass
y0 = 2.8368   # Initial Condition - Parsec 

#MACHOs
fdm=1.        #Dark Matter Density

#Initial Alpha and Beta
alpha= 0.375
beta = 37.5

#Settings
binsize=7200   #Binning Size


#Constants
G=6.67408e-11 *6.731e13  #SI to ODE units
sigunit = 1000  *1.0222    #kms^-1 to pcGyrs^-1


##################


#Read files
time=[]	
lag=[]	

inputfile =open("lagr.g01","r")
inputdata = inputfile.readlines()
for i in range (0,len(inputdata)):
	x = (float(inputdata[i][9:14]))*(10**float(inputdata[i][15:18]))   #Read time of step
	y = (float(inputdata[i][72:78]))                                   #Read Lagragia radius (L5), half-light radius
	time.append(x*0.013*(2**0.5))
	lag.append(y*3.5*1000)


#Binning
ctime=[]
clag=[]
cvar=[]

for i in range(int(len(time)/binsize)):
	yvalue=sum(lag[i*binsize:(i+1)*binsize ])/binsize
	yvar=np.var(lag[i*binsize:(i+1)*binsize ])
	ctime.append(time[int((i+0.5)*binsize)])
	clag.append(yvalue)
	cvar.append(math.sqrt(yvar + (0.175)**2 ))

### Theoretical Predictions


def lnlambda(r,mmac,sigma):
#	lnlam = math.log(r*sigma**2/(G*(1+mmac)))    #Couloumb Logarithm from Bradnt(2016)
	lnlam = 8.2-1.9                              #Couloumb Logarithm from Penarrubia(2019)
	return lnlam


def model(y,t,mmac,rho,sigma):                   #Differnetial Equation from Brandt (2016)
    k=math.sqrt(32*math.pi**3/3)*G*mmac/sigma*lnlambda(y,mmac,sigma)*fdm   
    dydt = k*((alpha*Ms/(rho*y**2))+2*beta*y)**-1
    return dydt

t = np.linspace(0,tfinal,2000)   #Time points

#Calcualte rho0 and sigma for Eridanus II, You can choose a different rho0 and sigma if you want
rho0 = 0.037*(1+504/(rhalo*1000))**3
sigma = 2.33*(1+(rhalo*1000)/504)**(3/2)*((rhalo*1000)/504)**(-1/2)*math.sqrt(3)*sigunit


y4 = odeint(model,y0,t,args=(20,rho0,sigma))   #Solving ODE

plt.errorbar(ctime,clag,cvar,label="Superbox Simulations")                     #Plotting Results
plt.plot(t,y4,ls='dashed',label="Theoretical Predictions",color="tab:blue")
plt.ylabel("  Half-light radius, r$_h$ (pc)", fontsize=13)
plt.xlabel(" Time, t (Gyrs)", fontsize=13)   
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.legend()
plt.show()


