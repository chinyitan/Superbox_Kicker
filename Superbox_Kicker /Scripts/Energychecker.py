import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import numpy as np
import struct



###################################  # Settings

#N-body units
Munit=5.6e10  #Solar Mass
runit=3.5  # kpc
tunit = 0.013 # Gyrs

#Input Settings
CONTfiletitle ="m1e3r05_30.CONT"      #Initial CONT file, not modified by superbox
filenametitle= "m1e3r05_3-g01."
step=14400         #Step Length		
istep= 14400           #Initial Step
fstep=72000+step  # FInal Step (Last step included)


npartot= 50000  #Number of particles
nvar=17         #Number of variables
dlen= npartot* nvar  #length of Snapshot



#Dark matter halo data (N-body units)
gamma=0       
Mh=1.4e9 /Munit    #Msol
rs = 1.6 /runit    # kpc
Mmac=20/Munit      #Msol

Simtime =1 /tunit  # Length of Evolution of Cluster (Gyrs)


################################### 


		

def binreader(filename1):  #Read Snapshot for energy
	data = open(filename1, 'rb').read(4*dlen)
	values = struct.unpack('f'*dlen, data)
	del data	
	for i in range (0,(int(len(values)/nvar))):
#		R2.append((values[nvar*i]**2+values[nvar*i+1]**2+values[nvar*i+2]**2)/2)   	   #Distance squaredof all particles
#		V2.append((values[nvar*i+3]**2+values[nvar*i+4]**2+values[nvar*i+5]**2)/2)     #Velocity Square of all particles 
		E.append((values[nvar*i+12]+values[nvar*i+13])/2)                              #Sum all energy contribution , converted to N-body Unit
		EU.append((values[nvar*i+13])/2)                                               #Sum of all potential,converted to N-body unit
		EK.append((values[nvar*i+12])/2)                                               #Sum of kinetic,converted to N-body unit
		Vk.append(values[nvar*i+16]/2)                                                 #Sum velocity kicks contrubution,converted to N-body Unit 


def binreaderCONT(filenameCONT):      #Read CONT file for initial velocities
	floatdata=[]
	r=[]
	Dat_file=open (filenameCONT,"rb") 
	Datdata = list(Dat_file.read())
	for i in range (0,(int(len(Datdata)/4))):
		hexval=((Datdata[4*i+3])*256**3)+((Datdata[4*i+2])*256**2)+((Datdata[4*i+1])*256) +((Datdata[4*i]))
		binval=f"{hexval:0>32b}"
		f = int(str(binval), 2)
		floatdata.append(struct.unpack('f', struct.pack('I', f))[0])	        
	for i in range (0,(int((len(floatdata)-120)/6))):
		#Position velocity 6-vectors
		r.append([])
		r[i].append(floatdata[6*i+120])                  #Append positions of particles
		r[i].append(floatdata[6*i+1+120])
		r[i].append(floatdata[6*i+2+120])
		r[i].append(floatdata[6*i+3+120]/math.sqrt(2))
		r[i].append(floatdata[6*i+4+120]/math.sqrt(2))   # Append velocities of particle and convert in Nbody unit
		r[i].append(floatdata[6*i+5+120]/math.sqrt(2))	
	timestep = floatdata[63]*math.sqrt(2)    #Read timesteps length in Nbody units
	return r,timestep



De={}       #Dictionaries
Deu={}
Dek={}
Dkv2={}
timestep = 0 #Calculated by Script

rCONT,timestep=binreaderCONT(CONTfiletitle )  #Read CONT file



for tp in range (istep,fstep,step):    # Read DAT file
	E=[]
	EU=[]
	EK=[]
	Vk=[] 
	filename1=str(filenametitle+(format(tp, '07d'))) 
	binreader(filename1)
	De[tp]=E       #Compiling data to Dictionaries
	Deu[tp]=EU
	Dek[tp]=EK	
	Dkv2[tp]=Vk      

TimeSI=[]  #Timesteps
TimeNb=[]

SupEngy=[]  #Simulations
SupVar=[] 
SupPot=[]
SupKin=[]

ThEngy=[]  #Predictions


#Obtain Superbox Energies

for tp in range(istep,fstep,step):   #Time for each snapshot
	TimeNb.append(tp*timestep)
	TimeSI.append(tp*timestep*tunit)
	summ=0
	summ2=0
	summ3=0
	summ4=0
	
	for i in range (len(De[tp])):
		summ+=(De[tp][i]-De[istep][i])        # Sum(E-E0)
		summ2+=(De[tp][i]-De[istep][i])**2    # Sum((E-E0)^2)  
		summ3+=(Deu[tp][i]-Deu[istep][i])     # Sum(U-U0)
		summ4+=(Dek[tp][i]-Dek[istep][i])     # Sum(T-T0)

	SupEngy.append(summ/len(De[tp]))	                             # Calculate <(E-E0)>
	SupVar.append((summ2/len(De[tp]))-(summ/len(De[tp]))**2)         # Calculate <(E-E0)^2> - (<(E-E0)>)^2
	SupPot.append(summ3/len(Deu[tp]))
	SupKin.append(summ4/len(Deu[tp]))



	
#Calculate Energy Gradient for 1st Snapshot  

sumv=0
for i in range(len(Dkv2[step])):       #Can use istep if istep contains kicks
	sumv+=Dkv2[step][i]
sumv=sumv/(len(Dkv2[step]))
Engra2=sumv/timestep
print("Semi-Analytic gradient  :  " +str(Engra2))



######   Calculate Theoretical Predictions


def rhor(r):      #Dehnen density profile
    rhoval = (3-gamma)*Mh*rs/(4*math.pi)/(r**gamma*(r+rs)**(4-gamma))
    return rhoval
 
def sigma2r(r):    #Dehnen velcoity dispertion (gamma=0)
    sigma2 = Mh/30*(rs+6*r)/(rs+r)**2 
    return sigma2 


#Create linspace for time
TimeNb2=np.linspace(0,Simtime,1000)
Time2SI=[]
for i in range(len(TimeNb2)):
	Time2SI.append(TimeNb2[i]*tunit)


#Calculate <dv2>
vdv2 = 3*sigma2r(0) 
dv2= math.sqrt((32*math.pi**3) /(3*vdv2))*(8.2-1.9)*Mmac*rhor(0) #t taken out   
print("Theoretical gradient    :  " +str(dv2))

#Calculate Theoretical Energy Variation
for i in range (len(TimeNb2)):
	Eval= (1/2*dv2  *(TimeNb2[i]-TimeNb[0])) 
	ThEngy.append(Eval)


ThPot=[2*x for x in ThEngy]      # Predicted potential
ThKin=[-x for x in ThEngy]


#Plot graphs		


plt.plot(TimeSI,SupEngy,color='k',label="Total Energy change (E-E0)")
plt.plot(TimeSI,SupPot,color='b',label=" Potential Energy change (U-U0)")
plt.plot(TimeSI,SupKin,color='r',label="Kinetic Energy change (U-U0)")


plt.plot(Time2SI,ThEngy,color='k',ls="--")
plt.plot(Time2SI,ThPot,color='b',ls="--")
plt.plot(Time2SI,ThKin,color='r',ls="--")


plt.xlabel("Time, t (Gyrs)", fontsize=13)
plt.ylabel("$E-E_0$ (N-body units) ", fontsize=13)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.legend()
plt.show()




#########################   #Looking at dE/dt(t) instead of E-E0(t)

DTimeSI=[]

SupDEngy=[]  #Simulations 
SupDPot=[]
SupDKin=[]



Tstep=timestep*step

for tp in range(istep+step,fstep,step):   #Time for each snapshot

	summ=0
	summ2=0
	summ3=0

	
	for i in range (len(De[tp])):
		summ+=(De[tp][i]-De[tp-step][i])      # E_i-E_(i-1)
		summ2+=(Deu[tp][i]-Deu[tp-step][i])  
		summ3+=(Dek[tp][i]-Dek[tp-step][i])

	SupDEngy.append(summ/(len(De[tp])*Tstep))	#Dividing by timestep and dividing by particle number
	SupDPot.append(summ2/(len(Deu[tp])*Tstep))
	SupDKin.append(summ3/(len(Deu[tp])*Tstep))	

for i in range(len(TimeSI)-1):              #Taking average for timestep
	DTimeSI.append((TimeSI[i+1]+TimeSI[i])/2)


#Calculate Energy Dispersion
ThDEngy = 1/2*dv2                  # Predicted potential difference
ThDPot  = 2*ThDEngy
ThDKin  = -ThDEngy



#Plot graphs		


plt.plot(DTimeSI,SupDEngy,color='k',label="Total Energy change (E-E0)")
plt.plot(DTimeSI,SupDPot,color='b',label="Potential Energy change (U-U0)")
plt.plot(DTimeSI,SupDKin,color='r',label="Kinetic Energy change (T-T0)")

plt.axhline(ThDEngy,color='k',ls="--")
plt.axhline(ThDPot,color='b',ls="--")
plt.axhline(ThDKin,color='r',ls="--")

plt.xlabel("Time t (Gyrs)", fontsize=13)
plt.ylabel(" $\\delta E / \\delta t $ (N-body units) ", fontsize=13)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

plt.legend()
plt.show()



