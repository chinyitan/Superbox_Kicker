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
G = 1

#Input Settings
CONTfiletitle ="m1e3r05_30.CONT"   #Initial CONT file, not modified by superbox
filenametitle= "m1e3r05_3-g01."
step=14400        #Step Length		
istep=14400        #Initial Step
fstep=72000 +step  # FInal Step (Last step included)


npartot= 50000  #Number of particles
nvar=17         #Number of variables
dlen= npartot* nvar  #length of Snapshot



# DM Halo Paramter
Mhalo= 1.4e9 /Munit   #Solar Mass
rhalo=  1.6  /runit   #kpc
rho0 = Mhalo/(4/3*math.pi*rhalo**3)   


#Cluster paramter
Ms = 1350. /Munit    #Solar Mass


###################################  

		

def binreader(filename1):  #Read Snapshot for energy
	data = open(filename1, 'rb').read(4*dlen)
	values = struct.unpack('f'*dlen, data)
	del data	
	E1val=0
	E2val=0    
	for i in range (0,(int(len(values)/nvar))): 
		r1.append(math.sqrt(values[nvar*i]**2+values[nvar*i+1]**2+values[nvar*i+2]**2))   #Append position of particles, to calculate r_half      
		E1val+=((values[nvar*i+13]-values[nvar*i+15])/2 )   #Calcluate U_self=U_tot-U_halo (Converted to N-body)
		E2val+=(values[nvar*i+15]/2 ) 
	E1.append(E1val)
	E2.append(E2val)

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
	timestep = floatdata[63]*math.sqrt(2)    #Read timestep length in Nbody units
	return r,timestep




Dr={}
De={}       #Dictionaries
De1={}
De2={}
Dv2={}
Dkv2={}
timestep = 0 #Calculated by Script

rCONT,timestep=binreaderCONT(CONTfiletitle )  #Read CONT file



for tp in range (istep,fstep,step):    # Read DAT file
	E1=[]
	E2=[]
	r1=[]     
	filename1=str(filenametitle+(format(tp, '07d'))) 
	binreader(filename1) 
	De1[tp]=E1    #Compiling data to Dictionaries
	De2[tp]=E2
	Dr[tp]=r1 

    

rhalflist=[]    
Alphalist=[]
Betalist=[]

for tp in range (istep,fstep,step):    
    rlist=Dr[tp]                     #Calcaulating r_half
    rlist2=sorted(rlist)
    rhalf=rlist2[25000]
    rhalflist.append(rhalf*3500)
    Alphalist.append(-De1[tp][0]*rhalf/(G*Ms) /50000)    #Calculating Alpha
    Betalist.append( ((De2[tp][0]/50000)-(-(G*Mhalo)/(2*rhalo)))/(G*rho0*rhalf**2))  #Calculating Beta


plt.scatter(rhalflist,Alphalist)
plt.plot(rhalflist,Alphalist,label="5 pc")
plt.xlabel("Radius,r, (pc)", fontsize=13)
plt.ylabel('Alpha,$\\alpha $', fontsize=13)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.show()


plt.scatter(rhalflist,Betalist)
plt.plot(rhalflist,Betalist,label="5 pc")
plt.xlabel("Radius,r, (pc)", fontsize=13)
plt.ylabel("Beta,$\\beta$ ", fontsize=13)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.show()

