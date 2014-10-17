from geotrans import *
from physics import *
from constants import *
fig=plt.figure(figsize=(8,8))
ax=fig.gca()
T="\t"
V=0

#########################################
#PRIMARY PARAMETERS
#########################################
#//////////////////////////////
#STAR
#//////////////////////////////
Mstar=1.0*MSUN
Rstar=1.0*RSUN
#//////////////////////////////
#PLANET
#//////////////////////////////
Mp=1.0*MSAT
Rp=1.0*RSAT
fp=1.0 #Oblateness
#//////////////////////////////
#RINGS
#//////////////////////////////
fe=RSAT_ARING/RSAT #Exterior ring (Rp)
fi=RSAT_BRING/RSAT #Interior ring (Rp)
ir=30.0*DEG #Ring inclination
phir=190.0*DEG #Ring roll angle
tau=1.0 #Opacity
#//////////////////////////////
#ORBIT
#//////////////////////////////
ap=1.0*AU
ep=0.6
iorb=89.9*DEG
wp=0.0*DEG

#########################################
#DERIVATIVE PARAMETERS
#########################################
#//////////////////////////////
#ROTATION MATRICES
#//////////////////////////////
#From orbit to sky
Mi=rotMat([1,0,0],-iorb)
Mw=rotMat([0,0,1],wp)
Mos=dot(Mi,Mw)
#From rings to sky
Mpr=rotMat([0,0,1],phir)
Mir=rotMat([1,0,0],-ir)
Mro=dot(Mpr,Mir)
Mrs=dot(Mi,Mro)

#//////////////////////////////
#STAR
#//////////////////////////////
rhostar=DENSITY(Mstar,Rstar)

#//////////////////////////////
#PLANET
#//////////////////////////////
Rp=Rp/Rstar

#//////////////////////////////
#RING
#//////////////////////////////
Ri=fi*Rp
Re=fe*Rp

#APARENT ORIENTATION
rx=dot(Mrs,[1.0,0.0,0.0])
ry=dot(Mrs,[0.0,1.0,0.0])
rz=dot(Mrs,[0.0,0.0,1.0])
ieff=arccos(abs(dot(rz,[0,0,1])))
teff=-sign(rz[0])*ARCTAN(abs(rz[0]),abs(rz[1]))

#//////////////////////////////
#ORBIT
#//////////////////////////////
#PERIOD
Porb=2*pi*sqrt(ap**3/(GCONST*(Mstar+Mp))) #Period (s)
norb=2*pi/Porb
#ANOMALIES
fcen=270*DEG-wp #Central true anomaly (rad)
Ecen=2*arctan(sqrt((1-ep)/(1+ep))*tan(fcen/2))
Mcen=Ecen-ep*sin(Ecen)
#print "fcen,Mcen,Ecen = ",fcen*RAD,Mcen*RAD,Ecen*RAD
#TIME AND POSITION
tcen=Mcen/norb
rcen=ellipseRadiusE(ap,ep,fcen) #r central (km)
rpcen=AR3(rcen*cos(fcen),rcen*sin(fcen),0)
#print "rpcen f = ",rpcen
#print ap*cos(Ecen)-ap*ep
rpcen=AR3(ap*cos(Ecen)-ap*ep,ap*sqrt(1-ep**2)*sin(Ecen),0)
#print "rpcen E = ",rpcen
#print "r orbit cen = ",rpcen
#print "Mos = ",Mos
#PROJECTED POSITION
Pcen=dot(Mos,rpcen)
#print "r sky cen = ",Pcen
#IMPACT PARAMETER
Borb=Pcen[1]/Rstar
#ESTIMATED TRANSIT DURATION
if abs(Borb)>1:
    print "No transit."
    exit(1)
df=Rstar*sqrt(1-Borb**2)/rcen
f15=fcen-df
r15=ellipseRadiusE(ap,ep,f15) 
P15=dot(Mos,AR3(r15*cos(f15),r15*sin(f15),0))/Rstar
t15=timeOrbit(ep,norb,f15)
f35=fcen+df
r35=ellipseRadiusE(ap,ep,f35) 
t35=timeOrbit(ep,norb,f35)
P35=dot(Mos,AR3(r35*cos(f35),r35*sin(f35),0))/Rstar
dtrans=t35-t15

#########################################
#REPORT
#########################################
if V:
    print "Star primary:"
    print T,"Ms = %e kg"%Mstar
    print T,"Rs = %e kg"%Rstar
    print "Planet primary:"
    print T,"Mp = %e kg = %e Mstar"%(Mp,Mp/Mstar)
    print T,"Rp = %e kg = %e Rstar"%(Rp,Rp/Rstar)
    print "Rings primary:"
    print T,"fi,fe = %e,%e Rp"%(fi,fe)
    print T,"Inclination (orbit) = %.1f deg"%(ir*RAD)
    print T,"Roll (orbit) = %.1f deg"%(phir*RAD)
    print T,"Opacity = %.2f"%(tau)
    print "Orbit primary:"
    print T,"ap = %e km = %e AU = %e Rstar"%(ap,ap/AU,ap/Rstar)
    print T,"Eccentricity = %.2f"%(ep)
    print T,"Inclination (visual) = %.2f deg"%(iorb*RAD)
    print T,"Periapsis argument = %.2f deg"%(wp*RAD)
    print
    print "Star derivative:"
    print T,"Star density = %e kg/m^3 = %e rhosun"%(rhostar,rhostar/RHOSUN)
    print "Planetary derivative:"
    print T,"Radius (relative) = %e Rstar"%(Rp)
    print "Rings derivative:"
    print T,"Internal ring (relative) = %.2f Rstar"%(Ri)
    print T,"External ring (relative) = %.2f Rstar"%(Re)
    print T,"Apparent inclination = %.2f deg"%(ieff*RAD)
    print T,"Apparent roll = %.2f deg"%(teff*RAD)
    print "Orbit derivative:"
    print T,"Period = %e s = %e h = %e d = %e yr"%(Porb,Porb/HOUR,Porb/DAY,Porb/YEAR)
    print T,"Mean Angular velocity = %e rad/s = %e Rstar/s = %e Rp/s"%(norb,norb*rcen/Rstar,norb*rcen/(Rp*Rstar))
    print T,"Central true anomaly = %e deg"%(fcen*RAD)
    print T,"Central eccentric anomaly = %e deg"%(Ecen*RAD)
    print T,"Central mean anomaly = %e deg"%(Mcen*RAD)
    print T,"Central radius = %e km = %e AU = %e Rstar"%(rcen,rcen/AU,rcen/Rstar)
    print T,"Impact parameter = %e Rstar"%(Borb)
    print T,"Central time = %e s = %e Porb"%(tcen,tcen/Porb)
    print T,"Estimated t1.5 = %e s = %e Porb"%(t15,t15/Porb)
    print T,"Estimated t3.5 = %e s = %e Porb"%(t35,t35/Porb)
    print T,"Estimated transit duration = %e s = %e h"%(dtrans,dtrans/3600.0)

"""
Roll angle is taken in such a way that phir=0 defines the summer
solstice for the northern hemisphere, phir=90 is the spring equinox
for the northern hemisphere
"""

#########################################
#FIGURES
#########################################
Star=Figure(AR(0,0),1,1,0,'Star')
C=AR(0.0,Borb)
Planet=Figure(C,Rp,Rp,0,'Planet')
Ringe=Figure(C,Re,Re*cos(ieff),teff,'Ringext')
Ringi=Figure(C,Ri,Ri*cos(ieff),teff,'Ringint')

#########################################
#TRANSIT
#########################################
POrbit=Orbit(ap/Rstar,ep,Porb,Mos)
"""
t=tcen+dt/2
At=transitAreaTime(t,POrbit,Planet,Ringe,Ringi)
print At
"""
Ar,Ps=ringedPlanetArea(Planet,Ringe,Ringi)
t1,A,dt,nfun=contactTime(tcen-dtrans/2,dtrans/(Rstar/Rp/2),POrbit,Planet,Ringe,Ringi,
                         Ar,Borb,Side=OUTSIDE,Phase=INGRESS,tola=1E-6)
t2,A,dt,nfun=contactTime(tcen-dtrans/2,dtrans/(Rstar/Rp/2),POrbit,Planet,Ringe,Ringi,
                         Ar,Borb,Side=INSIDE,Phase=INGRESS,tola=1E-6)
t3,A,dt,nfun=contactTime(tcen+dtrans/2,dtrans/(Rstar/Rp/2),POrbit,Planet,Ringe,Ringi,
                         Ar,Borb,Side=INSIDE,Phase=EGRESS,tola=1E-6)
t4,A,dt,nfun=contactTime(tcen+dtrans/2,dtrans/(Rstar/Rp/2),POrbit,Planet,Ringe,Ringi,
                         Ar,Borb,Side=OUTSIDE,Phase=EGRESS,tola=1E-6)

#########################################
#PLOT
#########################################
#PLOT FIGURES
plotEllipse(ax,Star,patch='true',fc='y',ec='none')
plotEllipse(ax,Planet,patch='true',fc='b',ec='none')
plotEllipse(ax,Ringe,color='k')
plotEllipse(ax,Ringi,color='k')
#plotPoint(ax,toPoint(AR(Pcen[0],Pcen[1])))
#plotPoint(ax,toPoint(AR(P15[0],P15[1])))
#plotPoint(ax,toPoint(AR(P35[0],P35[1])))

#PLOT ORBIT
fs=linspace(fcen-1*DEG,fcen+1*DEG,100)
rs=ellipseRadiusE(ap/Rstar,ep,fs)
xs=rs*cos(fs)
ys=rs*sin(fs)
zs=zeros_like(xs)
rs=array([dot(Mos,AR3(x,y,z)) for x,y,z in zip(xs,ys,zs)])
ax.plot(rs[:,0],rs[:,1],'b-')

#########################################
#DECORATION
#########################################
xrange=Pcen[0]/Rstar
yrange=Pcen[1]/Rstar
range=0.0*1.5+0.0*ap/Rstar+15*Rp
ax.set_xlim((-range+xrange,range+xrange))
ax.set_ylim((-range+yrange,range+yrange))
ax.grid()
fig.savefig("plots/transit-real.png")
