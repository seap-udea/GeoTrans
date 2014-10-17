from geotrans import *
from physics import *
from constants import *
fig=plt.figure(figsize=(8,8))
ax=fig.gca()
T="\t"
V=1

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
iorb=89.8*DEG
wp=45.0*DEG

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
fcen=mod(270*DEG-wp,360*DEG) #Central true anomaly (rad)
Ecen=mod(2*arctan(sqrt((1-ep)/(1+ep))*tan(fcen/2)),360*DEG)
Mcen=mod(Ecen-ep*sin(Ecen),360*DEG)
#TIME AND POSITION
tcen=Mcen/norb
rcen=ellipseRadiusE(ap,ep,fcen) #r central (km)
#PROJECTED POSITION
Pcen=dot(Mos,AR3(rcen*cos(fcen),rcen*sin(fcen),0))
#IMPACT PARAMETER
Borb=Pcen[1]/Rstar

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
#SCRIPT
#########################################
df=Rstar/rcen
print df

#########################################
#PLOT
#########################################
#PLOT FIGURES
plotEllipse(ax,Star,patch='true',fc='y',ec='none')
plotEllipse(ax,Planet,patch='true',fc='b',ec='none')
plotEllipse(ax,Ringe,color='k')
#plotEllipse(ax,Ringi,color='k')
plotPoint(ax,toPoint(AR(Pcen[0],Pcen[1])))

#PLOT ORBIT
fs=linspace(fcen-1*DEG,fcen+1*DEG,100)
rs=ellipseRadiusE(ap/Rstar,ep,fs)
xs=rs*cos(fs)
ys=rs*sin(fs)
zs=zeros_like(xs)
rs=array([dot(Mos,AR3(x,y,z)) for x,y,z in zip(xs,ys,zs)])
ax.plot(rs[:,0],rs[:,1],'b-')

#PLOT RINGS
fs=linspace(0,2*pi,100)
xs=Re*cos(fs)
ys=Re*sin(fs)
zs=zeros_like(xs)
rs=array([dot(Mrs,AR3(x,y,z))+Pcen/Rstar for x,y,z in zip(xs,ys,zs)])
for i in xrange(len(rs)):
    c='k'
    if (rs[i,2]-Pcen[2]/Rstar)<0:c='r'
    ax.plot([rs[i,0]],[rs[i,1]],'o',markeredgecolor='none',color=c,markersize=2)

#PLOT AXIS
rp=3*Rp*rz+Pcen/Rstar
ax.plot([Pcen[0]/Rstar,rp[0]],
        [Pcen[1]/Rstar,rp[1]],'k-',linewidth=3)

rp=3*Rp*ry+Pcen/Rstar
ax.plot([Pcen[0]/Rstar,rp[0]],
        [Pcen[1]/Rstar,rp[1]],'c-',linewidth=3)

rp=3*Rp*rx+Pcen/Rstar
ax.plot([Pcen[0]/Rstar,rp[0]],
        [Pcen[1]/Rstar,rp[1]],'g-',linewidth=3)

rp=3*Rp*AR3(cos(teff),sin(teff),0)+Pcen/Rstar
ax.plot([Pcen[0]/Rstar,rp[0]],
        [Pcen[1]/Rstar,rp[1]],'k-',linewidth=3)

#########################################
#DECORATION
#########################################
xrange=Pcen[0]/Rstar
yrange=Pcen[1]/Rstar
range=0.0*1.5+0.0*ap/Rstar+10*Rp
ax.set_xlim((-range+xrange,range+xrange))
ax.set_ylim((-range+yrange,range+yrange))
ax.grid()
fig.savefig("plots/transit-real.png")
