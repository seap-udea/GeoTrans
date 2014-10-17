from geotrans import *
from physics import *
from constants import *
fig=plt.figure(figsize=(8,8))
ax=fig.gca()
T="\t"

#########################################
#PRIMARY PARAMETERS
#########################################
#//////////////////////////////
#STAR
#//////////////////////////////
Mstar=1.0*MSUN
Rstar=1.0*RSUN
print "Star primary:"
print T,"Ms = %e kg"%Mstar
print T,"Rs = %e kg"%Rstar

#//////////////////////////////
#PLANET
#//////////////////////////////
Mp=1.0*MSAT
Rp=1.0*RSAT
fp=1.0 #Oblateness
print "Planet primary:"
print T,"Mp = %e kg = %e Mstar"%(Mp,Mp/Mstar)
print T,"Rp = %e kg = %e Rstar"%(Rp,Rp/Rstar)

#//////////////////////////////
#RINGS
#//////////////////////////////
fe=RSAT_ARING/RSAT #Exterior ring (Rp)
fi=RSAT_BRING/RSAT #Interior ring (Rp)
ir=10.0*DEG #Ring inclination
"""
Roll angle is taken in such a way that phir=0 defines the summer
solstice for the northern hemisphere, phir=90 is the spring equinox
for the northern hemisphere
"""
phir=270.0*DEG #Ring roll angle
tau=1.0 #Opacity
print "Rings primary:"
print T,"fi,fe = %e,%e Rp"%(fi,fe)
print T,"Inclination (orbit) = %.1f deg"%(ir*RAD)
print T,"Roll (orbit) = %.1f deg"%(phir*RAD)
print T,"Opacity = %.2f"%(tau)

#//////////////////////////////
#ORBIT
#//////////////////////////////
ap=1.0*AU
ep=0.8
iorb=80.0*DEG
wp=90.0*DEG
print "Orbit primary:"
print T,"ap = %e km = %e AU = %e Rstar"%(ap,ap/AU,ap/Rstar)
print T,"Eccentricity = %.2f"%(ep)
print T,"Inclination (visual) = %.2f deg"%(iorb*RAD)
print T,"Periapsis argument = %.2f deg"%(wp*RAD)

#########################################
#DERIVATIVE PARAMETERS
#########################################
print
#//////////////////////////////
#STAR
#//////////////////////////////
rhostar=DENSITY(Mstar,Rstar)
print "Star derivative:"
print T,"Star density = %e kg/m^3 = %e rhosun"%(rhostar,rhostar/RHOSUN)

#//////////////////////////////
#PLANET
#//////////////////////////////
Rp=Rp/Rstar
print "Planetary derivative:"
print T,"Radius (relative) = %e Rstar"%(Rp)

#//////////////////////////////
#RING
#//////////////////////////////
Ri=fi*Rp
Re=fe*Rp
print "Rings derivative:"
print T,"Internal ring (relative) = %.2f Rstar"%(Ri)
print T,"External ring (relative) = %.2f Rstar"%(Re)

#//////////////////////////////
#ORBIT
#//////////////////////////////
Porb=2*pi*sqrt(ap**3/(GCONST*(Mstar+Mp))) #Period (s)
fcen=270*DEG-wp #Central true anomaly (rad)
rcen=ellipseRadiusE(ap,ep,fcen) #r central (km)
Borb=rcen*cos(iorb)/Rstar #Impact parameter (Rstar)
print "Orbit derivative:"
print T,"Period = %e s = %e h = %e d = %e yr"%(Porb,Porb/HOUR,Porb/DAY,Porb/YEAR)
print T,"Central anomaly = %e deg"%(fcen*RAD)
print T,"Central radius = %e km = %e AU = %e Rstar"%(rcen,rcen/AU,rcen/Rstar)
print T,"Impact parameter = %e Rstar"%(Borb)

#//////////////////////////////
#ROTATION MATRICES
#//////////////////////////////
#From orbit to sky
Mi=rotMat([1,0,0],-iorb)
Mw=rotMat([0,0,1],wp)
Mos=dot(Mi,Mw)

fs=linspace(0,2*pi,100)
rs=ellipseRadiusE(ap/Rstar,ep,fs)
xs=rs*cos(fs)
ys=rs*sin(fs)
zs=zeros_like(xs)
rs=array([dot(Mos,AR3(x,y,z)) for x,y,z in zip(xs,ys,zs)])
for i in xrange(len(rs)):
    c='b'
    if rs[i,2]<0:c='r'
    ax.plot([rs[i,0]],[rs[i,1]],'o',markeredgecolor='none',color=c)

#From rings to sky
Mpr=rotMat([0,0,1],phir)
Mir=rotMat([1,0,0],-ir)
Mro=dot(Mpr,Mir)
Mrs=dot(Mi,Mro)

fs=linspace(0,2*pi,100)
rs=ellipseRadiusE(ap/Rstar,ep,fs)
R=100
xs=R*cos(fs)
ys=R*sin(fs)
zs=zeros_like(xs)

rp=dot(Mos,AR3(rcen/Rstar*cos(fcen),rcen/Rstar*sin(fcen),0))
ax.plot([rp[0]],[rp[1]],'go',markersize=10)

rs=array([dot(Mrs,AR3(x,y,z))+AR3(rp[0],rp[1],0) for x,y,z in zip(xs,ys,zs)])
for i in xrange(len(rs)):
    c='b'
    if rs[i,2]<0:c='r'
    ax.plot([rs[i,0]],[rs[i,1]],'o',markeredgecolor='none',color=c,markersize=2)

#########################################
#FIGURES
#########################################
Star=Figure(AR(0,0),1,1,0,'Star')
C=AR(0.0,Borb)
Planet=Figure(C,Rp,Rp,0,'Planet')
Ringe=Figure(C,Re,Re,0,'Ringext')
Ringi=Figure(C,Ri,Ri,0,'Ringint')

#########################################
#SCRIPT
#########################################

#########################################
#PLOT
#########################################
plotEllipse(ax,Star,patch='true',fc='y',ec='none')
plotEllipse(ax,Planet,patch='true',fc='b',ec='none')
plotEllipse(ax,Ringe,color='k')
plotEllipse(ax,Ringi,color='k')

range=1.0*1.5+0.0*ap/Rstar
ax.set_xlim((-range,range))
ax.set_ylim((-range,range))
ax.grid()
fig.savefig("plots/rotations.png")
