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
phir=30.0*DEG #Ring roll angle
tau=1.0 #Opacity
#//////////////////////////////
#ORBIT
#//////////////////////////////
ap=1.0*AU
ep=0.1
iorb=89.8*DEG
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

#APARENT INCLINATION
rax=dot(Mrs,[0.0,0.0,1.0])
ieff=mod(arccos(dot(rax,[0,0,1])),pi/2)
teff=ARCTAN(abs(rax[0]),abs(rax[1]))
print ieff*RAD
print teff*RAD
rlo=dot(Mrs,[0.0,1.0,0.0])
rtr=dot(Mrs,[1.0,0.0,0.0])

#//////////////////////////////
#ORBIT
#//////////////////////////////
Porb=2*pi*sqrt(ap**3/(GCONST*(Mstar+Mp))) #Period (s)
fcen=270*DEG-wp #Central true anomaly (rad)
rcen=ellipseRadiusE(ap,ep,fcen) #r central (km)
Pcen=dot(Mos,AR3(rcen*cos(fcen),rcen*sin(fcen),0))
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
    print "Orbit derivative:"
    print T,"Period = %e s = %e h = %e d = %e yr"%(Porb,Porb/HOUR,Porb/DAY,Porb/YEAR)
    print T,"Central anomaly = %e deg"%(fcen*RAD)
    print T,"Central radius = %e km = %e AU = %e Rstar"%(rcen,rcen/AU,rcen/Rstar)
    print T,"Impact parameter = %e Rstar"%(Borb)
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
Ringe=Figure(C,Re,Re*sin(ieff),teff,'Ringext')
Ringi=Figure(C,Ri,Ri*sin(ieff),teff,'Ringint')

#########################################
#SCRIPT
#########################################

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
rp=3*Rp*rax+Pcen/Rstar
ax.plot([Pcen[0]/Rstar,rp[0]],
        [Pcen[1]/Rstar,rp[1]],'k-',linewidth=3)

rp=3*Rp*rlo+Pcen/Rstar
ax.plot([Pcen[0]/Rstar,rp[0]],
        [Pcen[1]/Rstar,rp[1]],'c-',linewidth=3)

rp=3*Rp*rtr+Pcen/Rstar
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
range=0.0*1.5+0.0*ap/Rstar+5*Rp
ax.set_xlim((-range+xrange,range+xrange))
ax.set_ylim((-range+yrange,range+yrange))
ax.grid()
fig.savefig("plots/transit-real.png")