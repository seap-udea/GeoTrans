from geotrans import *
from os import system
from sys import argv
fig=plt.figure(figsize=(8,8))
ax=fig.gca()

#########################################
#INPUT PARAMETERS
#########################################
NP=1E3
Rp=0.2
B=0.5

#INCLINATION OF THE RINGS
i=85.0*DEG
#RING TILT
q=240.0*DEG

fe=2.0
fi=1.5
Re=fe*Rp
Ri=fi*Rp

#########################################
#INPUT PARAMETERS
#########################################
Star=Figure(AR(0.0,0.0),1.0,1.0,0.0,'Star')
C=AR(0.0,0.0)
Planet=Figure(C,Rp,Rp,0.0,'Planet')
Ringe=Figure(C,Re,Re*cos(i),q,'Ringe')
Ringi=Figure(C,Ri,Ri*cos(i),q,'Ringi')

#########################################
#INTERSECTION
#########################################
Ar,Ps=ringedPlanetArea(Planet,Ringe,Ringi)

TOL=1E-4

#////////////////////
#CONTACT 1
#////////////////////
x,A,dx,n=contactPosition(Planet,Ringe,Ringi,Ar,B,Side=OUTSIDE,Phase=INGRESS,tola=TOL)
print "C1: x,A,dx,n = ",x,A,dx,n
Planet.C[0]=x
style=dict(color='r',linestyle='-')
plotEllipse(ax,Planet,**style)
plotEllipse(ax,Ringe,**style)
plotEllipse(ax,Ringi,**style)
ax.text(Planet.C[0],Planet.C[1]+Rp,"x=%.2f"%x)

#////////////////////
#CONTACT 2
#////////////////////
x,A,dx,n=contactPosition(Planet,Ringe,Ringi,Ar,B,Side=INSIDE,Phase=INGRESS,tola=TOL)
print "C2: x,A,dx,n = ",x,A,dx,n
Planet.C[0]=x
style=dict(color='r',linestyle='--')
plotEllipse(ax,Planet,**style)
plotEllipse(ax,Ringe,**style)
plotEllipse(ax,Ringi,**style)
ax.text(Planet.C[0],Planet.C[1]+Rp,"x=%.2f"%x)

#////////////////////
#CONTACT 3
#////////////////////
x,A,dx,n=contactPosition(Planet,Ringe,Ringi,Ar,B,Side=INSIDE,Phase=EGRESS,tola=TOL)
print "C3: x,A,dx,n = ",x,A,dx,n
Planet.C[0]=x
style=dict(color='b',linestyle='--')
plotEllipse(ax,Planet,**style)
plotEllipse(ax,Ringe,**style)
plotEllipse(ax,Ringi,**style)
ax.text(Planet.C[0],Planet.C[1]+Rp,"x=%.2f"%x)

#////////////////////
#CONTACT 4
#////////////////////
x,A,dx,n=contactPosition(Planet,Ringe,Ringi,Ar,B,Side=OUTSIDE,Phase=EGRESS,tola=TOL)
print "C4: x,A,dx,n = ",x,A,dx,n
Planet.C[0]=x
style=dict(color='b',linestyle='-')
plotEllipse(ax,Planet,**style)
plotEllipse(ax,Ringe,**style)
plotEllipse(ax,Ringi,**style)
ax.text(Planet.C[0],Planet.C[1]+Rp,"x=%.2f"%x)

#########################################
#PLOT
#########################################
plotEllipse(ax,Star)

RANGE=1.5
ax.set_xlim((-RANGE,RANGE))
ax.set_ylim((-RANGE,RANGE))
ax.grid()
fig.savefig("plots/contacts.png")
