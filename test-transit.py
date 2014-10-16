from geotrans import *
fign=plt.figure(figsize=(8,8))
axn=fign.gca()

fige=plt.figure(figsize=(8,8))
axe=fige.gca()

#########################################
#INPUT PARAMETERS
#########################################
NP=1E3
C=AR(-0.829,0.5)
Rp=0.2

#INCLINATION OF THE RINGS
i=85.0*DEG
#RING TILT
q=240.0*DEG

fe=2.0
fi=1.5
Re=fe*Rp
Ri=fi*Rp

#########################################
#FIGURES
#########################################
Star=Figure(AR(0.0,0.0),1.0,1.0,0.0,'Star')
Planet=Figure(C,Rp,Rp,0.0,'Planet')
Ringe=Figure(C,Re,Re*cos(i),q,'Ringe')
Ringi=Figure(C,Ri,Ri*cos(i),q,'Ringi')

#########################################
#AREAS
#########################################
At,Ps,Feqs=transitArea(Planet,Ringe,Ringi)
print "Analytic Area = %f"%At
Am,dA,xs,ys=transitAreaMontecarlo(Planet,Ringe,Ringi,NP=NP)
print "Montecarlo Area = %f +/- %f"%(Am,dA)
if dA>0:ns=abs(At-Am)/dA
else:ns=0
print "Discrepance = %f sigmas"%(ns)
Am,dA,xsn,ysn=transitAreaMontecarlo(Feqs[0],Feqs[1],Feqs[2],NP=NP)

#########################################
#PLOT
#########################################
#PLOT ELLIPSES NORMAL
plotEllipse(axn,Star)
plotEllipse(axn,Planet)
plotEllipse(axn,Ringe)
plotEllipse(axn,Ringi)
#PLOT POINTS

#PLOT ELLIPSES EQUIVALENT
plotEllipse(axe,Star)
plotEllipse(axe,Feqs[0])
plotEllipse(axe,Feqs[1])
plotEllipse(axe,Feqs[2])
for P in Ps:
    if ('re' in P.name) or False:
        plotPoint(axe,P,label=True)

#MONTECARLO POINTS
axn.plot(xs,ys,'r.',markersize=2)
axe.plot(xsn,ysn,'r.',markersize=2)

#########################################
#DECORATION
#########################################
C=C-0*AR(Re,0)
RANGE=1.5
axn.set_xlim((-RANGE*Re+C[0],RANGE*Re+C[0]))
axn.set_ylim((-RANGE*Re+C[1],RANGE*Re+C[1]))
axn.grid()
fign.savefig("plots/transit.png")

C=Feqs[0].C
RANGE=1.5
axe.set_xlim((-RANGE*Re+C[0],RANGE*Re+C[0]))
axe.set_ylim((-RANGE*Re+C[1],RANGE*Re+C[1]))
axe.grid()
fige.savefig("plots/transit-equivalent.png")
