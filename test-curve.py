from geotrans import *
from os import system
from sys import argv

#########################################
#INPUT PARAMETERS
#########################################
NP=1E3
Rp=0.2
B=0.5

#INCLINATION OF THE RINGS
i=70.0*DEG
#RING TILT
q=20.0*DEG

fe=2.0
fi=1.5
Re=fe*Rp
Ri=fi*Rp

#########################################
#CURVE
#########################################
xmin=-1.4;xmax=1.4;N=50
xss=[]
xts=[]
Ats=[]
Aps=[]
mAs=[]
dAs=[]
Xs=linspace(xmin,xmax,N)
n=0
system("rm -rf confs/*")
for x in Xs:
    print "*"*50
    print "x = %e"%x
    C=AR(x,B)
    Star=Figure(AR(0.0,0.0),1.0,1.0,0.0,'Star')
    Planet=Figure(C,Rp,Rp,0.0,'Planet')
    Ringe=Figure(C,Re,Re*cos(i),q,'Ringe')
    Ringi=Figure(C,Ri,Ri*cos(i),q,'Ringi')
    At,Ps,Feqs=transitFiguresAreaAnalytic(Planet,Ringe,Ringi)
    Ap,Ps,Feqs=transitFiguresAreaAnalytic(Planet,FNULL,FNULL)
    FDIS.write("%e\n"%x)
    print "Analytic Area = %f"%At
    print "Analytic Area (only planet) = %f"%Ap
    xts+=[x];Ats+=[1-At/pi];Aps+=[1-Ap/pi]

    """
    Am,dA,xs,ys=transitFiguresAreaMontecarlo(Planet,Ringe,Ringi,NP=NP)
    xss+=[x];mAs+=[1-Am/pi];dAs+=[dA/pi]
    print "Montecarlo Area = %f +/- %f"%(Am,dA)
    if dA>0:ns=abs(At-Am)/dA
    else:ns=0
    print "Discrepance = %f sigmas"%(ns)
    Am,dA,xsn,ysn=transitFiguresAreaMontecarlo(Feqs[0],Feqs[1],Feqs[2],NP=NP)
    print "Montecarlo Area Normalized = %f +/- %f"%(Am,dA)
    if dA>0:ns=abs(At-Am)/dA
    else:ns=0
    print "Discrepance = %f sigmas"%(ns)
    xss+=[x];mAs+=[1-Am/pi];dAs+=[dA/pi]
    #"""

    #CONFIGURATION PLOT
    fig=plt.figure(figsize=(16,8))
    ax=fig.add_axes([0.05,0.1,0.45,0.8])
    axc=fig.add_axes([0.57,0.1,0.40,0.8])
    plotEllipse(ax,Star,color='k')
    plotEllipse(ax,Planet,color='b',patch=True)
    plotEllipse(ax,Ringe,color='b')
    plotEllipse(ax,Ringi,color='b')
    axc.plot(xts,Ats,'ko-',markersize=5)
    axc.plot(xts,Aps,'bo-',markersize=5)
    axc.set_xlim((xmin,xmax))
    maxd=(Rp**2+Re**2*cos(i)-Ri**2*cos(i))
    axc.set_ylim((1.0-maxd,1.0+maxd/10))
    C=C-0*AR(Re,0)
    RANGE=1.5
    ax.set_xlim((-RANGE*Re+C[0],RANGE*Re+C[0]))
    ax.set_ylim((-RANGE*Re+C[1],RANGE*Re+C[1]))
    ax.text(0.05,0.95,"At = %.2e"%(At),transform=ax.transAxes)
    ax.set_title("x = %.3f"%x)
    suffix="Rp_%e-i_%e-q_%e-b_%e-x_%e"%(Rp,i*RAD,q,B,x)
    fig.savefig("confs/%03d-conf-%s.png"%(n,suffix))

    n+=1
    
#########################################
#PLOT
#########################################
fign=plt.figure(figsize=(8,8))
axn=fign.gca()

axn.plot(Xs,Ats,'ko-',markersize=5)
axn.plot(Xs,Aps,'bo-',markersize=5)
#axn.plot(xss,mAs,'ro')
#axn.errorbar(xss,mAs,yerr=dAs,linestyle='none',color='r')

#########################################
#DECORATION
#########################################
axn.set_xticks(linspace(xmin,xmax,30))
xt=axn.get_xticks()
xl=[]
for x in xt:
    xl+=["%.2f"%x]
axn.set_xticklabels(xl,fontsize=8,rotation=90)

maxd=(Rp**2+Re**2*cos(i)-Ri**2*cos(i))
axc.set_xlim((xmin,xmax))
axc.set_ylim((1.0-maxd,1.0+maxd/10))
axn.grid()
suffix="Rp_%e-i_%e-q_%e"%(Rp,i*RAD,q*RAD)
fign.savefig("store/transit-curve-%s.png"%suffix)

#########################################
#ANIMATION
#########################################
try:
    print "Animating..."
    system("convert -delay %s confs/*.png store/anim-%s.gif"%(argv[1],suffix))
except:
    pass
