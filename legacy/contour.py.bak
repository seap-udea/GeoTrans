from geotrans import *
from system import *
Ringed=System

AR=lambda x,y:array([x,y])

class Figure(object):
    def __init__(self,C,a,b,ct,st,name):
        self.C=C
        self.a=a
        self.b=b
        self.cost=ct
        self.sint=st
        self.name=name

def plotEllipse(ax,F,patch=False,**args):
    """
    Plot ellipse.  Optional arguments are for "plot" command.
    """
    C=F.C
    a=F.a
    b=F.b
    if patch:
        cir=pat.Circle((F.C[0],F.C[1]),F.a,**args)
        ax.add_patch(cir)
    else:
        Es=linspace(0,2*pi,1000)
        xs=a*cos(Es)
        ys=b*sin(Es)
        rs=array([rotTrans(AR(x,y),F.cost,F.sint,C) for x,y in zip(xs,ys)])
        ax.plot(rs[:,0],rs[:,1],'-',**args)

def plotPlanets(ax,S,Nx=5,Ny=5,
                xmin=0.0,scalex=1.0,ymin=0,scaley=90,
                yoffset=0,
                fh=0,fv=0):
    if Nx>2:
        dc=scalex/(Nx-1)
        cieffs=arange(xmin,xmin+scalex+dc,dc)
    else:
        cieffs=array([xmin])
    if Ny>2:
        dt=scaley/(Ny-1)
        teffs=arange(ymin,ymin+scaley+dt,dt)*DEG
    else:
        teffs=array([ymin])

    if fh==0:
        fh=0.03/S.Rp
    if fv==0:
        fv=fh

    for cieff in cieffs:
        ieff=arccos(cieff)
        for teff in teffs:
            x=(cieff-xmin)/scalex
            y=(teff*RAD-ymin)/scaley+yoffset
            C=AR(x,y)
            Planet=Figure(C,fh*S.Rp,fv*S.Rp,1.0,0.0,'Planet')
            Ringe=Figure(C,fh*S.Re,fv*S.Re*cos(ieff),cos(teff),sin(teff),'Ringext')
            Ringi=Figure(C,fh*S.Ri,fv*S.Ri*cos(ieff),cos(teff),sin(teff),'Ringint')
            plotEllipse(ax,Planet,patch=True,zorder=10,color='k',transform=ax.transAxes)
            plotEllipse(ax,Ringe,zorder=10,color='k',transform=ax.transAxes)
            plotEllipse(ax,Ringi,zorder=10,color='k',transform=ax.transAxes)

def contourPhotoRing():

    S=System
    try:
        qcalc=int(argv[1])
    except:
        qcalc=1

    print BARL,"Calculating Photo-ring effect contours",RBAR

    #////////////////////////////////////////
    #SYSTEM
    #////////////////////////////////////////
    RingedC=copyObject(Ringed)
    ap=Ringed.ap/Ringed.Rstar
    P=Ringed.Porb/HOUR
    ip=Ringed.iorb
    B=RingedC.Borb
    Rp=RingedC.Rp
    rho_true=S.Mstar/(4*pi/3*S.Rstar**3)

    #////////////////////////////////////////
    #MAKE A MAP
    #////////////////////////////////////////
    cieffmin=0.01
    cieffmax=1.0
    Ncieffs=30
    teffmin=0.0*DEG
    teffmax=90.0*DEG
    Nteffs=30
    qcorrected=False

    if qcalc:
        cieffs=linspace(cieffmin,cieffmax,Ncieffs)
        teffs=linspace(teffmin,teffmax,Nteffs)
        #cieffs=[0.5];teffs=[45*DEG]
        
        IS,TS=meshgrid(cieffs,teffs)
        PR=zeros_like(IS)   
        PRN=zeros_like(IS)   
 
        ii=0
        for ci in cieffs:
            i=arccos(ci)
            print "Testing ieff (cos ieff) = ",i*RAD,ci
            jj=0
            for t in teffs:
                print TAB,"Testing teff = ",t*RAD
                RingedC.ieff=i

                RingedC.block=blockFactor(RingedC.tau,i)

                RingedC.Ringext.b=RingedC.Ringext.a*cos(i)
                RingedC.Ringext.cost=cos(t)
                RingedC.Ringext.sint=sin(t)

                RingedC.Ringint.b=RingedC.Ringint.a*cos(i)
                RingedC.Ringint.cost=cos(t)
                RingedC.Ringint.sint=sin(t)
                
                B=RingedC.Borb

                #NUMERICAL
                tcsp=contactTimes(RingedC)
                tT=(tcsp[-1]-tcsp[1])/HOUR
                tF=(tcsp[-2]-tcsp[2])/HOUR

                p=ringedPlanetArea(RingedC)/pi
                pRn=rhoObserved_Seager(p,
                                       RingedC.Rstar,
                                       tT,tF,
                                       P)/rho_true

                #ANALYTICAL
                xpa1=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=-1,sign=-1,
                                     qcorrected=qcorrected)
                xpa2=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=-1,sign=+1,
                                     qcorrected=qcorrected)
                xpa3=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=+1,sign=-1,
                                     qcorrected=qcorrected)
                xpa4=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=+1,sign=+1,
                                     qcorrected=qcorrected)
                tT=P*arcsin((xpa4-xpa1)/(ap*sin(ip)))/(2*pi)
                tF=P*arcsin((xpa3-xpa2)/(ap*sin(ip)))/(2*pi)

                p=analyticalTransitArea(RingedC.Rp,RingedC.block,
                                        RingedC.fi,RingedC.fe,
                                        i)/pi
                pR=rhoObserved_Seager(p,
                                      RingedC.Rstar,
                                      tT,tF,
                                      P)/rho_true

                print 2*TAB,"PR (Numerical) = %.6e, PR (Analytical) = %.6e"%(log10(pRn),log10(pR))

                PR[ii,jj]=log10(pR)
                PRN[ii,jj]=log10(pRn)

                jj+=1
            ii+=1

        savetxt("ISP.dat",IS)
        savetxt("TSP.dat",TS)
        savetxt("PR.dat",PR)
        savetxt("PRN.dat",PRN)

    IS=loadtxt("ISP.dat")
    TS=loadtxt("TSP.dat")
    PR=loadtxt("PR.dat")
    PRN=loadtxt("PRN.dat")
    #data=loadtxt("posterior-PhotoRing-lowtilt.dat") #EXPERIMENT
    
    cmap=plt.get_cmap("rainbow")

    #////////////////////////////////////////
    #CONTOUR
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()

    dmin=PRN.min()
    dmax=PRN.max()
    
    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(PRN),
                  levels=levels,cmap=cmap)

    #ax.plot(cos(data[::5,4]*DEG),abs(data[::5,7]),'k+',zorder=5) #EXP

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Photo-Ring Effect"%(B),position=(0.5,1.02),
                 fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$\log(\rho_{\rm obs}/\rho_\star)$",fontsize=16)
    yts=cbar.ax.get_yticks()
    yl=[]
    levels=[]
    for yt in yts:
        yv=(dmin+yt*(dmax-dmin))
        yl+=["%.3f"%((dmin+yt*(dmax-dmin))*1E0)]
        levels+=[yv]
    cbar.ax.set_yticklabels(yl)
    c=ax.contour(IS,TS*RAD,transpose(PRN),
               levels=levels,
               colors=['k'],linestyles=[':'])
    ax.clabel(c,inline=1,fontsize=10)
    levelsC=levels

    ax.contour(IS,TS*RAD,transpose(PRN),
               levels=[0.0],
               colors=['k'],linestyles=['-'],linewidths=['2'])
    
    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/PhotoRingContour.png")

    #////////////////////////////////////////
    #CONTOUR
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()

    dmin=PR.min()
    dmax=PR.max()
    
    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(PR),
                  levels=levels,cmap=cmap)

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Photo-Ring Effect"%(B),position=(0.5,1.02),
                 fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$\log(\rho_{\rm obs}/\rho_\star)$",fontsize=14)
    yts=cbar.ax.get_yticks()
    yl=[]
    levels=[]
    for yt in yts:
        yv=(dmin+yt*(dmax-dmin))
        yl+=["%.3f"%((dmin+yt*(dmax-dmin))*1E0)]
        levels+=[yv]
    cbar.ax.set_yticklabels(yl)
    c=ax.contour(IS,TS*RAD,transpose(PR),
                 levels=levelsC,
                 colors=['k'],linestyles=[':'])
    ax.clabel(c,inline=1,fontsize=10)

    ax.contour(IS,TS*RAD,transpose(PR),
               levels=[0.0],
               colors=['k'],linestyles=['-'],linewidths=['2'])
    
    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/PhotoRingContour-Analytical.png")

contourPhotoRing()
