"""Generate Photoring contour
"""

from geotrans import *
from system import *
Ringed=System

AR=lambda x,y:array([x,y])

def contourPhotoRing():

    S=System
    try:
        qcalc=int(argv[1])
    except:
        qcalc=1

    print(BARL,"Calculating Photo-ring effect contours",RBAR)

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
        
        IS,TS=meshgrid(cieffs,teffs)
        PRN=zeros_like(IS)   
 
        ii=0
        for ci in cieffs:
            i=arccos(ci)
            print("Testing ieff (cos ieff) = ",i*RAD,ci)
            jj=0
            for t in teffs:
                print(TAB,"Testing teff = ",t*RAD)
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

                PRN[ii,jj]=log10(pRn)

                jj+=1
            ii+=1

        savetxt("tmp/ISP.dat",IS)
        savetxt("tmp/TSP.dat",TS)
        savetxt("tmp/PRN.dat",PRN)

    IS=loadtxt("tmp/ISP.dat")
    TS=loadtxt("tmp/TSP.dat")
    PRN=loadtxt("tmp/PRN.dat")
    
    #////////////////////////////////////////
    #CONTOUR
    #////////////////////////////////////////
    cmap=plt.get_cmap("rainbow")

    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()

    dmin=PRN.min()
    dmax=PRN.max()
    
    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(PRN),
                  levels=levels,cmap=cmap)

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Photo-Ring Effect (numerical)"%(B),position=(0.5,1.02),
                 fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$\log_{10}(\rho_{\rm obs}/\rho_\star)$",fontsize=16)

    levels=linspace(dmin,dmax,10)
    c=ax.contour(IS,TS*RAD,transpose(PRN),
               levels=levels,
               colors=['k'],linestyles=[':'])
    ax.clabel(c,inline=1,fontsize=10)
    levelsC=levels

    ax.contour(IS,TS*RAD,transpose(PRN),
               levels=[0],
               colors=['k','k'],linestyles=['-','-'],linewidths=[2])
    
    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)
    
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/PhotoRingContour.png")
    plt.show()
    
contourPhotoRing()
