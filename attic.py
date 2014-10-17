fs=linspace(0,2*pi,100)
rs=ellipseRadiusE(ap/Rstar,ep,fs)
xs=rs*cos(fs)
ys=rs*sin(fs)
zs=zeros_like(xs)
rs=array([dot(Morb,AR3(x,y,z)) for x,y,z in zip(xs,ys,zs)])
ax.plot(rs[:,0],rs[:,1])


e=0.001

from scipy.optimize import newton
def EoK(E,M=0.0):
    return E-e*sin(E)-M

Ms=linspace(0,pi,100)
Es=[]
Ers=[]
dEs=[]
for M in Ms:
    E=eccentricAnomaly2(e,M)
    Er=newton(EoK,M,args=(M,),tol=1E-9)
    Es+=[E]
    Ers+=[Er]
    dEs+=[(E-Er)/pi]

print "Maximum difference: ",max(abs(array(dEs)))
#ax.plot(Ms,Es,'b-')
#ax.plot(Ms,Ers,'r-')
ax.plot(Ms,dEs,'r-')
fig.savefig("plots/EoK.png")
exit(0)

def eccentricAnomaly3(e,M):
    """
    Using a method based on: Mikhola
    """
    b=0.5*M/(4*e+0.5)
    a=(1-e)/(4*e+0.5)
    z=(b+(b*b+a*a*a)**0.5)**(1./3)
    s=z-a/z
    E=M+e*(3*s-4*s**3)
    return E

def eccentricAnomaly(e,M):
    """
    Using a method based on: Mikhola
    """
    p=(1-e)/(4*e+0.5)
    q=0.5*M/(4*e+0.5)
    z=((p*p*p+q*q)**0.5+q)**(1./3)
    w=((p*p*p+q*q)**0.5-q)**(1./3)
    s=z-w
    s=2*q/(z*z+p+p*p/z)

    E=M+e*(3*s-4*s**3)
    return E

print "fcen = ",fcen*RAD
print "tcen = ",tcen/Porb
print "df = ",df*RAD
print "f15 = ",f15*RAD
print "f35 = ",f35*RAD
print "t15 = ",t15/Porb
print "t35 = ",t35/Porb
print "dt = ",dtrans/Porb

#PLOT POINTS
plotPoint(ax,toPoint(AR(Pcen[0],Pcen[1])))
plotPoint(ax,toPoint(AR(P15[0],P15[1])))
plotPoint(ax,toPoint(AR(P35[0],P35[1])))

