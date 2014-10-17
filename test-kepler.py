from geotrans import *
from physics import *
from constants import *
from numpy import *
fig=plt.figure(figsize=(8,8))
ax=fig.gca()

e=0.2
def EoK(E,M=0.0):
    return E-e*sin(E)-M

Ms=linspace(pi/4,pi,100)
Es=[]
Ers=[]
dEs=[]
for M in Ms:
    #E=eccentricAnomaly(e,M)
    E=eccentricAnomalyFast(e,M)
    Er=newton(EoK,M,args=(M,),tol=1E-15)
    Es+=[E]
    Ers+=[Er]
    dEs+=[(E-Er)/pi]

print "Maximum difference: ",max(abs(array(dEs)))
#ax.plot(Ms,Es,'b-')
#ax.plot(Ms,Ers,'r-')
ax.plot(Ms,dEs,'r-')
fig.savefig("plots/EoK.png")
exit(0)
