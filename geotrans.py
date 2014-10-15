###################################################
#MODULES REQUIRED
###################################################
from matplotlib import pyplot as plt
from matplotlib import patches as pat
from cmath import sqrt as csqrt,phase
from numpy import *
from sys import argv,exit
FIG=plt.figure(figsize=(8,8))
AX=FIG.add_axes([0.1,0.1,0.8,0.8])

###################################################
#MACROS
###################################################
#//////////////////////////////
#MACROS
#//////////////////////////////
RAD=180/pi
DEG=pi/180
MAG=lambda P:sqrt(sum(P**2))
AR=lambda x,y:array([x,y])
RAND=random.rand
ARCTAN=lambda num,den:mod(arctan2(num,den),2*pi)
ZERO=finfo(float).eps
IMAGTOL=1E-5
FIGTOL=1E-3
EQUAL=lambda x,y:abs(x-y)<=ZERO
BARL="*"*50+"\n"
RBAR="\n"+"*"*50
def VERB(routine):print BARL,routine,RBAR
INGRESS=-1
EGRESS=+1
OUTSIDE=+1
INSIDE=-1

#//////////////////////////////
#BEHAVIOR
#//////////////////////////////
VERBOS=0
VERBOSE=[0]*10
VERBOSE[0]=0
VERBOSE[1]=1
VERBOSE[2]=1
VERBOSE[3]=1

#//////////////////////////////
#DEBUGGING
#//////////////////////////////
FDBG=open("roots.txt","w")
FDIS=open("discriminant.txt","w")
DISCRIMINANT=0

###################################################
#DATA TYPES
###################################################
for i in xrange(10):
    if VERBOSE[i]==0:
        VERBOSE[i:]=[0]*(10-i)
        break

class Point(object):
    def __init__(self,pos,fig1,fig2):
        self.pos=pos
        self.fig1=fig1
        self.fig2=fig2
        self.name='P'
    def __str__(self):
        s="("
        s+="Pos = "+str(self.pos)+","
        s+="Fig1 = "+str(self.fig1)+","
        s+="Fig2 = "+str(self.fig2)
        s+=")"
        return s

class Figure(object):
    def __init__(self,C,a,b,t,name):
        self.C=C
        self.a=a
        self.b=b
        self.t=t
        self.name=name
    def __str__(self):
        s="("
        s+="C="+str(self.C)+","
        s+="a="+str(self.a)+","
        s+="b="+str(self.b)+","
        s+="t="+str(self.t)+","
        s+="name="+str(self.name)
        s+=")"
        return s
FNULL=Figure(AR(0,0),0,0,0,'')
FONE=Figure(AR(0,0),1.0,1.0,0,'')

###################################################
#CONFIGURATION
###################################################
#//////////////////////////////
#DATA RELATED
#//////////////////////////////
def commonFigs(P1,P2):
    Fs=list(set([P1.fig1,P1.fig2])&set([P2.fig1,P2.fig2]))
    return Fs

def figNames(Fs):
    s=""
    for F in Fs:
        s+=F.name+","
    s=s.strip(",")
    return s

def pointNames(Ps):
    s=""
    for P in Ps:
        s+=P.name+","
    s=s.strip(",")
    return s

def toPoint(v):    
    P=Point(v,FNULL,FNULL)
    return P

#//////////////////////////////
#BASIC
#//////////////////////////////
def bisectFun(fun,x1,x2,tol=1E-3,maxfun=10,**args):
    f1=fun(x1,**args)
    if EQUAL(f1,0):return x1,1
    f2=fun(x2,**args)
    if EQUAL(f2,0):return x2,1
    if f1*f2>0:return -1,-1
    fc=1
    ifun=2
    while True:
        if VERBOSE[1]:print "x1,f1,x2,f2 = ",x1,f1,x2,f2
        #x=x2-f2*(x2-x1)/(f2-f1) #False position
        x=(x1+x2)/2.0 #Bisection
        if (abs(x2-x1)<=tol) or (ifun>=2+maxfun):break
        fc=fun(x,**args)
        ifun+=1
        if EQUAL(fc,0):return x,ifun
        if f1*fc<0:
            x2=x
            f2=fc
        else:
            f1=fc
            x1=x
    return x,ifun

def rotTrans(r,t,b):
    M=array([[cos(t),-sin(t)],[sin(t),cos(t)]])
    rp=dot(M,r)+b
    return rp

def realArray(zs):
    ims=imag(zs)
    #return real(zs[abs(ims)<=ZERO])
    return real(zs[abs(ims)<=IMAGTOL])

def roundComplex(zs):
    for i in xrange(len(zs)):
        if abs(zs[i].imag)<=IMAGTOL:zs[i]=real(zs[i])
    return zs

def sortPolygonVertices(Ps):
    Np=len(Ps)
    if Np==0:return Ps
    C=AR(0,0)
    for P in Ps:C=C+P.pos
    C=C/Np
    i=0
    qs=[]
    for P in Ps:
        dr=P.pos-C
        q=ARCTAN(dr[1],dr[0])
        #if VERBOSE[0]:print q*RAD
        qs+=[q]
        i+=1
    ies=argsort(qs)
    return Ps[ies]

def figureArea(F):
    return pi*F.a*F.b

def thirdPoly(a,b,c,d,x):
    y=a*x**3+b*x**2+c*x+d
    return y

def fourthPoly(a,b,c,d,e,x):
    y=a*x**4+b*x**3+c*x**2+d*x+e
    return y

def ellipseCoefficients(F):
    if VERBOSE[3]:VERB("ellipseCoefficients")
    C=F.C
    x=C[0];x2=x**2
    y=C[1];y2=y**2
    a=F.a;a2=a**2
    b=F.b;b2=b**2

    v1=1-b2/a2
    v2=-2*x*(b2/a2)
    v3=-2*y
    v4=x2*(b2/a2)+y2-b2+(b2/a2)
    v8=-2*x*(b2/a2)
    v10=-2*x*(b2/a2)
    vs=array([v1,v2,v3,v4,v8,v10])
    if VERBOSE[3]:print "Auxiliar variables: ",vs
    
    e=v2*v10-v4**2
    d=-2*v3*v4
    c=-v2*v8-v3**2-2*v1*v4
    b=-2*v1*v3
    a=-v1**2
    if VERBOSE[3]:print "Coefficients: ",a,b,c,d,e

    return a,b,c,d,e

def ellipseDiscriminant(F):
    a,b,c,d,e=ellipseCoefficients(F)
    a0=float(e)/a
    a1=float(d)/a
    a2=float(c)/a
    a3=float(b)/a
    a4=1.0
    a=1.0
    b=-a2
    c=a1*a3-4*a0
    d=4*a2*a0-a1**2-a3**2*a0
    D=18*a*b*c*d-4*b**3*d+b**2*c**2-4*a*c**3-27*a**2*d**2
    return D

def exc(F):
    Star=Figure(AR(0.0,0.0),1.0,1.0,0.0,'Star')
    #ROOT FUNCTION
    def fun(x):
        if VERBOSE[1]:print "Fun @ x = ",x
        F.C[0]=x
        if VERBOSE[1]:print "Center = ",F.C
        if EQUAL(F.t,0):
            Cg=abs(F.C)
        else:
            Ca=rotTrans(F.C,-F.t,AR(0,0))
            if VERBOSE[1]:print "Rotated Center = ",Ca
            Cg=abs(Ca)
        if VERBOSE[1]:print "Effective Center = ",Cg
        Faux=Figure(Cg,F.a,F.b,0.0,'Auxiliar')
        Ps=eIc(Faux,Star)
        Pposs=array([P.pos for P in Ps])
        Pposs=array([abs(P.pos[0]) for P in Ps])
        l=Pposs[Pposs==123]
        if len(l)==0:return -1
        else:return +1

    #FINDING EXTREME
    C=F.C
    a=F.a
    b=F.b
    t=F.t
    ya1=C[1]+a*sin(t)
    ya2=C[1]-a*sin(t)
    yb1=C[1]+b*cos(t)
    yb2=C[1]-b*cos(t)
    ys=[ya1,ya2,yb1,yb2]
    if VERBOSE[1]:print "Sets of y: ",ys
    ym=min(ys)
    if VERBOSE[1]:print "Minimum: ",ym
    if ym<1:
        xm=sqrt(1-ym**2)
        if VERBOSE[1]:print "Center xm = ",xm
        xm=sqrt(1-ym**2)+a*cos(t)
        if VERBOSE[1]:print "Best guess x: ",xm
    else:
        if VERBOSE[1]:print "No intersection."
        return -1.0

    xmin=max(xm-b,0.0)
    xmax=max(xm+b,0.0)
    if VERBOSE[1]:print "Interval: ",xmin,xmax

    #FINDING ROOT
    x,n=bisectFun(fun,xmin,xmax,tol=1E-3,maxfun=10)
    if VERBOSE[1]:print "Number of evaluations: ",n
    return x

def cxc(F):
    b=F.C[1]
    R=1.0
    r=F.a
    x=sqrt((R+r)**2-b**2)
    return x

def thirdPolyRoot(a,b,c,d):
    #ROOT OF: a x^3 + b x^2 + c x + d = 0
    if VERBOSE[3]:VERB("thirdPolyRoot")
    if VERBOSE[3]:FDBG.write("Third polynomial:%e,%e,%e,%e\n"%(a,b,c,d))

    #AUXILIAR VARIABLES
    D=18*a*b*c*d-4*b**3*d+b**2*c**2-4*a*c**3-27*a**2*d**2
    if VERBOSE[3]:FDBG.write("Classical discriminant:%e\n"%(D))
    if VERBOSE[3] or True:FDIS.write("%e %d "%(D,sign(D)))
    
    D0=b**2-3*a*c
    D1=2*b**3-9*a*b*c+27*a**2*d
    D2=D1**2-4*D0**3
    C0=((D1+csqrt(D2))/2)**(1./3)
    u1=1;u2=(-1+3**0.5*1j)/2;u3=(-1-3**0.5*1j)/2
    
    #ROOTS
    xs=[]

    #RETURN ROOTS
    xs+=[-1/(3*a)*(b+u1*C0+D0/(u1*C0))]
    xs+=[-1/(3*a)*(b+u2*C0+D0/(u2*C0))]
    xs+=[-1/(3*a)*(b+u3*C0+D0/(u3*C0))]

    if VERBOSE[3]:print "Roots poly-3:",xs
    if VERBOSE[3]:FDBG.write("Roots poly-3:"+str(xs)+"\n")
    return array(xs)

def fourthPolyRoots(a,b,c,d,e):
    #ROOTS OF: a x^4 + b x^3 + c x^2 + d x + e = 0
    if VERBOSE[3]:VERB("fourthPolyRoots")
    if VERBOSE[3]:FDBG.write("Fourth polynomial:%e,%e,%e,%e,%e\n"%(a,b,c,d,e))

    #NORMALIZE
    a0=float(e)/a
    a1=float(d)/a
    a2=float(c)/a
    a3=float(b)/a
    a4=1.0
    if VERBOSE[3]:print "Polynomial (%f) x^4 + (%f) x^3 + (%f) x^2 + (%f) x + (%f) = 0"%(a4,a3,a2,a1,a0)

    #RELATED THIRD POLYNOMIAL
    a=1.0
    b=-a2
    c=a1*a3-4*a0
    d=4*a2*a0-a1**2-a3**2*a0
    xs=thirdPolyRoot(a,b,c,d)
    xs=realArray(xs)
    if VERBOSE[3]:print "Real roots poly-3= ",xs
    if VERBOSE[3]:FDBG.write("Real roots poly-3: %s\n"%(str(xs)))

    #CHOOSE REAL ROOTS
    Rs=array([csqrt(0.25*a3**2-a2+x) for x in xs])
    if VERBOSE[3]:FDBG.write("Auxiliar Rs: %s\n"%(str(Rs)))

    if VERBOSE[3]:print "Auxiliar Rs = ",Rs
    Rs=roundComplex(Rs)
    if VERBOSE[3]:print "Rounded Auxiliar Rs = ",Rs
    
    #SOLVE POLYNOMIAL
    Ds=[];Es=[]
    for x,R in zip(xs,Rs):
        if abs(R)<IMAGTOL:
            if VERBOSE[3]:print "Auxiliar R is effectively 0."
            Ds+=[csqrt(0.75*a3**2-2*a2+2*csqrt(x**2-4*a0))]
            Es+=[csqrt(0.75*a3**2-2*a2-2*csqrt(x**2-4*a0))]
        else:
            if VERBOSE[3]:print "Auxiliar R is effectively different from 0."
            Ds+=[csqrt(0.75*a3**2-R**2-2*a2+0.25*(4*a3*a2-8*a1-a3**3)/R)]
            Es+=[csqrt(0.75*a3**2-R**2-2*a2-0.25*(4*a3*a2-8*a1-a3**3)/R)]
    Ds=array(Ds);Es=array(Es)

    z1s=-0.25*a3+0.5*Rs+0.5*Ds
    z2s=-0.25*a3+0.5*Rs-0.5*Ds
    z3s=-0.25*a3-0.5*Rs+0.5*Es
    z4s=-0.25*a3-0.5*Rs-0.5*Es

    zs=array(zip(z1s,z2s,z3s,z4s))
    if VERBOSE[3]:FDBG.write("Roots poly-4:"+str(zs)+"\n")
    if VERBOSE[3]:print "Root sets poly-4: ",zs

    return zs

def ellipseRadius(F,cost,sint):
    #Ellipse radius vector
    return F.a*F.b/((F.b*cost)**2+(F.a*sint)**2)**0.5

def ellipsePoint(F,cost,sint):
    #Ellipse point
    r=ellipseRadius(F,cost,sint)
    return AR(r*cost,r*sint)+F.C

def pointInFigure(F,P):
    #If point inside return positive number
    t=F.t
    d=P.pos-F.C
    d=rotTrans(d,-t,AR(0,0))
    r=MAG(d)
    if EQUAL(r,ZERO):return +F.a
    cost=d[0]/r
    sint=d[1]/r
    rf=ellipseRadius(F,cost,sint)
    #if VERBOSE[0]:print "cost,rf,r = ",cost,rf,r
    return (rf-r)/r

#//////////////////////////////
#INTERSECTION
#//////////////////////////////
def qCircleCircle(F1,F2):
    #1:out,0:intersection,-1:in
    R1=F1.a;R2=F2.a
    D=MAG(F2.C-F1.C)
    if D>=(R1+R2):return -1
    if D<=abs(R1-R2):return 1
    return 0

def cIc(F1,F2):
    #NON-CONCENTRIC CIRCLES (F1<F2)
    q=qCircleCircle(F1,F2) or 0
    if q:return \
            Point(AR(123*q,123*q),F1,F2),\
            Point(AR(123*q,123*q),F1,F2)

    r=F1.a
    C=F1.C-F2.C
    R=F2.a
    x=C[0];y=C[1]
    D=MAG(C)
    r2=r**2;R2=R**2;D2=D**2
    dis=sqrt(\
        -D2**2-2*D2*r2+2*D2*R2-\
             r2**2+2*r2*R2+4*r2*x**2+\
             4*r2*y**2-R2**2)
    det=-2*r*y
    den=D2+r2-2*r*C[0]-R2
    t1=mod(2*arctan2(+dis+det,den),2*pi)
    t2=mod(2*arctan2(-dis+det,den),2*pi)
    return \
        Point(C+AR(r*cos(t1),r*sin(t1)),F1,F2),\
        Point(C+AR(r*cos(t2),r*sin(t2)),F1,F2)

def eIc(F1,F2):
    #NON-CONCENTRIC ELLIPSE AND UNITARY CIRCLE (F1<F2)
    if VERBOSE[3]:VERB("eIc")
    if VERBOSE[3]:print "Figures: ",F1,F2
    if VERBOSE[3]:FDBG.write("eIc Figures: %s, %s\n"%(F1,F2))
    
    C=F1.C
    x=C[0]
    y=C[1]
    a=F1.a
    b=F1.b
    qin=sign(pointInFigure(F2,toPoint(C)))
    if VERBOSE[3]:print "qin = ",qin

    #CHECK IF ELLIPSE IS ACTUALLY A CIRCLE AND 
    if EQUAL(a,b):
        if VERBOSE[3]:print "Ellipse is a Circle."
        return cIc(F1,F2)

    ac,bc,cc,dc,ec=ellipseCoefficients(F1)

    ys=fourthPolyRoots(ac,bc,cc,dc,ec)
    ys=realArray(ys)

    if VERBOSE[3]:FDBG.write("Real ys: %s\n"%(str(ys)))
    if VERBOSE[3]:print "Intersetion in y (Real) = ",ys
    ys=ys[abs(ys)<1]
    if VERBOSE[3]:print "Intersection in y (In range) = ",ys

    if len(ys)==0:
        return \
            Point(AR(123*qin,123*qin),F1,F2),\
            Point(AR(123*qin,123*qin),F1,F2)
    
    #SOLUTION
    qtrad=0
    try:
        xs=[]
        alpha1=0
        alpha2=1
        beta1=-2*x/a**2
        beta2=1/a**2
        det=(alpha1*beta2-alpha2*beta1)
        if VERBOSE[3]:print "Determinant: ",det
        if abs(det)>=1E-13:
            for yp in ys:
                alpha0=yp**2-1
                beta0=yp**2/b**2-2*y*yp/b**2+x**2/a**2+y**2/b**2-1
                xs+=[(alpha2*beta0-alpha0*beta2)/det]
            xs=array(xs)
        else:qtrad=1
    except ValueError as error:
        print error
        qtrad=2

    if qtrad:
        if VERBOSE[3]:print "Using traditional formula for x (reason = %d)"%qtrad
        xs=[]
        for y in ys:
            x=sqrt(1-y**2)
            xs+=[x,-x]
        xs=array(xs)
        
    if VERBOSE[3]:print "Intersection in x (qtrad = %d) = "%qtrad,xs

    qps=[]
    i=0
    for x,y in zip(xs,ys):
        if VERBOSE[3]:print "Testing couple: ",x,y
        qps+=[pointInFigure(F1,toPoint(AR(x,y)))]
        if VERBOSE[3]:print "Value: ",qps[i]
        i+=1
    qps=abs(array(qps))
    iargs=qps.argsort()
    
    return \
        Point(AR(xs[iargs[0]],ys[iargs[0]]),F1,F2),\
        Point(AR(xs[iargs[1]],ys[iargs[1]]),F1,F2)

def cIe(F1,F2):
    #Concentric circle (F2) and ellipse (F1)
    C=F1.C;a=F1.a;b=F1.b;R=F2.a
    if R<b:return \
            Point(AR(-123,-123),F1,F2),\
            Point(AR(-123,-123),F1,F2),\
            Point(AR(-123,-123),F1,F2),\
            Point(AR(-123,-123),F1,F2)

    a2=a**2;b2=b**2;R2=R**2
    x=a*sqrt((R2-b2)/(a2-b2))
    y=b*sqrt((a2-R2)/(a2-b2))
    return\
        Point(AR(x,y)+C,F1,F2),\
        Point(AR(-x,y)+C,F1,F2),\
        Point(AR(-x,-y)+C,F1,F2),\
        Point(AR(x,-y)+C,F1,F2)
    
#//////////////////////////////
#AREAS
#//////////////////////////////
def ellipseSegment(F,P1,P2):
    C=F.C
    a=F.a
    b=F.b
    dr1=P1.pos-C
    dr2=P2.pos-C
    t1=ARCTAN(dr1[1]/b,dr1[0]/a)
    t2=ARCTAN(dr2[1]/b,dr2[0]/a)
    dt=abs(t2-t1)
    if dt>pi:dt=2*pi-dt
    #if VERBOSE[0]:print "Ellipse segment: ",F.name
    #if VERBOSE[0]:print "t1,t2,dt = ",t1*RAD,t2*RAD,dt*RAD
    A=a*b/2*(dt-sin(dt))
    #if VERBOSE[0]:print "A = ",A
    return A

def ellipseSegmentOriented(F,P1,P2,sgn=+1):
    
    if VERBOSE[2]:VERB("ellipseSegmentOriented")
    if VERBOSE[2]:print "Ellipse Segment between points (sgn = %d): "%sgn,P1,P2
    q1=pointInFigure(F,P1)
    q2=pointInFigure(F,P2)
    if VERBOSE[2]:print "Condition point 1:",q1
    if VERBOSE[2]:print "Condition point 2:",q2
    if q1<-FIGTOL or q2<-FIGTOL:
        return 0

    C=F.C
    a=F.a
    b=F.b
    #if VERBOSE[0]:print "Figure: ",F
    #if VERBOSE[0]:print "Sign: ",sgn
    
    #ANGLES
    dr1=P1.pos-C
    dr2=P2.pos-C
    t1=ARCTAN(dr1[1]/b,dr1[0]/a)
    t2=ARCTAN(dr2[1]/b,dr2[0]/a)
    dt=abs(t2-t1)
    if VERBOSE[0]:print "t1,t2,dt = ",t1*RAD,t2*RAD,dt*RAD

    #FAKE SEGMENT
    if t1==t2:
        if pointInFigure(F,P1)<ZERO:
            if VERBOSE[3]:print "Segment closed."
            if sgn<0:
                if VERBOSE[3]:print "Minor curve: All area."
                return pi*a*b
            if sgn>0:
                if VERBOSE[3]:print "Major curve: No area."
                return 0.0
        if t1>pi:return -123
        if t1<pi:return +123

    if dt>pi:dt=2*pi-dt
    A=a*b/2*(dt-sin(dt))
    #if VERBOSE[0]:print "Minimum area= ",A

    if sgn<0:
        #SMALL FIGURE CONDITION
        R=P2.pos-P1.pos
        det=(C[0]*R[1]-C[1]*R[0])
        if EQUAL(det,ZERO):det=1
        lam=(P1.pos[0]*R[1]-P1.pos[1]*R[0])/det
        #if VERBOSE[0]:print "Lambda = ",lam
        if lam>1:A=pi*a*b-A

    #if VERBOSE[0]:print "Area Segment = ",A

    return A

def planeTriangle(Ps):
    #HERON'S FORMULA
    P1,P2,P3=Ps
    S=[MAG(P1.pos-P2.pos),MAG(P2.pos-P3.pos),MAG(P1.pos-P3.pos)]
    s=sum(S)/2
    A=sqrt(s*prod(s-S))
    return A

def planeQuad(Ps):
    #BRETSCHNEIDER'S FORMULA
    P1,P2,P3,P4=Ps
    a=MAG(P1.pos-P2.pos)
    b=MAG(P2.pos-P3.pos)
    c=MAG(P3.pos-P4.pos)
    d=MAG(P4.pos-P1.pos)
    S=array([a,b,c,d])
    p=MAG(P2.pos-P4.pos)
    q=MAG(P1.pos-P3.pos)
    s=sum(S)/2
    Aq=sqrt(prod(s-S)-0.25*(a*c+b*d+p*q)*(a*c+b*d-p*q))
    return Aq

def montecarloArea(Fs,oper,Npoints=1E3):
    #GET MONTECARLO AREA
    #if VERBOSE[0]:print "Calculating Montecarlo Area with N = %d"%Npoints
    Es=linspace(0,2*pi,100)
    cosEs=cos(Es)
    sinEs=sin(Es)
    xs=array([])
    ys=array([])
    for F in Fs[1:]:
        C=F.C
        cost=cos(F.t);sint=sin(F.t)
        x=C[0];y=C[1]
        a=F.a
        b=F.b
        xs=concatenate((xs,(a*cosEs*cost-b*sinEs*sint)+x))
        ys=concatenate((ys,(a*cosEs*sint+b*sinEs*cost)+y))
    xmin=xs.min();xmax=xs.max()
    ymin=ys.min();ymax=ys.max()
    if VERBOSE[0]:print "Montecarlo xmin,xmax = ",xmin,xmax
    if VERBOSE[0]:print "Montecarlo ymin,ymax = ",ymin,ymax

    xs=[];ys=[]
    i=0
    c=0
    while i<Npoints:
        #GENERATE POINT
        P=Point(AR(xmin+(xmax-xmin)*RAND(),ymin+(ymax-ymin)*RAND()),FNULL,FNULL)
        cond=True
        o=0
        for F in Fs:
            cond=(cond and ((oper[o]*pointInFigure(F,P))>=0))
            o+=1

        if cond:
            c+=1
            xs+=[P.pos[0]]
            ys+=[P.pos[1]]
        i+=1
        
    A=(xmax-xmin)*(ymax-ymin)*(float(c)/i)
    if c>0:dA=(1/sqrt(c))*A
    else:dA=0
    return A,dA,xs,ys

def convexTriangle(Ps,shapes=[+1,+1,+1]):
    #shape:+1=curved out,0=plane,-1=curved in
    P1,P2,P3=Ps

    #Side 12
    if shapes[0]==0:
        #if VERBOSE[0]:print "P1,P2 is plane."
        A1=0
    else:
        Fs=commonFigs(P1,P2)
        #if VERBOSE[0]:print "P1,P2 common figures: ",figNames(Fs)
        try:A1=min([ellipseSegment(F,P1,P2) for F in Fs])
        except IndexError:A1=0
        #if VERBOSE[0]:print "A1 = ",A1
    #Side 23
    if shapes[1]==0:
        #if VERBOSE[0]:print "P2,P3 is plane."
        A2=0
    else:
        Fs=commonFigs(P2,P3)
        #if VERBOSE[0]:print "P2,P3 common figures: ",figNames(Fs)
        try:A2=min([ellipseSegment(F,P2,P3) for F in Fs])
        except IndexError:A2=0
        #if VERBOSE[0]:print "A2 = ",A2
    #Side 31
    if shapes[2]==0:
        #if VERBOSE[0]:print "P3,P1 is plane."
        A3=0
    else:
        Fs=commonFigs(P3,P1)
        #if VERBOSE[0]:print "P3,P1 common figures: ",figNames(Fs)
        try:A3=min([ellipseSegment(F,P3,P1) for F in Fs])
        except IndexError:A3=0
        #if VERBOSE[0]:print "A2 = ",A3
        
    Ac=A1+A2+A3
    #if VERBOSE[0]:print "Curved Area = ",Ac
    At=planeTriangle(Ps)
    #if VERBOSE[0]:print "Plane Area = ",At

    return At+Ac
    
def convexQuad(Ps,shapes=[+1,+1,+1,+1]):
    P1,P2,P3,P4=Ps
    #if VERBOSE[0]:print "Convex Quadrangle:"
    #Side 12
    if shapes[0]==0:
        #if VERBOSE[0]:print "P1,P2 is plane."
        A12=0
    else:
        Fs=commonFigs(P1,P2)
        #if VERBOSE[0]:print "P1,P2 common figures: ",figNames(Fs)
        try:A12=min([ellipseSegment(F,P1,P2) for F in Fs])
        except IndexError:A12=0
        #if VERBOSE[0]:print "A12 = ",A12
    #Side 23
    if shapes[1]==0:
        #if VERBOSE[0]:print "P2,P3 is plane."
        A23=0
    else:
        Fs=commonFigs(P2,P3)
        #if VERBOSE[0]:print "P2,P3 common figures: ",figNames(Fs)
        try:A23=min([ellipseSegment(F,P2,P3) for F in Fs])
        except IndexError:A23=0
        #if VERBOSE[0]:print "A23 = ",A23
    #Side 34
    if shapes[2]==0:
        #if VERBOSE[0]:print "P3,P4 is plane."
        A34=0
    else:
        Fs=commonFigs(P3,P4)
        #if VERBOSE[0]:print "P3,P1 common figures: ",figNames(Fs)
        try:A34=min([ellipseSegment(F,P3,P4) for F in Fs])
        except IndexError:A34=0
        #if VERBOSE[0]:print "A34 = ",A34
    #Side 41
    if shapes[3]==0:
        #if VERBOSE[0]:print "P4,P1 is plane."
        A41=0
    else:
        Fs=commonFigs(P4,P1)
        #if VERBOSE[0]:print "P4,P1 common figures: ",figNames(Fs)
        try:A41=min([ellipseSegment(F,P4,P1) for F in Fs])
        except IndexError:A41=0
        #if VERBOSE[0]:print "A41 = ",A41
        
    Ac=A12+A23+A34+A41
    #if VERBOSE[0]:print "Curved Area = ",Ac
    Aq=planeQuad(Ps)
    #if VERBOSE[0]:print "Plane Area = ",Aq
    return Aq+Ac

def convexPolygon(Ps):
    nP=len(Ps)
    #if VERBOSE[0]:print "Sides of the polygon: ",nP
    if nP==0:
        A=0.0
    elif nP==2:
        A=leafArea(Ps)
    elif nP==3:
        A=convexTriangle(Ps)
    elif nP==4:
        A=convexQuad(Ps)
    elif nP==5:
        A1=convexQuad(Ps[:4],shapes=[+1,+1,+1,0])
        A2=convexTriangle((Ps[0],Ps[3],Ps[4]),shapes=[0,1,1])
        A=A1+A2
    elif nP==6:
        A1=convexQuad((Ps[:4]),shapes=[+1,+1,+1,0])
        A2=convexQuad((Ps[0],Ps[3],Ps[4],Ps[5]),shapes=[0,+1,+1,+1])
        A=A1+A2
    return A
        
def leafArea(Ps):
    if VERBOSE[1]:VERB("leafArea")

    #BY DEFAULT FIGURE 1 > FIGURE2
    F1=Ps[0].fig1;AF1=figureArea(F1)
    F2=Ps[1].fig2;AF2=figureArea(F2)
    sgn=+1

    if AF1<AF2:
        sgn=-1
        AF=AF1
    else:AF=AF2

    if VERBOSE[1]:print "Figures: ",F1,F2
    if VERBOSE[1]:print "Larger figure (+1 if F1): ",sgn

    #CHECK IF POINTS ARE SPECIAL
    if Ps[0].pos[0]==123:
        #SMALLER FIGURE INSIDE
        if VERBOSE[1]:print "Small figure completely inside."
        return AF
    elif Ps[0].pos[0]==-123:
        #SMALLER FIGURE OUTSIDE
        if VERBOSE[1]:print "Small figure completely outside."
        return 0

    Pspos=[P.pos for P in Ps]
    if VERBOSE[1]:print "Points: ",pointNames(Ps)
    if VERBOSE[1]:print "Point positions: ",Pspos

    A1=ellipseSegmentOriented(F1,Ps[0],Ps[1],sgn=sgn)
    A2=ellipseSegmentOriented(F2,Ps[0],Ps[1],sgn=-sgn)
    if A1==0 and A2==0:return 0
    if VERBOSE[1]:print "Segment 1 Area: ",A1
    if VERBOSE[1]:print "Segment 2 Area: ",A2

    if sgn>0:
        if A2==-123:return figureArea(F2)
        if A2==123:return 0
    else:
        if A1==-123:return figureArea(F1)
        if A1==123:return 0

    Al=A1+A2
    return Al

#//////////////////////////////
#TRANSIT AREA
#//////////////////////////////
def transitAreaAnalytic(C,Rp,Re,Ri,i):
    #FIGURES
    Cg=abs(C)
    Star=Figure(AR(0.0,0.0),1.0,1.0,'Star')
    Planet=Figure(Cg,Rp,Rp,'Planet')
    Ringe=Figure(Cg,Re,Re*cos(i),'Ringe')
    Ringi=Figure(Cg,Ri,Ri*cos(i),'Ringi')

    #INTERSECTION POINTS
    Psa=[]

    #STAR AND PLANET
    Psp1,Psp2=cIc(Planet,Star)
    Psp1.name='Psp1';Psp2.name='Psp2';
    Psa+=[Psp1,Psp2]
    Asp=leafArea((Psp1,Psp2))
    if VERBOSE[0]:print "Area planet inside star: ",Asp

    #IF NO RINGS (i=90) USE ONLY PLANET
    if abs(Ringe.b)<ZERO:
        if VERBOSE[0]:print "Using only planet."
        return Asp,Psa,Feqs

    #STAR AND EXTERNAL RINGS
    Psre1,Psre2=eIc(Ringe,Star)
    Psre1.name='Psre1';Psre2.name='Psre2';
    Psa+=[Psre1,Psre2]
    Psn=[str(P) for P in Psa]
    if VERBOSE[0]:print "Points Star-External Ring: ",Psre1.pos,Psre2.pos
    
    Psri1,Psri2=eIc(Ringi,Star)
    Psri1.name='Psri1';Psri2.name='Psri2';
    Psa+=[Psri1,Psri2]
    Psn=[str(P) for P in Psa]
    if VERBOSE[0]:print "Points Star-Internal Ring: ",Psri1.pos,Psri2.pos

    Ppre1,Ppre2,Ppre3,Ppre4=cIe(Ringe,Planet)
    Ppre1.name='Ppre1';Ppre2.name='Ppre2';
    Ppre3.name='Ppre3';Ppre4.name='Ppre4';
    Psa+=[Ppre1,Ppre2,Ppre3,Ppre4]
    Ppri1,Ppri2,Ppri3,Ppri4=cIe(Ringi,Planet)
    Ppri1.name='Ppri1';Ppri2.name='Ppri2';
    Ppri3.name='Ppri3';Ppri4.name='Ppri4';
    Psa+=[Ppri1,Ppri2,Ppri3,Ppri4]

    #RING AREAS
    Asre=leafArea((Psre1,Psre2))
    if VERBOSE[0]:print "Area external ring inside star: ",Asre
    Asri=leafArea((Psri1,Psri2))
    if VERBOSE[0]:print "Area internal ring inside star: ",Asri

    #COMMON POINTS
    #EXTERNAL RING
    Ps=array([Ppre1,Psre1,Ppre2,Psp2,Ppre3,Psre2,Ppre4,Psp1])
    Fs=[Star,Planet,Star,Ringe,Star,Planet,Star,Ringe]
    Qs=array([pointInFigure(F,P) for F,P in zip(Fs,Ps)])
    Pine=sortPolygonVertices(Ps[Qs>0])
    if VERBOSE[0]:print "External Ring Exclusion Points: ",pointNames(Pine)

    #INTERNAL RING
    Ps=array([Ppri1,Psri1,Ppri2,Psp2,Ppri3,Psri2,Ppri4,Psp1])
    Fs=[Star,Planet,Star,Ringi,Star,Planet,Star,Ringi]
    Qs=array([pointInFigure(F,P) for F,P in zip(Fs,Ps)])
    Pini=sortPolygonVertices(Ps[Qs>0])
    if VERBOSE[0]:print "Internal Ring Exclusion Points: ",pointNames(Pini)

    #COMMON AREAS
    if len(Pine)==0 or len(Pini)==0:
        Pcusp=toPoint(AR(Planet.C[0],Planet.C[1]+Planet.b))
        qcusp=pointInFigure(Star,Pcusp)
        if VERBOSE[0]:print "Pcusp = ",Pcusp.pos
        if VERBOSE[0]:print "Condition = ",qcusp

    if len(Pine)==0:
        if qcusp<0:Asrec=0.0
        else:Asrec=Asp
    else:Asrec=convexPolygon(Pine)
    if len(Pini)==0:
        if qcusp<0:Asric=0.0
        else:Asric=Asp
    else:Asric=convexPolygon(Pini)

    if VERBOSE[0]:print "External ring common area: ",Asrec
    if VERBOSE[0]:print "Internal ring common area: ",Asric

    #TRANSIT AREA
    Atrans=Asp+Asre-Asri-Asrec+Asric

    return Atrans,Psa

def transitFiguresAreaAnalytic(Planet,Ringe,Ringi):

    #BASIC PROPERTIES
    C=Planet.C
    if VERBOSE[1]:print "Original Center: ",C
    t=Ringe.t
    Rp=Planet.a
    Rea=Ringe.a
    Reb=Ringe.b
    Ria=Ringi.a
    Rib=Ringi.b

    #PUTTING RINGS HORIZONTAL
    Ca=rotTrans(C,-t,AR(0,0))
    if VERBOSE[1]:print "Center for rings horizontal: ",Ca

    #NORMALIZING POSITION
    Cg=abs(Ca)
    if VERBOSE[1]:print "Equivalent center of figures:",Cg

    #FIGURES
    Star=Figure(AR(0.0,0.0),1.0,1.0,0.0,'Star')
    Planet=Figure(Cg,Rp,Rp,0.0,'Planet')
    Ringe=Figure(Cg,Rea,Reb,0.0,'Ringe')
    Ringi=Figure(Cg,Ria,Rib,0.0,'Ringi')
    Feqs=[Planet,Ringe,Ringi]

    #INTERSECTION POINTS
    Psa=[]

    #STAR AND PLANET
    Psp1,Psp2=cIc(Planet,Star)
    Psp1.name='Psp1';Psp2.name='Psp2';
    Psa+=[Psp1,Psp2]
    Asp=leafArea((Psp1,Psp2))
    if VERBOSE[0]:print "Area planet inside star: ",Asp

    #IF NO RINGS (i=90) USE ONLY PLANET
    if abs(Ringe.b)<ZERO:
        if VERBOSE[0]:print "Using only planet."
        return Asp,Psa,Feqs

    #STAR AND EXTERNAL RINGS
    Psre1,Psre2=eIc(Ringe,Star)
    Psre1.name='Psre1';Psre2.name='Psre2';
    Psa+=[Psre1,Psre2]
    Psn=[str(P) for P in Psa]
    if VERBOSE[0]:print "Points Star-External Ring: ",Psre1.pos,Psre2.pos
    #exit(0)
    
    Psri1,Psri2=eIc(Ringi,Star)
    Psri1.name='Psri1';Psri2.name='Psri2';
    Psa+=[Psri1,Psri2]
    Psn=[str(P) for P in Psa]
    if VERBOSE[0]:print "Points Star-Internal Ring: ",Psri1.pos,Psri2.pos
    #exit(0)

    Ppre1,Ppre2,Ppre3,Ppre4=cIe(Ringe,Planet)
    Ppre1.name='Ppre1';Ppre2.name='Ppre2';
    Ppre3.name='Ppre3';Ppre4.name='Ppre4';
    Psa+=[Ppre1,Ppre2,Ppre3,Ppre4]
    Ppri1,Ppri2,Ppri3,Ppri4=cIe(Ringi,Planet)
    Ppri1.name='Ppri1';Ppri2.name='Ppri2';
    Ppri3.name='Ppri3';Ppri4.name='Ppri4';
    Psa+=[Ppri1,Ppri2,Ppri3,Ppri4]

    #RING AREAS
    Asre=leafArea((Psre1,Psre2))
    if VERBOSE[0]:print "Area external ring inside star: ",Asre
    Asri=leafArea((Psri1,Psri2))
    if VERBOSE[0]:print "Area internal ring inside star: ",Asri

    #COMMON POINTS
    #EXTERNAL RING
    Ps=array([Ppre1,Psre1,Ppre2,Psp2,Ppre3,Psre2,Ppre4,Psp1])
    Fs=[Star,Planet,Star,Ringe,Star,Planet,Star,Ringe]
    Qs=array([pointInFigure(F,P) for F,P in zip(Fs,Ps)])
    Pine=sortPolygonVertices(Ps[Qs>0])
    if VERBOSE[0]:print "External Ring Exclusion Points: ",pointNames(Pine)

    #INTERNAL RING
    Ps=array([Ppri1,Psri1,Ppri2,Psp2,Ppri3,Psri2,Ppri4,Psp1])
    Fs=[Star,Planet,Star,Ringi,Star,Planet,Star,Ringi]
    Qs=array([pointInFigure(F,P) for F,P in zip(Fs,Ps)])
    Pini=sortPolygonVertices(Ps[Qs>0])
    if VERBOSE[0]:print "Internal Ring Exclusion Points: ",pointNames(Pini)

    #COMMON AREAS
    if len(Pine)==0 or len(Pini)==0:
        Pcusp=toPoint(AR(Planet.C[0],Planet.C[1]+Planet.b))
        qcusp=pointInFigure(Star,Pcusp)
        if VERBOSE[0]:print "Pcusp = ",Pcusp.pos
        if VERBOSE[0]:print "Condition = ",qcusp

    if len(Pine)==0:
        if qcusp<0:Asrec=0.0
        else:Asrec=Asp
    else:Asrec=convexPolygon(Pine)
    if len(Pini)==0:
        if qcusp<0:Asric=0.0
        else:Asric=Asp
    else:Asric=convexPolygon(Pini)

    if VERBOSE[0]:print "External ring common area: ",Asrec
    if VERBOSE[0]:print "Internal ring common area: ",Asric

    #TRANSIT AREA
    Atrans=Asp+Asre-Asri-Asrec+Asric

    return Atrans,Psa,Feqs

def ringedPlanetBox(Planet,Ringe):
    C=Planet.C
    B=C[1]
    Rp=Planet.a
    a=Ringe.a
    b=Ringe.b
    t=Ringe.t
    qs=linspace(0,2*pi,100)
    E=Figure(AR(0.0,0.0),a,b,0,'Dummy')
    Ps=[rotTrans(ellipsePoint(E,cos(q),sin(q)),t,C) for q in qs]
    ys=[P[1] for P in Ps]
    yemin=min(ys)
    yemax=max(ys)
    ymin=min(yemin,B-Rp)
    ymax=max(yemax,B+Rp)
    return ymin,ymax

def contactPosition(Planet,Ringe,Ringi,Ar,B,Phase=EGRESS,Side=OUTSIDE,tola=1E-3,tolx=0.0,maxfun=10):
    """
    Ar: Area ringed planet
    B: Impact parameter (|B|<~1+Rp+Re)
    Side: Side of the transit, +1:OUTSIDE, -1:INSIDE
    Phase: Phase of the transit, +1:EGRESS, -1:INGRESS
    tola: Tolerance in fraction of ringed planet area (Ar)
    tolx: Tolerance in fraction of star radius
    maxfun: maximum number of callings to Area

    Contacts are:
    T1: INGRESS,OUTSIDE
    T2: INGRESS,INSIDE
    T3: EGRESS,INSIDE
    T4: EGRESS,OUTSIDE
    """
    #COUNTERS
    ifun=0
    
    #INPUT PARAMETERS
    Rp=Planet.a
    a=Ringe.a
    b=Ringe.b
    t=Ringe.t

    #IMPACT PARAMETER
    Planet.C[1]=B
    
    #CHECK EXTREMES IF BOXED ANALYSIS IS REQUIRED
    if abs(B)>(1-Rp-a):
        ymin,ymax=ringedPlanetBox(Planet,Ringe)
        ys=[abs(ymin),abs(ymax)]
        ylow=min(ys);yup=max(ys)
        if abs(ylow)>1:return -1,-1,-1,-1
        if abs(yup)>1:
            if Side==INSIDE:return 0,0,0,0

    #COMPARISON AREA
    if Side>0:Ac=0
    else:Ac=Ar

    #CENTER OF INGRESS/EGRESS
    if abs(B)<1:x1=Phase*sqrt(1-B**2)
    else:x1=0
    Planet.C[0]=x1
    A1,Ps,Fs=transitFiguresAreaAnalytic(Planet,Ringe,Ringi);A1-=Ac
    ifun+=1
    if VERBOS:print "x1,A1=",x1,A1

    #ADVANCE
    x2=x1+Side*Rp/2
    Planet.C[0]=x2
    A2,Ps,Fs=transitFiguresAreaAnalytic(Planet,Ringe,Ringi);A2-=Ac
    ifun+=1
    if VERBOS:print "x2,A2=",x2,A2

    if VERBOS:print "Direction of Transit: A2 - A1 = ",A2-A1

    #FIRST EXTRAPOLATED POSITION
    while True:
        if VERBOS:print "x1,x2 = ",x1,x2
        x=x1-(x2-x1)/(A2-A1)*A1
        Planet.C[0]=x
        A,Ps,Fs=transitFiguresAreaAnalytic(Planet,Ringe,Ringi);A-=Ac
        ifun+=1
        if VERBOS:print "x,A = ",x,A
        
        if A==0:
            if VERBOS:print "Zero area."
            while A2*A==0:
                if VERBOS:print "Searching, x2,A2,x,A = ",x2,A2,x,A
                x=(x2+x)/2
                Planet.C[0]=x
                A,Ps,Fs=transitFiguresAreaAnalytic(Planet,Ringe,Ringi);A-=Ac
                if VERBOS:print "x,A = ",x,A
                ifun+=1
                if VERBOS:raw_input()
            if VERBOS:print "After search, x2,A2,x,A = ",x2,A2,x,A

        x1=x2
        A1=A2
        x2=x
        A2=A
        dx=abs(x2-x1)/2
        if abs(A)<tola*abs(Ar) or dx<tolx or ifun>maxfun:break
        if VERBOS:raw_input()

    if VERBOS:print "x,A=",x,A
    return x,abs(A),dx,ifun

def ringedPlanetArea(Planet,Ringe,Ringi):
    C=Planet.C
    t=Ringe.t
    Rp=Planet.a
    Rea=Ringe.a
    Reb=Ringe.b
    Ria=Ringi.a
    Rib=Ringi.b

    if EQUAL(Reb,ZERO):return pi*Rp**2,[]
    
    #PUTTING RINGS HORIZONTAL
    Ca=rotTrans(C,-t,AR(0,0))

    #NORMALIZING POSITION
    Cg=abs(Ca)

    #FIGURES
    Star=Figure(AR(0.0,0.0),1.0,1.0,0.0,'Star')
    Planet=Figure(Cg,Rp,Rp,0.0,'Planet')
    Ringe=Figure(Cg,Rea,Reb,0.0,'Ringe')
    Ringi=Figure(Cg,Ria,Rib,0.0,'Ringi')
    Feqs=[Planet,Ringe,Ringi]

    #AREAS
    Asp=pi*Rp**2
    Asre=pi*Rea*Reb
    Asri=pi*Ria*Rib
    if VERBOSE[0]:print "Area planet: ",Asp
    if VERBOSE[0]:print "Area external ring: ",Asre
    if VERBOSE[0]:print "Area internal ring: ",Asri

    #INTERSECTION POINTS
    Psa=[]
    Ppre1,Ppre2,Ppre3,Ppre4=cIe(Ringe,Planet)
    Ppre1.name='Ppre1';Ppre2.name='Ppre2';
    Ppre3.name='Ppre3';Ppre4.name='Ppre4';
    Psa+=[Ppre1,Ppre2,Ppre3,Ppre4]
    Ppri1,Ppri2,Ppri3,Ppri4=cIe(Ringi,Planet)
    Ppri1.name='Ppri1';Ppri2.name='Ppri2';
    Ppri3.name='Ppri3';Ppri4.name='Ppri4';
    Psa+=[Ppri1,Ppri2,Ppri3,Ppri4]

    #EXTERNAL RING
    Ps=array([Ppre1,Ppre2,Ppre3,Ppre4])
    Fs=[Planet,Planet,Planet,Planet]
    Qs=array([pointInFigure(F,P) for F,P in zip(Fs,Ps)])
    Pine=sortPolygonVertices(Ps[Qs>-FIGTOL])
    if VERBOSE[0]:print "External Ring Exclusion Points: ",pointNames(Pine)

    #INTERNAL RING
    Ps=array([Ppri1,Ppri2,Ppri3,Ppri4]) 
    Fs=[Planet,Planet,Planet,Planet]
    Qs=array([pointInFigure(F,P) for F,P in zip(Fs,Ps)])
    Pini=sortPolygonVertices(Ps[Qs>-FIGTOL])
    if VERBOSE[0]:print "Internal Ring Exclusion Points: ",pointNames(Pini)

    Pcusp=toPoint(AR(Planet.C[0],Planet.C[1]+Planet.b))
    if VERBOSE[0]:print "Pcusp = ",Pcusp.pos
    if len(Pine)==0:
        qcusp=pointInFigure(Ringe,Pcusp)
        if VERBOSE[0]:print "Condition Ringe= ",qcusp
        if qcusp<-FIGTOL:Asrec=0.0
        else:Asrec=Asp
    else:Asrec=convexPolygon(Pine)
    if len(Pini)==0:
        qcusp=pointInFigure(Ringi,Pcusp)
        if VERBOSE[0]:print "Condition Ringi= ",qcusp
        if qcusp<-FIGTOL:Asric=0.0
        else:Asric=Asp
    else:Asric=convexPolygon(Pini)

    if VERBOSE[0]:print "External ring common area: ",Asrec
    if VERBOSE[0]:print "Internal ring common area: ",Asric

    #RINGED PLANET AREA
    Aringed=Asp+Asre-Asri-Asrec+Asric

    for P in Psa:
        P.pos=rotTrans(P.pos,t,AR(0.0,0.0))

    return Aringed,Psa

def transitAreaMontecarlo(C,Rp,Re,Ri,i,NP=1E3):
    Star=Figure(AR(0.0,0.0),1.0,1.0,'Star')
    Planet=Figure(C,Rp,Rp,'Planet')
    Ringe=Figure(C,Re,Re*cos(i),'Ringe')
    Ringi=Figure(C,Ri,Ri*cos(i),'Ringi')

    mA1,dA1,xs1,ys1=montecarloArea([Star,Ringe,Planet,Ringi],
                                   [+1,+1,-1,-1],Npoints=NP)
    mA2,dA2,xs2,ys2=montecarloArea([Star,Planet],
                                   [+1,+1],Npoints=NP)
    mA=mA1+mA2
    dA=sqrt(dA1**2+dA2**2)
    if dA>0:
        nd=abs(log10(dA))+1
        fmA="%%.%df"%nd
        fmE="%%.%df"%(nd+1)
        mA=float(fmA%mA)
        dA=float(fmE%dA)
    return mA,dA,concatenate((xs1,xs2)),concatenate((ys1,ys2))

def transitFiguresAreaMontecarlo(Planet,Ringe,Ringi,NP=1E3):
    Star=Figure(AR(0.0,0.0),1.0,1.0,0.0,'Star')

    mA1,dA1,xs1,ys1=montecarloArea([Star,Ringe,Planet,Ringi],
                                   [+1,+1,-1,-1],Npoints=NP)
    mA2,dA2,xs2,ys2=montecarloArea([Star,Planet],
                                   [+1,+1],Npoints=NP)
    mA=mA1+mA2
    dA=sqrt(dA1**2+dA2**2)
    if dA>0:
        nd=abs(log10(dA))+1
        fmA="%%.%df"%nd
        fmE="%%.%df"%(nd+1)
        mA=float(fmA%mA)
        dA=float(fmE%dA)
    return mA,dA,concatenate((xs1,xs2)),concatenate((ys1,ys2))

#//////////////////////////////
#PLOT
#//////////////////////////////
def plotEllipse(ax,F,patch=False,**args):
    C=F.C
    a=F.a
    b=F.b
    t=F.t
    if patch:
        cir=pat.Circle((F.C[0],F.C[1]),F.a)
        ax.add_patch(cir)
    else:
        Es=linspace(0,2*pi,1000)
        xs=a*cos(Es)
        ys=b*sin(Es)
        rs=array([rotTrans(AR(x,y),t,C) for x,y in zip(xs,ys)])
        ax.plot(rs[:,0],rs[:,1],'-',**args)
        ax.plot([C[0]],[C[1]],'ko')
        
def plotPoint(ax,P,label=False,**args):
    ax.plot([P.pos[0]],[P.pos[1]],'o',
            markeredgecolor='none',**args)
    if label:
        x=P.pos[0]+0.005
        y=P.pos[1]+0.005
        ax.text(x,y,P.name,fontsize=10)

#//////////////////////////////
#DEBUGGING
#//////////////////////////////
