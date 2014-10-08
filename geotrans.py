###################################################
#MODULES REQUIRED
###################################################
from matplotlib import pyplot as plt
from matplotlib import patches as pat
from cmath import sqrt as csqrt,phase
from numpy import *
from sys import exit

###################################################
#MACROS
###################################################
RAD=180/pi
DEG=pi/180
MAG=lambda P:sqrt(sum(P**2))
AR=lambda x,y:array([x,y])
ARCTAN=lambda num,den:mod(arctan2(num,den),2*pi)
RAND=random.rand

###################################################
#DATA TYPES
###################################################
class Point(object):
    def __init__(self,pos,fig1,fig2):
        self.pos=pos
        self.fig1=fig1
        self.fig2=fig2
    def __str__(self):
        s="("
        s+=str(self.pos)+","
        s+=str(self.pos)+","
        s+=str(self.fig1)+","
        s+=str(self.fig2)+")"
        return s

class Figure(object):
    def __init__(self,C,a,b,name):
        self.C=C
        self.a=a
        self.b=b
        self.name=name
    def __str__(self):
        s="("
        s+="C="+str(self.C)+","
        s+="a="+str(self.a)+","
        s+="b="+str(self.b)+","
        s+="name="+str(self.name)+")"
        return s

###################################################
#CONFIGURATION
###################################################
#VERBOSITY
VERBOSE=0

#BEHAVIOR
PHASETOL=1E-9

#NULL FIGURE
FNULL=Figure(AR(0,0),0,0,'')

###################################################
#ROUTINES
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

def toPoint(v):    
    P=Point(v,FNULL,FNULL)
    return P

#//////////////////////////////
#BASIC
#//////////////////////////////
def figureArea(F):
    return pi*F.a*F.b

def thirdPoly(a,b,c,d,x):
    y=a*x**3+b*x**2+c*x+d
    return y

def fourthPoly(a,b,c,d,e,x):
    y=a*x**4+b*x**3+c*x**2+d*x+e
    return y

def ellipseDiscriminant(F1,F2):
    C=F1.C
    x=C[0]
    y=C[1]
    a=F1.a
    b=F1.b

    v1=1/b**2-1/a**2
    v2=-2*x/a**2
    v3=-2*y/b**2
    v4=x**2/a**2+y**2/b**2-1+1/a**2
    v8=-2*x/a**2
    v10=-2*x/a**2
    #if VERBOSE:print "Auxiliar: ",v1,v2,v3,v4,v8,v10

    e=v2*v10-v4**2
    d=-2*v3*v4
    c=-v2*v8-v3**2-2*v1*v4
    b=-2*v1*v3
    a=-v1**2
    #if VERBOSE:print "Coefficients Fourth: ",a,b,c,d,e

    a0=float(e)/a
    a1=float(d)/a
    a2=float(c)/a
    a3=float(b)/a
    a4=1.0
    #if VERBOSE:print "Auxiliar 2: ",a0,a1,a2,a3
    
    aa=1.0
    bb=-a2
    cc=a1*a3-4*a0
    dd=4*a2*a0-a1**2-a3**2*a0
    #if VERBOSE:print "Coefficients Third: ",aa,bb,cc,dd

    D=18-4*bb**2/(aa*cc)+bb*cc/(aa*dd)-4*cc**2/(bb*dd)-27*aa*dd/(bb*cc)

    #if VERBOSE:print "Original Discriminant = ",D

    if abs(D)<1E-10:D=0
    return a,b,c,d,e,D*aa*bb*cc*dd


def thirdPolyRoot(a,b,c,d):
    """
    REAL ROOTS OF THE THIRD DEGREE POLYNOMIAL
    See wikipedia
    Polynomial of the form:
      a x^3 + b x^2 + c x + d = 0
    """
    #if VERBOSE:print "Third Polynomial Coefficients: ",a,b,c,d

    #DISCRIMINANT
    D=18-4*b**2/(a*c)+b*c/(a*d)-4*c**2/(b*d)-27*a*d/(b*c)
    D*=a*b*c*d
    #if VERBOSE:print "Discriminant: ",D

    D0=b**2-3*a*c
    D1=2*b**3-9*a*b*c+27*a**2*d
    D2=D1**2-4*D0**3
    C0=((D1+csqrt(D2))/2)**(1./3)
    u1=1
    u2=(-1+3**0.5*1j)/2
    u3=(-1-3**0.5*1j)/2
    
    #ROOTS
    xs=[]

    #RETURN ROOTS
    x=-1/(3*a)*(b+u1*C0+D0/(u1*C0))
    #if VERBOSE:print "x1,phase(x1)",x,phase(x)
    if abs(sin(phase(x)))<PHASETOL:xs+=[abs(x.real)]
    else:xs+=[-123]
    x=-1/(3*a)*(b+u2*C0+D0/(u2*C0))
    #if VERBOSE:print "x2,phase(x2)",x,phase(x)
    if abs(sin(phase(x)))<PHASETOL:xs+=[abs(x.real)]
    else:xs+=[-123]
    x=-1/(3*a)*(b+u3*C0+D0/(u3*C0))
    #if VERBOSE:print "x3,phase(x3)",x,phase(x)
    if abs(sin(phase(x)))<PHASETOL:xs+=[abs(x.real)]
    else:xs+=[-123]

    xs=array(xs)
    if len(xs[xs!=-123])>1:D=abs(D)
    #if VERBOSE:print "Effective Discriminant: ",D

    return array(xs),D

def fourthPolyRoots(a,b,c,d,e):
    """
    See: 
    http://mathworld.wolfram.com/QuarticEquation.html

    Polynomial of the form:
      a x^4 + b x^3 + c x^2 + d x + e = 0
    """
    #if VERBOSE:print "Fourth Polynomial Coefficients: ",a,b,c,d,e

    #TRANSFORM TO FORM x^4 + a_3 x^3 + a_2 x^2 + a_1 x + a_0
    a0=float(e)/a
    a1=float(d)/a
    a2=float(c)/a
    a3=float(b)/a
    a4=1.0

    a=1.0
    b=-a2
    c=a1*a3-4*a0
    d=4*a2*a0-a1**2-a3**2*a0
    xs,D=thirdPolyRoot(a,b,c,d)
    #if VERBOSE:print "Solution 3rd Poly: xs,D = ",xs,D

    if D>=0:
        q=-1
        if abs(abs(xs[0])-abs(xs[2]))/abs(abs(xs[0])+abs(xs[2]))<1E-3:q=+1
        return array([123*q,123*q,123*q,123*q]) 

    R=csqrt(0.25*a3**2-a2+xs[xs!=-123])
    D=csqrt(0.75*a3**2-R**2-2*a2+0.25*(4*a3*a2-8*a1-a3**3)/R)
    E=csqrt(0.75*a3**2-R**2-2*a2-0.25*(4*a3*a2-8*a1-a3**3)/R)

    zs=[]
    z=-0.25*a3+0.5*R+0.5*D
    #if VERBOSE:print "z1,phase(z1) = ",z,phase(z)
    if abs(sin(phase(z)))<1E-3:zs+=[abs(real(z))]
    z=-0.25*a3+0.5*R-0.5*D
    #if VERBOSE:print "z2,phase(z2) = ",z,phase(z)
    if abs(sin(phase(z)))<1E-3:zs+=[abs(real(z))]
    z=-0.25*a3-0.5*R+0.5*E
    #if VERBOSE:print "z3,phase(z3) = ",z,phase(z)
    if abs(sin(phase(z)))<1E-3:zs+=[abs(real(z))]
    z=-0.25*a3-0.5*R-0.5*E
    #if VERBOSE:print "z4,phase(z4) = ",z,phase(z)
    if abs(sin(phase(z)))<1E-3:zs+=[abs(real(z))]
    
    #if VERBOSE:print "Solution = ",zs
    return array(zs)

def ellipseRadius(F,cost,sint):
    #Ellipse radius vector
    return F.a*F.b/((F.b*cost)**2+(F.a*sint)**2)**0.5

def ellipsePoint(F,cost,sint):
    #Ellipse point
    r=ellipseRadius(F,cost,sint)
    return AR(r*cost,r*sint)+F.C

def pointInFigure(F,P):
    #If point inside return positive number
    d=P.pos-F.C
    r=MAG(d)
    cost=d[0]/r
    sint=d[1]/r
    rf=ellipseRadius(F,cost,sint)
    ##if VERBOSE:print "cost,rf,r = ",cost,rf,r
    return rf-r

#//////////////////////////////
#INTERSECTION
#//////////////////////////////
def qCircleCircle(F1,F2):
    #1:out,0:intersection,-1:in
    R1=F1.a;R2=F2.a
    D=MAG(F2.C-F1.C)
    if D>=(R1+R2):return 1
    if D<=abs(R1-R2):return -1
    return 0

def cIc(F1,F2):
    #Non-concentric circles R_F1>R_F2
    q=qCircleCircle(F1,F2) or 0
    if q:return \
            Point(AR(123*q,123*q),F1,F2),\
            Point(AR(123*q,123*q),F1,F2)

    R=F1.a
    r=F2.a
    C=F2.C-F1.C
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
    return Point(C+AR(r*cos(t1),r*sin(t1)),F1,F2),Point(C+AR(r*cos(t2),r*sin(t2)),F1,F2)

def cIe(F1,F2):
    #Concentric circle (F2) and ellipse (F1)
    C=F1.C;a=F1.a;b=F1.b;R=F2.a
    if R<b:return \
            Point(AR(-1,-1),F1,F2),\
            Point(AR(-1,-1),F1,F2),\
            Point(AR(-1,-1),F1,F2),\
            Point(AR(-1,-1),F1,F2)

    a2=a**2;b2=b**2;R2=R**2
    x=a*sqrt((R2-b2)/(a2-b2))
    y=b*sqrt((a2-R2)/(a2-b2))
    return\
        Point(AR(x,y)+C,F1,F2),\
        Point(AR(-x,y)+C,F1,F2),\
        Point(AR(-x,-y)+C,F1,F2),\
        Point(AR(x,-y)+C,F1,F2)
    
def eIc(F1,F2):
    #Non-concentric ellipse (F1) and unitary circle (F2)
    C=F1.C
    x=C[0]
    y=C[1]
    a=F1.a
    b=F1.b

    if abs(x)<1E-3:
        a2=a**2;b2=b**2
        c=y;c2=y**2
        D=a2**2+a2*(-b2+c2-1)+b2
        if D<0:
            return \
                Point(AR(-123,-123),F1,F2),\
                Point(AR(-123,-123),F1,F2)
        yi=(a2*c-sqrt(b2*D))/(a2-b2)
        if yi>1:
            return \
                Point(AR(123,123),F1,F2),\
                Point(AR(123,123),F1,F2)
        xi=sqrt(1-yi**2)
        return Point(AR(xi,yi),F1,F2),Point(AR(-xi,yi),F1,F2)

    a,b,c,d,e,D=ellipseDiscriminant(F1,F2)
    #if VERBOSE:print "Ellipse Discriminant: ",D

    if D>=0:
        q=1
        #if VERBOSE:print "No intersection."
        if MAG(C)<1:q=-1
        return \
            Point(AR(123*q,123*q),F1,F2),\
            Point(AR(123*q,123*q),F1,F2)

    if abs(y)<1E-8:
        a=F1.b;b=F1.a
        a2=a**2;b2=b**2
        c=x;c2=x**2
        D=a2**2+a2*(-b2+c2-1)+b2
        if D<0:return Point(AR(-123,-123),F1,F2),Point(AR(-123,-123),F1,F2)
        xi=(a2*c-sqrt(b2*D))/(a2-b2)
        yi=sqrt(1-xi**2)
        return Point(AR(xi,yi),F1,F2),Point(AR(xi,-yi),F1,F2)
    
    ys=fourthPolyRoots(a,b,c,d,e)
    #if VERBOSE:print "Ys = ",ys

    if len(ys[abs(ys)==123])<4:
        x=C[0]
        y=C[1]
        a=F1.a
        b=F1.b
        alpha0=ys**2-1
        alpha1=0*AR(1,1)
        alpha2=1*AR(1,1)
        beta0=ys**2/b**2-2*y*ys/b**2+x**2/a**2+y**2/b**2-1
        beta1=-2*x/a**2*AR(1,1)
        beta2=1/a**2*AR(1,1)
        xs=(alpha2*beta0-alpha0*beta2)/\
            (alpha1*beta2-alpha2*beta1)
    else:
        q=sign(ys[0])
        xs=AR(123*q,123*q)
        ys=AR(123*q,123*q)

    #if VERBOSE:print "Xs = ",xs

    return Point(AR(xs[0],ys[0]),F1,F2),Point(AR(xs[1],ys[1]),F1,F2)

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
    #if VERBOSE:print "Ellipse segment: ",F.name
    #if VERBOSE:print "t1,t2,dt = ",t1*RAD,t2*RAD,dt*RAD
    A=a*b/2*(dt-sin(dt))
    #if VERBOSE:print "A = ",A
    return A

def ellipseSegmentOriented(F,P1,P2,FR):
    C=F.C
    a=F.a
    b=F.b
    #if VERBOSE:print "Oriented Segment"

    #if VERBOSE:print "Figure: ",F
    #if VERBOSE:print "Target Figure: ",FR

    At=figureArea(F)
    #if VERBOSE:print "Figure Area: ",At

    dr1=P1.pos-C
    dr2=P2.pos-C
    t1=ARCTAN(dr1[1]/b,dr1[0]/a)
    t2=ARCTAN(dr2[1]/b,dr2[0]/a)
    dt=abs(t2-t1)

    #if VERBOSE:print "t1,t2,dt = ",t1*RAD,t2*RAD,dt*RAD
    if t1==t2:
        if t1>pi:return -123
        if t1<pi:return +123

    if dt>pi:dt=2*pi-dt
    A=a*b/2*(dt-sin(dt))
    #if VERBOSE:print "Area minimum = ",A

    if t2>270.0*DEG and t1<90.0*DEG:t2=t2-2*pi
    if t1>270.0*DEG and t2<90.0*DEG:t1=t1-2*pi
    tm=mod((t1+t2)/2,2*pi)
    #if VERBOSE:print "tm = ",tm*RAD
    P=ellipsePoint(F,cos(tm),sin(tm))
    cP=pointInFigure(FR,toPoint(P))
    #if VERBOSE:print "Condition = ",cP

    if cP>0:As=A
    else:As=At-A
    #if VERBOSE:print "Area Segment = ",As

    return As

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
    print "Calculating Montecarlo Area with N = %d"%Npoints
    Es=linspace(0,2*pi,100)
    xs=array([])
    ys=array([])
    for F in Fs[1:]:
        C=F.C
        x=C[0];y=C[1]
        a=F.a
        b=F.b
        xs=concatenate((xs,a*cos(Es)+x))
        ys=concatenate((ys,b*sin(Es)+y))
    xmin=xs.min();xmax=xs.max()
    ymin=ys.min();ymax=ys.max()

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
        #if VERBOSE:print "P1,P2 is plane."
        A1=0
    else:
        Fs=commonFigs(P1,P2)
        #if VERBOSE:print "P1,P2 common figures: ",figNames(Fs)
        try:A1=min([ellipseSegment(F,P1,P2) for F in Fs])
        except IndexError:A1=0
        #if VERBOSE:print "A1 = ",A1
    #Side 23
    if shapes[1]==0:
        #if VERBOSE:print "P2,P3 is plane."
        A2=0
    else:
        Fs=commonFigs(P2,P3)
        #if VERBOSE:print "P2,P3 common figures: ",figNames(Fs)
        try:A2=min([ellipseSegment(F,P2,P3) for F in Fs])
        except IndexError:A2=0
        #if VERBOSE:print "A2 = ",A2
    #Side 31
    if shapes[2]==0:
        #if VERBOSE:print "P3,P1 is plane."
        A3=0
    else:
        Fs=commonFigs(P3,P1)
        #if VERBOSE:print "P3,P1 common figures: ",figNames(Fs)
        try:A3=min([ellipseSegment(F,P3,P1) for F in Fs])
        except IndexError:A3=0
        #if VERBOSE:print "A2 = ",A3
        
    Ac=A1+A2+A3
    #if VERBOSE:print "Curved Area = ",Ac
    At=planeTriangle(Ps)
    #if VERBOSE:print "Plane Area = ",At

    return At+Ac
    
def convexQuad(Ps,shapes=[+1,+1,+1,+1]):
    P1,P2,P3,P4=Ps
    #if VERBOSE:print "Convex Quadrangle:"
    #Side 12
    if shapes[0]==0:
        #if VERBOSE:print "P1,P2 is plane."
        A12=0
    else:
        Fs=commonFigs(P1,P2)
        #if VERBOSE:print "P1,P2 common figures: ",figNames(Fs)
        try:A12=min([ellipseSegment(F,P1,P2) for F in Fs])
        except IndexError:A12=0
        #if VERBOSE:print "A12 = ",A12
    #Side 23
    if shapes[1]==0:
        #if VERBOSE:print "P2,P3 is plane."
        A23=0
    else:
        Fs=commonFigs(P2,P3)
        #if VERBOSE:print "P2,P3 common figures: ",figNames(Fs)
        try:A23=min([ellipseSegment(F,P2,P3) for F in Fs])
        except IndexError:A23=0
        #if VERBOSE:print "A23 = ",A23
    #Side 34
    if shapes[2]==0:
        #if VERBOSE:print "P3,P4 is plane."
        A34=0
    else:
        Fs=commonFigs(P3,P4)
        #if VERBOSE:print "P3,P1 common figures: ",figNames(Fs)
        try:A34=min([ellipseSegment(F,P3,P4) for F in Fs])
        except IndexError:A34=0
        #if VERBOSE:print "A34 = ",A34
    #Side 41
    if shapes[3]==0:
        #if VERBOSE:print "P4,P1 is plane."
        A41=0
    else:
        Fs=commonFigs(P4,P1)
        #if VERBOSE:print "P4,P1 common figures: ",figNames(Fs)
        try:A41=min([ellipseSegment(F,P4,P1) for F in Fs])
        except IndexError:A41=0
        #if VERBOSE:print "A41 = ",A41
        
    Ac=A12+A23+A34+A41
    #if VERBOSE:print "Curved Area = ",Ac
    Aq=planeQuad(Ps)
    #if VERBOSE:print "Plane Area = ",Aq
    return Aq+Ac

def convexPolygon(Ps):
    nP=len(Ps)
    #if VERBOSE:print "Sides of the polygon: ",nP
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
        A1=convexQuad((Ps[0],Ps[1],Ps[3],Ps[4]),shapes=[+1,+1,+1,0])
        A2=convexQuad((Ps[0],Ps[4],Ps[5],Ps[2]),shapes=[0,+1,+1,+1])
        A=A1+A2
    return A
        
def leafArea(Ps):
    #FIGURE 2 IS THE SMALLEST FIGURE
    F1=Ps[0].fig1
    F2=Ps[1].fig2

    #F1 IN THE PROGRADE DIRECTION
    A1=ellipseSegmentOriented(F1,Ps[0],Ps[1],F2)
    #if VERBOSE:print "Segment 1 Area: ",A1

    if A1==-123:
        At1=figureArea(F1)
        At2=figureArea(F2)
        #if VERBOSE:print "Areas: ",At1,At2
        if At1>At2:return At2
        else:return At1
    elif A1==123:
        return 0

    A2=ellipseSegmentOriented(F2,Ps[0],Ps[1],F1)
    #if VERBOSE:print "Segment 2 Area: ",A2

    Al=A1+A2
    return Al

#//////////////////////////////
#PLOT
#//////////////////////////////
def plotEllipse(ax,F,**args):
    C=F.C
    a=F.a
    b=F.b
    Es=linspace(0,2*pi,100)
    xs=a*cos(Es)+C[0]
    ys=b*sin(Es)+C[1]
    ax.plot(xs,ys,'-',**args)
    ax.plot([C[0]],[C[1]],'ko')

def plotPoint(ax,P,label='none',**args):
    ax.plot([P.pos[0]],[P.pos[1]],'o',
            markeredgecolor='none',**args)
    if label!='none':
        ax.text(P.pos[0],P.pos[1],label,fontsize=10)
