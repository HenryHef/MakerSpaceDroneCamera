import numpy as np
#from shapely.geometry.polygon import LinearRing

#centered at 0,0
def arccot(x):
    if(x==0):return np.pi/2.0
    return np.arctan(1/x)

def ABCD_elipse_to_ThetaAB_elpise(A,B,C,D):
    A=A*1.0
    B=B*1.0
    C=C*1.0
    D=D*1.0
    A=A
    B=B
    C=C
    theta=arccot((A-C)/B)/2
    M=np.cos(theta)**2
    N=np.sin(theta)**2
    T=np.tan(theta)
    T2=T**2
    U=np.tan(theta+np.pi/2)
    U2=U**2
    axisA=np.sqrt((D*(1+T2))/(A+B*T+C*T2))
    axisB=np.sqrt((D*(1+U2))/(A+B*U+C*U2))
    return (axisA,axisB,theta)

def ellipse_polyline(elipse_with_a_b_angle, n=100):
    a=elipse_with_a_b_angle[0]
    b=elipse_with_a_b_angle[1]
    angle=elipse_with_a_b_angle[2]
    ax=a*np.cos(angle)
    ay=a*np.sin(angle)
    bx=-b*np.sin(angle)
    by=b*np.cos(angle)

    return [(ax*np.cos(i*1.0/n*np.pi*2)+bx*np.sin(i*1.0/n*np.pi*2),ay*np.cos(i*1.0/n*np.pi*2)+by*np.sin(i*1.0/n*np.pi*2)) for i in range(n)]

def intersections(a, b):
    intersections = []
    for i in range(len(a)):
        for j in range(len(b)):
            if intersectsQ(a[i],a[(i+1)%len(a)],b[j],b[(j+1)%len(b)]):
                intersections.append(intersectionPoint(a[i],a[(i+1)%len(a)],b[j],b[(j+1)%len(b)]))
    return intersections


#!/usr/bin/python

def ccw(A,B,C):
    return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])

def intersectsQ(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)
def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersectionPoint(p1, p2, p3, p4):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    L1 = (A, B, -C)
    A = (p3[1] - p4[1])
    B = (p4[0] - p3[0])
    C = (p3[0]*p4[1] - p4[0]*p3[1])
    L2 = (A, B, -C)
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False

# B is the bisector of AC
def getAO_BO_CO_fromAngleBOA_BOC_AC(angleBOA, angleBOC, AC):
    d=AC/2
    alpha=angleBOC
    beta=angleBOA
    C_A=np.cos(alpha)
    C_B=np.cos(beta)
    C_AB=np.cos(alpha+beta)
    #we find using math that the intersection of 2 certen elipses is the solution to the problem
    elipse1 = ABCD_elipse_to_ThetaAB_elpise(1,-2*C_AB,1,4*(d**2))
    elipse2 = ABCD_elipse_to_ThetaAB_elpise(C_AB**2,2*C_AB-4*C_AB*(C_A**2),1,4*(d**2)*(C_A**2))

    points=intersections(ellipse_polyline(elipse1),ellipse_polyline(elipse2)) #returns points of the form AO,CO
    goodPoints = []
    for point in points:
        if(point[0]>=0 and point[1]>=0):
            OA = point[0]
            OC = point[1]
            OB = np.sqrt((2*OA*OA+2*OC*OC-AC*AC)/4)
            goodPoints.append([OA,OB,OC])
    return goodPoints
