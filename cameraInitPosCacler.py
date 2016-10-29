import numpy as np
import random
import itertools
from kMeansAlgoritham import find_centers
from elipseAproximateSolver import getAO_BO_CO_fromAngleBOA_BOC_AC
from PIL import Image
#import matplotlib.pyplot as plt

#from sklearn.cluster import KMeans
#from sklearn.datasets.samples_generator import make_blobs

#defining world axis to be x to the right, y is up and z is out of the wall

PARELLEL_SENSITIVITY = 10**-1
COPLANER_ANGLE_THRESHHOLD = 10/180.0*np.pi
PANGLE_EQUAL_SENSITIVITY=10**-2
DEBUG=False

def clusterTo6Sets(blackPointList):
	re=find_centers(blackPointList,6)[0]

	return re

def get_CRef_DirectionVectors(points,W_max_angle,H_max_angle,pixH,pixW):
	#not that the direction vectors lengths have no meaning and are all normalized to 1
	if(DEBUG):print "points = "+str(points)
	h=np.tan(H_max_angle/180.0*np.pi)
	w=np.tan(W_max_angle/180.0*np.pi)
	if(DEBUG):print "h="+str(h)+"  "+str(H_max_angle)
	vectorsUN = []
	vectors = []
	for point in points:
		k_w = 2.0*point[0]/pixH-1.0
		k_h = 2.0*point[1]/pixW-1.0
		#TODO check this step
		vectorsUN.append((1,w*k_w,h*k_h))
	for vectorUN in vectorsUN:
		norm = np.sqrt(vectorUN[0]**2+vectorUN[1]**2+vectorUN[2]**2)
		vectors.append((vectorUN[0]/norm,vectorUN[1]/norm,vectorUN[2]/norm))
	if(DEBUG):print "vectorsUN="+str(vectorsUN)
	if(DEBUG):print "vectors="+str(vectors)
	return vectors

def closeToParellel(v1,v2):
	if(DEBUG):print "v1="+str(v1)
	if(DEBUG):print "v2="+str(v2)
	x = v1[0]/v2[0]
	y = v1[1]/v2[1]
	z = v1[2]/v2[2]
	if(DEBUG):print "ratioVec="+str([x,y,z])
	re = np.abs(x-y)<np.max([x,y])*(PARELLEL_SENSITIVITY) and np.abs(x-z)<np.max([x,z])*(PARELLEL_SENSITIVITY)
	if(DEBUG):print "re="+str(re)
	return re

def closeToCoplanerVectors(vSet):
	v1 = np.cross(vSet[0],vSet[1])
	v2 = np.cross(vSet[0],vSet[2])
	v3 = np.cross(vSet[1],vSet[2])
	a12 = angle(v1,v2)
	a13 = angle(v1,v3)
	a23 = angle(v2,v3)
	if(a12>np.pi/2):a12=np.pi-a12
	if(a13>np.pi/2):a13=np.pi-a13
	if(a23>np.pi/2):a23=np.pi-a23
	if(DEBUG):print "a12,a13,a23="+str((a12/np.pi*180,a13/np.pi*180,a23/np.pi*180))+"      v1,v2,v3="+str((v1,v2,v3))
	return a12<COPLANER_ANGLE_THRESHHOLD and a13<COPLANER_ANGLE_THRESHHOLD and a23<COPLANER_ANGLE_THRESHHOLD

def findCloseToCoplanerVectors(vectors):
	result = []
	for vSet in itertools.combinations(vectors,3):
		if(DEBUG):print "vSet="+str(vSet)
		if(closeToCoplanerVectors(vSet)):
			result.append((vSet))
	if(DEBUG):print "clost to coplaner="+str(result)
	return result

def dotproduct(v1, v2):
	return sum((a*b) for a, b in zip(v1, v2))

def length(v):
	return np.sqrt(dotproduct(v, v))

def normalize(v):
	l=length(v)
	return [v[i]/l for i in range(0,len(v))]

def angle(v1, v2):
	v=dotproduct(v1, v2) / (length(v1) * length(v2))
	#if(DEBUG):print "v1,v2="+str((v1,v2))
	#if(DEBUG):print v
	if(np.abs(v+1.0)<.0001): return np.pi
	if(np.abs(v-1.0)<.0001): return 0
	r=np.arccos(v)
	#if(DEBUG):print "r="+str(r)
	return r

def BOA_BOC_OAdir_OBdir_OCdir_from3Coplaner(vecs):
	dot01 = np.dot(vecs[0],vecs[1])
	dot02 = np.dot(vecs[0],vecs[2])
	dot12 = np.dot(vecs[1],vecs[2])
	if dot01<dot02 and dot01<dot12:
		#v 0 and 1 are the furthest apart ---> A=[0], B=[2], C=[1]
		return (np.abs(angle(vecs[0],vecs[2])),np.abs(angle(vecs[1],vecs[2])),vecs[0],vecs[2],vecs[1])
	elif dot02<dot12:
		#v 0 and 2 are the furthest apart ---> A=[0], B=[1], C=[2]
		return (np.abs(angle(vecs[0],vecs[1])),np.abs(angle(vecs[1],vecs[2])),vecs[0],vecs[1],vecs[2])
	else:
		#v 1 and 2 are the furthest apart ---> A=[1], B=[0], C=[2]
		return (np.abs(angle(vecs[0],vecs[1])),np.abs(angle(vecs[0],vecs[2])),vecs[1],vecs[0],vecs[2])

def getVectorOrderedList(dirVectors,cpVectorsSets,cpVecPairs):
	#return [(0,0) , (0,-1) , (0,1) , (1,0) , (1,1)a , (1,1)b]
	if(DEBUG):print""
	if(DEBUG):print ""
	if(DEBUG):print "big new method"
	if(DEBUG):print "dirVectors="+str(dirVectors)
	if(DEBUG):print "cpVectorsSets="+str(cpVectorsSets)
	if(DEBUG):print "cpVecPairs="+str(cpVecPairs)
	setOfCPVectorSetVectors = set()
	for se in cpVectorsSets:
		for v in se:
			setOfCPVectorSetVectors.add(v)

	setDirVectors = set()
	for v in dirVectors:
		setDirVectors.add(v)

	setOfCPVecPairEnds = set()
	setOfCPVecPairEnds.add(cpVecPairs[0][2])
	setOfCPVecPairEnds.add(cpVecPairs[0][4])
	setOfCPVecPairEnds.add(cpVecPairs[1][2])
	setOfCPVecPairEnds.add(cpVecPairs[1][4])

	vec00 = list(setOfCPVectorSetVectors-setOfCPVecPairEnds)
	if(DEBUG):print "vec00(s)="+str(vec00)
	vec00=vec00[0]
	if(DEBUG):print "vec00="+str(vec00)

	if(DEBUG):print "len(dirVectorSet)="+str(len(setDirVectors))
	if(DEBUG):print "len(setOfCPVectorSetVectors)="+str(len(setOfCPVectorSetVectors))
	if(DEBUG):print "dirVectorSet="+str(setDirVectors)
	if(DEBUG):print "setOfCPVectorSetVectors="+str(setOfCPVectorSetVectors)

	vec10=list(setDirVectors-setOfCPVectorSetVectors)
	if(DEBUG):print "vec10(s)="+str(vec10)
	vec10=vec10[0]
	if(DEBUG):print "vec10="+str(vec10)

	for se in cpVectorsSets:
		if vec00 in se:
			vec00CPSet = se
		else:
			notVec00CPSet = se
	setOfVectors_vec00CPSet = set(v for v in vec00CPSet)
	setOfVectors_notVec00CPSet = set(v for v in notVec00CPSet)

	if(DEBUG):print "len(setOfCPVecPairEnds)="+str(len(setOfCPVecPairEnds))
	if(DEBUG):print "len(setOfVectors_notVec00CPSet)="+str(len(setOfVectors_notVec00CPSet))
	if(DEBUG):print "setOfCPVecPairEnds="+str(setOfCPVecPairEnds)
	if(DEBUG):print "setOfVectors_notVec00CPSet="+str(setOfVectors_notVec00CPSet)

	vec0n1 = list(setOfCPVecPairEnds-setOfVectors_notVec00CPSet)
	if(DEBUG):print "vec0n1(s)="+str(vec0n1)
	vec0n1=vec0n1[0]
	if(DEBUG):print "vec0n1="+str(vec0n1)

	vec01 = list(setOfVectors_vec00CPSet-set([vec00])-set([vec0n1]))
	if(DEBUG):print set([vec00])
	if(DEBUG):print vec00
	if(DEBUG):print "vec01(s)="+str(vec01)
	vec01=vec01[0]
	if(DEBUG):print "vec01="+str(vec01)

	vec11s = list(setDirVectors-set([vec00])-set([vec0n1])-set([vec10])-set([vec01]))
	vec11a=vec11s[0]
	vec11b=vec11s[1]

	if(DEBUG):print""
	if(DEBUG):print ""

	return [vec00 , vec0n1 , vec01 , vec10 , vec11a , vec11b]
	#0,0 vector is setOfCPVectorSetVectors-setOfCPVecPairEnds

def get_1_0_vector(dirVectors,cpVectorsSets,cpVecPairs):
	listVecs = []
	for i in range(len(cpVectorsSets)):
		for j in range(len(cpVectorsSets[i])):
			listVecs.append(cpVectorsSets[i][j])
	return list(set(dirVectors)-set(listVecs))[0]

def get_0_0_vector(dirVectors,cpVectorsSets,cpVecPairs):
	getVectorOrderedList(dirVectors,cpVectorsSets,cpVecPairs)
	if(DEBUG):print "get_0_0_vector"
	if(DEBUG):print dirVectors
	if(DEBUG):print cpVectorsSets
	if(DEBUG):print cpVecPairs
	vs = []
	for se in cpVectorsSets:
		for v in se:
			vs.append(v) 
	return list(set(vs)-set([cpVecPairs[0][2],cpVecPairs[0][4],cpVecPairs[1][2],cpVecPairs[1][4]]))[0]

def get_0_n1_vector(dirVectors,cpVectorsSets,cpVecPairs):
	endVectors = [cpVecPairs[0][2],cpVecPairs[0][4],cpVecPairs[1][2],cpVecPairs[1][4]]
	if(DEBUG):print "ENDVEC="+str(endVectors)
	for vec in endVectors:
		lis = []
		for v in endVectors:
			if v!= vec:
				lis.append(v)
		if(DEBUG):print "LIST="+str(lis)
		if(closeToCoplanerVectors(lis)):
			return vec

def get_0_1_vector(dirVectors,cpVectorsSets,cpVecPairs):
	if (get_0_0_vector(dirVectors,cpVectorsSets,cpVecPairs)in cpVectorsSets[0]):
		return list(set(cpVectorsSets[0])-set([get_0_0_vector(dirVectors,cpVectorsSets,cpVecPairs),get_0_n1_vector(dirVectors,cpVectorsSets,cpVecPairs)]))[0]
	else:
		return list(set(cpVectorsSets[1])-set([get_0_0_vector(dirVectors,cpVectorsSets,cpVecPairs),get_0_n1_vector(dirVectors,cpVectorsSets,cpVecPairs)]))[0]

def get_1_1_vectors(dirVectors,cpVectorsSets,cpVecPairs):
	listVecs = []
	for i in range(len(cpVectorsSets)):
		for j in range(len(cpVectorsSets[i])):
			listVecs.append(cpVectorsSets[i][j])
	return list(i for i in (set(listVecs)-
			set([get_0_0_vector(dirVectors,cpVectorsSets,cpVecPairs)
				,get_0_n1_vector(dirVectors,cpVectorsSets,cpVecPairs)
				,get_0_1_vector(dirVectors,cpVectorsSets,cpVecPairs)])))


def getVectorDistancePairMap(cpVecPairs,AC):
	dataMap={}
	for vPair in cpVecPairs:
		if(DEBUG):print "cpVecPairs="+str(cpVecPairs)
		AO_BO_CO = getAO_BO_CO_fromAngleBOA_BOC_AC(vPair[0],vPair[1],AC)
		if(DEBUG):print "AO_BO_CO="+str(AO_BO_CO)
		dataMap[vPair[2]]=AO_BO_CO[0][0]
		dataMap[vPair[3]]=AO_BO_CO[0][1]
		dataMap[vPair[4]]=AO_BO_CO[0][2]
	return dataMap

def getLocation(centers,AC,pixW,pixH,max_angle_W,max_angle_H):
	if(DEBUG):print "step4.1"

	# these are in camera frame of reference
	dirVectors = get_CRef_DirectionVectors(centers,max_angle_W,max_angle_H,pixW,pixH)
	cpVectorsSets = findCloseToCoplanerVectors(dirVectors)
	cpVecPairs = [BOA_BOC_OAdir_OBdir_OCdir_from3Coplaner(cpVectorsSet) for cpVectorsSet in cpVectorsSets] #coplanerAngleDirVecParis

	if(DEBUG):print "LENTHD DIR VEC RET = "+str(len(dirVectors))

	# for now xy is wall, z is distance from wall. these are in camera frame of reference
	vecs = getVectorOrderedList(dirVectors,cpVectorsSets,cpVecPairs)#[vec00 , vec0n1 , vec01 , vec10 , vec11a , vec11b]
	vec_0_0 = vecs[0]
	vec_0_m1 = vecs[1]
	vec_0_1 = vecs[2]
	vec_1_0 = vecs[3]
	vec_1_1a= vecs[4]
	vec_1_1b= vecs[5]
	if(DEBUG):print "step4.2"

	vecDistMap = getVectorDistancePairMap(cpVecPairs,AC)
	dist_0_0=vecDistMap[vec_0_0]
	dist_0_1=vecDistMap[vec_0_1]
	y=(1-dist_0_1**2+dist_0_0**2)/2
	dist_1_1a=vecDistMap[vec_1_1a]
	dist_1_1b=vecDistMap[vec_1_1b]
	if(DEBUG):print "dist_1_1b="+str(dist_1_1b)
	if(DEBUG):print "dist_1_1a="+str(dist_1_1a)
	if(DEBUG):print "dist_0_1="+str(dist_0_1)
	if(DEBUG):print "dist_0_0="+str(dist_0_0)
	if(DEBUG):print "y="+str(y)
	xa = (1+dist_0_1**2-dist_1_1a**2)/2
	xb = (1+dist_0_1**2-dist_1_1b**2)/2
	xc=-xa
	xd=-xb
	if(DEBUG):print "xa="+str(xa)
	if(DEBUG):print "xb="+str(xb)
	if(DEBUG):print "dist_0_0**2-xa**2-y**2="+str((dist_0_0**2)-(xa**2)-(y**2))
	if(DEBUG):print "dist_0_0**2-xb**2-y**2="+str((dist_0_0**2)-(xb**2)-(y**2))
	za = np.sqrt((dist_0_0**2)-(xa**2)-(y**2))
	zb = np.sqrt((dist_0_0**2)-(xb**2)-(y**2))
	zc = np.sqrt((dist_0_0**2)-(xc**2)-(y**2))
	zd = np.sqrt((dist_0_0**2)-(xd**2)-(y**2))
	#these are in world space
	vec000a = (xa,y,za)
	vec000b = (xb,y,zb)
	vec000c = (xc,y,zc)
	vec000d = (xd,y,zd)
	vec100a = (xa-1,y,za)
	vec100b = (xb-1,y,zb)
	vec100c = (xc-1,y,zc)
	vec100d = (xd-1,y,zd)
	if(DEBUG):print "(vec000a,vec100a,vec000b,vec100b)="+str((vec000a,vec100a,vec000b,vec100b))
	angle_a = angle(vec000a,vec100a)
	angle_b = angle(vec000b,vec100b)
	angleReal000To100 = angle(vec_1_0,vec_0_0)
	x=-100000000
	z=-100000000
	if(DEBUG):print "(angle_a,angle_b,angleReal000To100)="+str((angle_a,angle_b,angleReal000To100))
	vec_1_1=[]
	vec_n1_1=[]
	if(np.abs(angle_a-angleReal000To100)<angleReal000To100*(PANGLE_EQUAL_SENSITIVITY)):#a is corect
		x=xa
		z=za
		vec_1_1=vec_1_1a
		vec_n1_1=vec_1_1b
	elif(np.abs(angle_b-angleReal000To100)<angleReal000To100*(PANGLE_EQUAL_SENSITIVITY)):#b is corect
		x=xb
		z=zb
		vec_1_1=vec_1_1b
		vec_n1_1=vec_1_1a
	else:
		if(DEBUG):print "NO SUCCESS!!!!"

	if(DEBUG):print "step4.3"
	#T=WC^-1 where C is colum matrix camera space and W is matrix world space 
	vec000 = normalize([x,y,z])
	vec100 = normalize([x-1,y,z])
	vec010 = normalize([x,y-1,z])

	if(DEBUG):print "vec000 len = "+str(length(vec000)) + "    vec_0_0 len = "+str(length(vec_0_0))
	if(DEBUG):print "vec000 len = "+str(length(vec100)) + "    vec_0_0 len = "+str(length(vec_1_0))
	if(DEBUG):print "vec000 len = "+str(length(vec010)) + "    vec_0_0 len = "+str(length(vec_0_1))
	W=columnMat3(vec_0_0,vec_1_0,vec_0_1)
	C=columnMat3(vec000,vec100,vec010)
	W_n1=np.linalg.inv(W)
	C_n1=np.linalg.inv(C)
	T=np.dot(C,W_n1)
	T_T=np.matrix.transpose(T)
	if(DEBUG):print "!!!!!!!!!! T_T="+str(T_T)
	forward = T_T[0]
	right = T_T[1]
	up = T_T[2]

	if(DEBUG):print "step4.5"
	if(DEBUG):print "map'"+str(vecDistMap)
	vecXYDistMap = {}
	vecXYDistMap[(0,0)]=vecDistMap[vec_0_0]
	vecXYDistMap[(0,1)]=vecDistMap[vec_0_1]
	if(DEBUG):print "vec_0_m1="+str(vec_0_m1)
	vecXYDistMap[(0,-1)]=vecDistMap[vec_0_m1]
	if(DEBUG):print "vec_1_1"+str(vec_1_1)
	vecXYDistMap[(1,1)]=vecDistMap[vec_1_1]
	vecXYDistMap[(-1,1)]=vecDistMap[vec_n1_1]
	checkLoc(x,y,z,vecXYDistMap)
	return (x,y,z,forward,up)

def getBitmapData(fileName):
	im = Image.open(fileName) #Can be many different formats.
	return (im.size[0],im.size[1],im.load())

def getBlackPointList(bitmapData):#TODO add considering chunks to make k-means faster and to remove noise.
	points = []
	pix = bitmapData[2]
	rMaxSq=20**2
	img = Image.open('output1.jpg')
	outputData = img.load()
	for i in range(bitmapData[0]):
		for j in range(bitmapData[1]):
			rgb = pix[i,j]
			if (rgb[0])**2+(rgb[1])**2+(255-rgb[2])**2<rMaxSq:
				points.append([i,j])
				outputData[i,j]=(0,255,0)
			else:
				outputData[i,j]=rgb
	img.save("output1.jpg", "JPEG")
	return points

def findCameraLocation(fileName,max_angle_W,max_angle_H):
	if(DEBUG):print "step1"
	data = getBitmapData(fileName)
	if(DEBUG):print "step2"
	blackPointList = getBlackPointList(data)
	if(DEBUG):print "step3"
	centersXYPixSpace = clusterTo6Sets(blackPointList)
	if(DEBUG):print "step4"

	img = Image.open('output2.png')
	outputData = img.load()
	intXYCens = [(np.round(x),np.round(y)) for x,y in centersXYPixSpace]
	if(DEBUG):print "intcenters = "+str(intXYCens)
	for i in range(data[0]):
		for j in range(data[1]):
			rgb = data[2][i,j]
			if (i,j) in intXYCens:
				if(DEBUG):print "it is in"
				outputData[i,j]=(0,255,0)
			else:
				outputData[i,j]=rgb
	img.save("output2.png", "PNG")


	location = getLocation(centersXYPixSpace,2,data[0],data[1],max_angle_W,max_angle_H)
	if(DEBUG):print "step5"+str(location)
	if(DEBUG):print location+(data[0],data[1])
	return location+(data[0],data[1])

def checkLoc(x,y,z,vecXYDistMap):
	vecs=[(0,0),(0,1),(0,-1),(1,1),(-1,1)]
	for vec in vecs:
		distXYZ = np.sqrt((x-vec[0])**2+(y-vec[1])**2+z**2)
		distMap = vecXYDistMap[vec]
		if(DEBUG):print "vec="+str(vec)+": ditsXYZ="+str(distXYZ)+" mapDist = "+str(distMap)

def columnMat3(v1,v2,v3):
	return [[v1[i],v2[i],v3[i]] for i in range(0,3)]

#if(DEBUG):print "LOC="+str(getLocation([(10,10),(10,15),(10,5),(5,5),(15,5),(15,10)],2,30,30,45,45))