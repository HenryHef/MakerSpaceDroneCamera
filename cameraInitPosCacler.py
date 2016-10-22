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

def clusterTo6Sets(blackPointList):
	return find_centers(blackPointList,6)[0]

def get_CRef_DirectionVectors(points,W_max_angle,H_max_angle,pixH,pixW):
	#not that the direction vectors lengths have no meaning and are all normalized to 1
	print points
	h=np.tan(H_max_angle)
	w=np.tan(W_max_angle)
	vectorsUN = []
	vectors = []
	for point in points:
		k_w = 2*point[0]/pixW-1
		k_h = 2*point[1]/pixH-1
		#TODO check this step
		vectorsUN.append((1,w*k_w,h*k_h))
	for vectorUN in vectorsUN:
		norm = np.sqrt(vectorUN[0]**2+vectorUN[1]**2+vectorUN[2]**2)
		vectors.append((vectorUN[0]/norm,vectorUN[1]/norm,vectorUN[2]/norm))
	return vectors

def closeToParellel(v1,v2):
	x = v1[0]/v2[0]
	y = v1[1]/v2[1]
	z = v1[2]/v2[2]
	return np.abs(x-y)<x*(10**-2) and np.abs(x-z)<x*(10**-4)

def closeToCoplanerVectors(vSet):
	v1 = np.cross(vSet[0],vSet[1])
	v2 = np.cross(vSet[0],vSet[2])
	v3 = np.cross(vSet[1],vSet[2])
	return closeToParellel(v1,v2) and closeToParellel(v1,v3)

def findCloseToCoplanerVectors(vectors):
	result = []
	for vSet in itertools.combinations(vectors,3):
		if(closeToCoplanerVectors(vSet)):
			result.append((v1,v2,v3))
	return result

def dotproduct(v1, v2):
	return sum((a*b) for a, b in zip(v1, v2))

def length(v):
	return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
	return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

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

def get_1_0_Vector(dirVectors,cpVectorsSets,cpVecPairs):
	return list(set(dirVectors)-set(v for v in se in cpVectorsSets))[0]

def get_0_0_Vector(dirVectors,cpVectorsSets,cpVecPairs):
	vs = []
	for se in cpVectorsSets:
		for v in se:
			vs.append(v) 
	return list(set(vs)-set([cpVecPairs[0][2],cpVecPairs[0][4],cpVecPairs[1][2],cpVecPairs[1][4]]))[0]

def get_0_n1_vector(dirVectors,cpVectorsSets,cpVecPairs):
	endVectors = [cpVecPairs[0][2],cpVecPairs[0][4],cpVecPairs[1][2],cpVecPairs[1][4]]
	for vec in endVectors:
		if(closeToCoplanerVectors(list(set(endVectors)-set(vec)))):
			return vec

def get_0_1_vector(dirVectors,cpVectorsSets,cpVecPairs):
	if cpVectorsSets[0].contains(get_0_0_Vector(dirVectors,cpVectorsSets,cpVecPairs)):
		return list(set(cpVectorsSets[0])-set([get_0_0_Vector(dirVectors,cpVectorsSets,cpVecPairs),get_0_n1_Vector(dirVectors,cpVectorsSets,cpVecPairs)]))[0]
	else:
		return list(set(cpVectorsSets[1])-set([get_0_0_Vector(dirVectors,cpVectorsSets,cpVecPairs),get_0_n1_Vector(dirVectors,cpVectorsSets,cpVecPairs)]))[0]

def get_1_1_vectors(dirVectors,cpVectorsSets,cpVecPairs):
	return list(set(v for v in se in cpVectorsSets)-
			set([get_0_0_Vector(dirVectors,cpVectorsSets,cpVecPairs)
				,get_0_n1_Vector(dirVectors,cpVectorsSets,cpVecPairs)
				,get_0_1_Vector(dirVectors,cpVectorsSets,cpVecPairs)]))


def getVectorDistancePairMap(cpVecPairs,AC):
	dataList = []
	for vPair in cpVecPairs:
		AO_BO_CO = getAO_BO_CO_fromAngleBOA_BOC_AC(cpVecPairs[0],cpVecPairs[1],AC)
		dataList.append((cpVecPairs[0],AO_BO_CO[0]))
		dataList.append((cpVecPairs[1],AO_BO_CO[1]))
		dataList.append((cpVecPairs[2],AO_BO_CO[2]))
	return map(dataList)

def getLocation(centers,AC,pixW,pixH,max_angle_W,max_angle_H):
	print "step4.1"

	# these are in camera frame of reference
	dirVectors = get_CRef_DirectionVectors(centers,pixW,pixH,max_angle_W,max_angle_H)
	cpVectorsSets = findCloseToCoplanerVectors(dirVectors)
	cpVecPairs = [BOA_BOC_OAdir_OBdir_OCdir_from3Coplaner(cpVectorsSet) for cpVectorsSet in cpVectorsSets] #coplanerAngleDirVecParis

	# for now xy is wall, z is distance from wall. these are in camera frame of reference
	vec_0_0 = get_0_0_Vector(dirVectors,cpVectorsSets,cpVecPairs)
	vec_0_m1 = get_0_n1_vector(dirVectors,cpVectorsSets,cpVecPairs)
	vec_0_1 = get_0_1_vector(dirVectors,cpVectorsSets,cpVecPairs)
	vec_1_1a= get_1_1_vectors(dirVectors,cpVectorsSets,cpVecPairs)[0]
	vec_1_1b= get_1_1_vectors(dirVectors,cpVectorsSets,cpVecPairs)[1]
	vec_1_0 = get_1_0_Vector(dirVectors,cpVectorsSets,cpVecPairs)
	print "step4.2"

	vecDistMap = getVectorDistancePairMap(cpVecPairs,AC)
	dist_0_0=vecDistMap[vec_0_0]
	dist_0_01=vecDistMap[vec_0_1]
	y=(1-dist_0_1**2+n_0_0**2)/2
	dist_1_1a=vecDistMap[vec_1_1a]
	dist_1_1b=vecDistMap[vec_1_1b]
	xa = (1-dist_0_1**2+dist_1_1a**2)/2
	xb = (1-dist_0_1**2+dist_1_1b**2)/2
	za = np.sqrt(dist_0_0**2-xa**2-y**2)
	zb = np.sqrt(dist_0_0**2-xb**2-y**2)
	#these are in world space
	vec000a = (xa,y,za)
	vec000b = (xb,y,zb)
	vec100a = (xb-1,y,zb)
	vec100b = (xb-1,y,zb)
	angle_a = angle(vec000a,vec100a)
	angle_b = angle(vec000b,vec100b)
	angleReal000To100 = angle(vec_1_0,vec_0_0)
	x=-100000000
	z=-100000000
	if(np.abs(angle_a-angleReal000To100)<angleReal000To100*(10**-4)):#a is corect
		x=xa
		z=za
	elif(np.abs(angle_b-angleReal000To100)<angleReal000To100*(10**-4)):#a is corect
		x=xb
		z=zb
	else:
		print "NO SUCCESS!!!!"

	print "step4.3"
	#F is vector in camFramOfReference forward. It can be writen as av1+bv2+cv3, where v1,v2,v3 are 3 vectors towrds difference points.
	#find ABC in terms of v1,v2,v3 
	
	ABCVec=np.linalg.inv([vec_0_0,vec_1_0,vec_0_1]).dot([1,0,0])
	vec000=(x,y,z)
	vec100=(x-1,y,z)
	vec010=(x,y-1,z)
	forwardCarry=[vec000[i]*ABCVec[0]+vec100[i]*ABCVec[1]+vac010[i]*ABCVec[2] for i in range(0,3)]
	forward=[forwardCarry[i]/length(forwardCarry) for i in range(0,3)]

	print "step4.4"
	#same thing as above, but with top dir vec
	ABCVec=np.linalg.inv([vec_0_0,vec_1_0,vec_0_1]).dot([0,1,0])
	vec000=(x,y,z)
	vec100=(x-1,y,z)
	vec010=(x,y-1,z)
	topCarry=[vec000[i]*ABCVec[0]+vec100[i]*ABCVec[1]+vac010[i]*ABCVec[2] for i in range(0,3)]
	top=[forwardCarry[i]/length(forwardCarry) for i in range(0,3)]

	print "step4.5"
	return (x,y,z,forward,top)

def getBitmap(fileName):
	im = Image.open(fileName) #Can be many different formats.
	return (im.size[0],im.size[1],im.load())

def getBlackPointList(bitmap):#TODO add considering chunks to make k-means faster and to remove noise.
	points = []
	pix = bitmap[2]
	rMaxSq=20**2
	for i in range(bitmap[0]):
		for j in range(bitmap[1]):
			rgb = pix[i,j]
			if (rgb[0])**2+(rgb[1])**2+(255-rgb[2])**2<rMaxSq:
				points.append([i,j])
	return points

def findCameraLocation(fileName,pixW,pixH,max_angle_W,max_angle_H):
	print "step1"
	bitmap = getBitmap(fileName)
	print "step2"
	blackPointList = getBlackPointList(bitmap)
	print "step3"
	centersXYPixSpace = clusterTo6Sets(blackPointList)
	print "step4"
	location = getLocation(centersXYPixSpace,2,pixW,pixH,max_angle_W,max_angle_H)
	print "step5"
	return location

