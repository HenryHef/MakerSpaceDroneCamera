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
PANGLE_EQUAL_SENSITIVITY=2*10**-1


def clusterTo6Sets(blackPointList):
	re=find_centers(blackPointList,6)[0]

	return re

def get_CRef_DirectionVectors(points,W_max_angle,H_max_angle,pixH,pixW):
	#not that the direction vectors lengths have no meaning and are all normalized to 1
	print "points = "+str(points)
	h=np.tan(H_max_angle/180.0*np.pi)
	w=np.tan(W_max_angle/180.0*np.pi)
	print "h="+str(h)+"  "+str(H_max_angle)
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
	print "vectorsUN="+str(vectorsUN)
	return vectors

def closeToParellel(v1,v2):
	print "v1="+str(v1)
	print "v2="+str(v2)
	x = v1[0]/v2[0]
	y = v1[1]/v2[1]
	z = v1[2]/v2[2]
	print "ratioVec="+str([x,y,z])
	re = np.abs(x-y)<np.max([x,y])*(PARELLEL_SENSITIVITY) and np.abs(x-z)<np.max([x,z])*(PARELLEL_SENSITIVITY)
	print "re="+str(re)
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
	print "a12,a13,a23="+str((a12/np.pi*180,a13/np.pi*180,a23/np.pi*180))
	return a12<COPLANER_ANGLE_THRESHHOLD and a13<COPLANER_ANGLE_THRESHHOLD and a23<COPLANER_ANGLE_THRESHHOLD

def findCloseToCoplanerVectors(vectors):
	result = []
	for vSet in itertools.combinations(vectors,3):
		if(closeToCoplanerVectors(vSet)):
			result.append((vSet))
	print "clost to coplaner="+str(result)
	return result

def dotproduct(v1, v2):
	return sum((a*b) for a, b in zip(v1, v2))

def length(v):
	return np.sqrt(dotproduct(v, v))

def angle(v1, v2):
	return np.arccos(dotproduct(v1, v2) / (length(v1) * length(v2)))

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

def get_1_0_vector(dirVectors,cpVectorsSets,cpVecPairs):
	listVecs = []
	for i in range(len(cpVectorsSets)):
		for j in range(len(cpVectorsSets[i])):
			listVecs.append(cpVectorsSets[i][j])
	return list(set(dirVectors)-set(listVecs))[0]

def get_0_0_vector(dirVectors,cpVectorsSets,cpVecPairs):
	print "get_0_0_vector"
	print dirVectors
	print cpVectorsSets
	print cpVecPairs
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
		print "cpVecPairs="+str(cpVecPairs)
		AO_BO_CO = getAO_BO_CO_fromAngleBOA_BOC_AC(vPair[0],vPair[1],AC)
		print "AO_BO_CO="+str(AO_BO_CO)
		dataMap[vPair[2]]=AO_BO_CO[0][0]
		dataMap[vPair[3]]=AO_BO_CO[0][1]
		dataMap[vPair[4]]=AO_BO_CO[0][2]
	return dataMap

def getLocation(centers,AC,pixW,pixH,max_angle_W,max_angle_H):
	print "step4.1"

	# these are in camera frame of reference
	dirVectors = get_CRef_DirectionVectors(centers,pixW,pixH,max_angle_W,max_angle_H)
	cpVectorsSets = findCloseToCoplanerVectors(dirVectors)
	cpVecPairs = [BOA_BOC_OAdir_OBdir_OCdir_from3Coplaner(cpVectorsSet) for cpVectorsSet in cpVectorsSets] #coplanerAngleDirVecParis

	# for now xy is wall, z is distance from wall. these are in camera frame of reference
	vec_0_0 = get_0_0_vector(dirVectors,cpVectorsSets,cpVecPairs)
	vec_0_m1 = get_0_n1_vector(dirVectors,cpVectorsSets,cpVecPairs)
	vec_0_1 = get_0_1_vector(dirVectors,cpVectorsSets,cpVecPairs)
	vec_1_1a= get_1_1_vectors(dirVectors,cpVectorsSets,cpVecPairs)[0]
	vec_1_1b= get_1_1_vectors(dirVectors,cpVectorsSets,cpVecPairs)[1]
	vec_1_0 = get_1_0_vector(dirVectors,cpVectorsSets,cpVecPairs)
	print "step4.2"

	vecDistMap = getVectorDistancePairMap(cpVecPairs,AC)
	dist_0_0=vecDistMap[vec_0_0]
	dist_0_1=vecDistMap[vec_0_1]
	y=(1-dist_0_1**2+dist_0_0**2)/2
	dist_1_1a=vecDistMap[vec_1_1a]
	dist_1_1b=vecDistMap[vec_1_1b]
	xa = (1-dist_0_1**2+dist_1_1a**2)/2
	xb = (1-dist_0_1**2+dist_1_1b**2)/2
	za = np.sqrt(dist_0_0**2-xa**2-y**2)
	zb = np.sqrt(dist_0_0**2-xb**2-y**2)
	#these are in world space
	vec000a = (xa,y,za)
	vec000b = (xb,y,zb)
	vec100a = (xb-1,y,za)
	vec100b = (xb-1,y,zb)
	print "(vec000a,vec100a,vec000b,vec100b)="+str((vec000a,vec100a,vec000b,vec100b))
	angle_a = angle(vec000a,vec100a)
	angle_b = angle(vec000b,vec100b)
	angleReal000To100 = angle(vec_1_0,vec_0_0)
	x=-100000000
	z=-100000000
	print "(angle_a,angle_b,angleReal000To100)="+str((angle_a,angle_b,angleReal000To100))
	if(np.abs(angle_a-angleReal000To100)<angleReal000To100*(PANGLE_EQUAL_SENSITIVITY)):#a is corect
		x=xa
		z=za
	elif(np.abs(angle_b-angleReal000To100)<angleReal000To100*(PANGLE_EQUAL_SENSITIVITY)):#b is corect
		x=xb
		z=zb
	else:
		print "NO SUCCESS!!!!"

	print "step4.3"
	#F is vector in camFramOfReference forward. It can be writen as av1+bv2+cv3, where v1,v2,v3 are 3 vectors towrds difference points.
	#find ABC in terms of v1,v2,v3 
	mat = [vec_0_0,vec_1_0,vec_0_1]
	print "mat="+str(mat)
	ABCVec=np.linalg.inv(mat).dot([1,0,0])
	vec000=(x,y,z)
	vec100=(x-1,y,z)
	vec010=(x,y-1,z)
	forwardCarry=[vec000[i]*ABCVec[0]+vec100[i]*ABCVec[1]+vec010[i]*ABCVec[2] for i in range(0,3)]
	forward=[forwardCarry[i]/length(forwardCarry) for i in range(0,3)]

	print "step4.4"
	#same thing as above, but with top dir vec
	ABCVec=np.linalg.inv([vec_0_0,vec_1_0,vec_0_1]).dot([0,1,0])
	vec000=(x,y,z)
	vec100=(x-1,y,z)
	vec010=(x,y-1,z)
	topCarry=[vec000[i]*ABCVec[0]+vec100[i]*ABCVec[1]+vec010[i]*ABCVec[2] for i in range(0,3)]
	top=[forwardCarry[i]/length(forwardCarry) for i in range(0,3)]

	print "step4.5"
	return (x,y,z,forward,top)

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
	print "step1"
	data = getBitmapData(fileName)
	print "step2"
	blackPointList = getBlackPointList(data)
	print "step3"
	centersXYPixSpace = clusterTo6Sets(blackPointList)
	print "step4"

	img = Image.open('output2.png')
	outputData = img.load()
	intXYCens = [(np.round(x),np.round(y)) for x,y in centersXYPixSpace]
	print "intcenters = "+str(intXYCens)
	for i in range(data[0]):
		for j in range(data[1]):
			rgb = data[2][i,j]
			if (i,j) in intXYCens:
				print "it is in"
				outputData[i,j]=(0,255,0)
			else:
				outputData[i,j]=rgb
	img.save("output2.png", "PNG")


	location = getLocation(centersXYPixSpace,2,max_angle_W,max_angle_H,data[0],data[1])
	print "step5"+str(location)
	print location+(data[0],data[1])
	return location+(data[0],data[1])

