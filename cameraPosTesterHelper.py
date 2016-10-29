import numpy as np
import random
import numpy as np
from cameraInitPosCacler import getLocation
from cameraInitPosCacler import normalize,angle,length
#import matplotlib.pyplot as plt

#from sklearn.cluster import KMeans
#from sklearn.datasets.samples_generator import make_blobs

#defining world axis to be x to the right, y is up and z is out of the wall
def uniformCircPoint(radMax):
	radSqrt = np.sqrt(random.random())*radMax
	theta = random.random()*2*np.pi
	return [radSqrt*np.cos(theta),radSqrt*np.sin(theta)]
def genTestDataSet(x,y,z):#for now, forward and up are defult
	pointsBase = [(0.0,0.0),(0.0,1.0),(0.0,-1.0),(1.0,0.0),(1.0,1.0),(-1.0,1.0)]
	rNoise = .003
	points = [np.add(point,uniformCircPoint(rNoise)) for point in pointsBase]
	vecsW = [(x-p[0],y-p[1],z) for p in points]
	#print "vecsW="+str(vecsW)

	forwardXY=[.2*(2*random.random()-1.0),0]#.2*(2*random.random()-1.0)]
	forward=[forwardXY[0],forwardXY[1],np.sqrt(1-forwardXY[0]**2-forwardXY[1]**2)]
	theta = random.random()*2*np.pi
	right=normalize(np.cross(forward,[0,np.cos(theta),np.sin(theta)]))
	up=np.cross(forward,right)
	#print angle(forward,right)
	#print length(forward)
	#print length(right)
	V=np.matrix.transpose(np.asarray([right,up,forward]))
	V_n1 = np.linalg.inv(V)

	#print forward
	#print right
	#print up

	vecs = [np.dot(V_n1,w) for w in vecsW]
	#print "vecs Transformed="+str(vecs)

	vecsUnitZ = [(v[0]/v[2],v[1]/v[2],1) for v in vecs]
	#print "vecsUnitZ="+str(vecsUnitZ)
	
	maxX = max(abs(v[0]) for v in vecsUnitZ)
	maxY = max(abs(v[1]) for v in vecsUnitZ)
	sizeX = maxX*(5.0+random.random())
	sizeY = maxY*(5.0+random.random())
	angleX = np.arctan(sizeX)
	angleY = np.arctan(sizeY)
	pixX = 4000.0
	pixY = 4800.0
	pixPos = [(pixX/2.0+v[0]/sizeX*pixX/2,pixY/2.0+v[1]/sizeY*pixY/2) for v in vecsUnitZ]
	#print pixPos
	return (pixX,pixY,angleX*180.0/np.pi,angleY*180.0/np.pi,x,y,z,pixPos,vecsUnitZ)

def genTest():#for now, forward and up are defult
	return genTestDataSet(random.random()*3-1.5,random.random()*3-1.5,random.random()*8+2)

test = genTest()
print "test was = "+str(test)
print getLocation(test[7],2,test[0],test[1],test[2],test[3])
print "test was = "+str((test[4],test[5],test[6]))
