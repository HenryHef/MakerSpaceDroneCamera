from cameraInitPosCacler import findCameraLocation
#Note that both forwardVector and upVectro are normalized and (should be)
#perpendicular, so the rightVector is forwardVector cross upVector=rightVector
class Camera:
    x=0
    y=0
    z=0
    forwardVector=[0,0,0]
    upVector=[0,0,0]
    pixW=0
    pixH=0
    max_angle_W=0
    max_angle_H=0
    def __init__(self,pixW,pixH,max_angle_W,max_angle_H):
        self.pixW=pixW
        self.pixH=pixH
        self.max_angle_W=max_angle_W
        self.max_angle_H=max_angle_H
    def loadPositionFromImage(self,fileName):
        data = findCameraLocation(fileName,self.pixW,self.pixH,self.max_angle_W,self.max_angle_H)
        self.x=data[0]
        self.y=data[1]
        self.z=data[2]
        self.forwardVector=data[3]
        self.upVector=data[4]
    def mapPointToDirWorldFromCamera(x,y):
        w=np.tan(self.W_max_angle)
        h=np.tan(self.H_max_angle)
        k_w = 2*x/self.pixW-1
        k_h = 2*y/self.pixH-1

        hVector=[upVector[i]*h*k_h for i in range(0,len(upVector))]
        wVector=[np.cross(forwardVector,upVector)[i]*w*k_w for i in range(0,3)]
        
            #TODO check this step
        vectorUN = [forwardVector[i]+hVector[i]+wVector[i] for i in range(0,3)]#forwardVector+(1,w*k_w,h*k_h)

        norm = np.sqrt(sum(vectorUN[i]**2 for i in range(0,3)))
        vector=[vectorUN[i]/norm for i in range (0,3)]
        return vector

    def __str__(self):
        return '(x,y,z)=('+str(self.x)+","+str(self.y)+","+str(self.z)+")    forwardVector="+str(self.forwardVector)+"    upVector="+str(self.upVector)\
            +"    pixW="+str(self.pixW) +"    pixH="+str(self.pixH)+"    max_angle_W="+str(self.max_angle_W)+"    max_angle_H="+str(self.max_angle_H)

print Camera(1000,1000,30,30).loadPositionFromImage("testBackup.jpg")  #field angles are 65 and 50