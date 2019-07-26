from abaqus import *
from abaqusConstants import *
from numpy import *
from math import *
import os
os.chdir(r'C:\AbaTemp\VoroRock')
from  subrutine import *
from math import sqrt
import string

# Define variables
#Domain of the model BOUNDSIZE=(width,height)
domain=(1,1)
# Grain interface thickness
gt=0.01
# Voronoi seed file
seedname='seed.dat'
# Small edge length
error=0.05;
# Abaqus job name
jobname='GBM'
import numpy as np
points=np.loadtxt(seedname).tolist()
points=[tuple(i) for i in points]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#BOUNDSIZE=(width,height) C:\Users\MQX\Desktop\test\
EP=gt;BOUNDSIZE=domain;
bound = [(0., 0.), (0., BOUNDSIZE[1]), (BOUNDSIZE[0], BOUNDSIZE[1]), (BOUNDSIZE[0], 0.)]
center=(BOUNDSIZE[0]/2.,BOUNDSIZE[0]/2.)
nbound=bound

point2=[]
for i in range(0,len(points)):
  newpoint=points[i]
  for ix in (-1,1,0):
    for iy in (1,-1,0): 
      pointadd=(newpoint[0]+float(ix)*bound[2][0],newpoint[1]+float(iy)*bound[2][1])
      point2.append(pointadd)

nbound2 = [(-bound[2][0],-bound[2][1]), (-bound[2][0],2.*bound[2][1]), (2.*bound[2][0],2.*bound[2][1]), (2.*bound[2][0],-bound[2][1])]
voronoi2 = bounded_voronoi(nbound2, point2)
#voronoit=voronoi2
print 'Bound Voronoi Done'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getCurrentModel():
    vpName = session.currentViewportName
    modelName = session.sessionState[vpName]['modelName']
    return mdb.models[modelName]

myModel = getCurrentModel()
s = myModel.ConstrainedSketch(name='Sketch A', sheetSize=200.0)
s0 = myModel.ConstrainedSketch(name='Sketch Ini', sheetSize=200.0)
sr = myModel.ConstrainedSketch(name='Sketch Rect', sheetSize=200.0)
mask=1

sr.rectangle(point1=(0.0, 0.0), point2=(BOUNDSIZE[0], BOUNDSIZE[1]))
geomcarac=[BOUNDSIZE[0],BOUNDSIZE[1]]
center=(BOUNDSIZE[0]/2.,BOUNDSIZE[0]/2.)
#voronoi = verifMask(voronoi,points,mask,geomcarac)
#voronoi = verifVoronoi(voronoi,points,1*0.03)
#drawVoronoi(voronoi,bound,points,s0)
def regulization(voronoi, error):
    # Voronoi delete small edge
    #print str(len(voronoi))
    points=list(voronoi.keys())
    for i in range(0,len(points)):
        poly=voronoi[points[i]]
        for j in range(0,len(poly)):
          poly[j]=(round(poly[j][0],5),round(poly[j][1],5))
        voronoi[points[i]]=poly
    #print str(len(points))
    node_origin=points[0:0]
    node_new=points[0:0]
    for i in range(0,len(points)):
        poly=voronoi[points[i]]
        #print str(len(poly))
        poly.append(poly[0])
        for j in range(0,len(poly)-1):
            vect=(poly[j+1][0]-poly[j][0],poly[j+1][1]-poly[j][1])
            edge_length=sqrt(vect[0]**2+vect[1]**2)
            center=((poly[j+1][0]+poly[j][0])/2,(poly[j+1][1]+poly[j][1])/2)
            if edge_length<error:
                if poly[j] not in node_origin:
                    #print 'ok'
                    node_origin.append(poly[j]);node_new.append(center);
                if poly[j+1] not in node_origin:
                    node_origin.append(poly[j+1]);node_new.append(center);
        #print 'node_origin length is '+str(len(node_origin))
    for i in range(0,len(points)):
        poly=voronoi[points[i]]
        poly2=poly[0:0]
        for j in range(0,len(poly)):
            if poly[j] in node_origin:
                pt=node_new[node_origin.index(poly[j])]
                if pt not in poly2:
                  poly2.append(pt)
            else:
                if poly[j] not in poly2:
                  poly2.append(poly[j])
        voronoi[points[i]]=poly2
    return node_origin,node_new,voronoi

#error=0.05;
[node_origin,node_new,voronoi2]=regulization(voronoi2, error)
print 'node_origin length is '+str(len(node_origin))
while len(node_origin)>0:
	[node_origin,node_new,voronoi2]=regulization(voronoi2, error)
	print 'node_origin length is '+str(len(node_origin))

drawVoronoi(voronoi2,nbound2,point2,s0)
print 'Draw Voronoi Done'

myModel.convertAllSketches()
p = myModel.Part(name="IniGBM", dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
p.BaseShell(sketch=sr)
s.unsetPrimaryObject()
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s0)

# ADD Grain Interface
refpoint=[];	cellpoint=[]
for i in range(0,len(point2)):
	if (i%20==0):
		print 'Grain-num '+str(i)+'/'+str(len(point2))
	poly=[]
	poly = voronoi2[point2[i]] 
	lg=[]
	lg1=[]
	lg2=[]
	for j in range(0,len(poly)-1):
		if (point_on_bound(poly[j],nbound2) and point_on_bound(poly[j+1],nbound2)):	
			scalar=0.0
		else:
			scalar=1.0	
		
		pointmil=midpoint(poly[j],poly[j+1])
		vect=(poly[j+1][0]-poly[j][0], poly[j+1][1]-poly[j][1])
		vect=normalize(vect);		ml1=slope(poly[j],poly[j+1])
		ptl=(pointmil[0]+scalar*(EP/2.0)*vect[1],pointmil[1]-scalar*(EP/2.0)*vect[0] )
		lig=(ptl,ml1)
		lig1=(pointmil[0],pointmil[1],scalar)
		lg.append(lig)
		lg1.append(lig1)
	
	if(point_on_bound(poly[len(poly)-1],nbound2) and point_on_bound(poly[0],nbound2)): 
		scalar=0.0
	else:
		scalar=1.0
	
	pointmil=midpoint(poly[len(poly)-1],poly[0])
	vect=(poly[len(poly)-1][0]-poly[0][0], poly[len(poly)-1][1]-poly[0][1]) 
	vect=normalize(vect)
	ml1=slope(poly[len(poly)-1],poly[0]) 
	ptl=(pointmil[0]-scalar*(EP/2.0)*vect[1], pointmil[1]+scalar*(EP/2.0)*vect[0])
	lig=(ptl,ml1)
	lg.append(lig)
	
	newpointj=[]
	newpointf=[]
	for j in range(0,len(lg)-1):
		ap=line_intersect(lg[j], lg[j+1])
		if(ap!=None):		newpointj.append(ap)
	
	ap=line_intersect(lg[len(lg)-1], lg[0])
	if(ap!=None):
		newpointj.append(ap)
	# Judge the self intersect
	newpointf=check_poly(newpointj)
	if(len(newpointf)>0):
		for j in range(0,len(newpointf)-1):
		 s.Line(newpointf[j],newpointf[j+1])
		s.Line(newpointf[len(newpointf)-1],newpointf[0])

# After ADD Thickness
myModel.convertAllSketches()
p = myModel.Part(name='GBM', dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
p.BaseShell(sketch=sr)
s.unsetPrimaryObject()
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getCurrentViewport():

    vpName = session.currentViewportName
    return session.viewports[vpName]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getNewModel():

    vpName = session.currentViewportName
    modelName = session.sessionState[vpName]['modelName']
    if(len(mdb.models[modelName].parts)==0):
      return mdb.models[modelName]    
    modelName=modelName+"-1"  
    return mdb.Model(modelName)

vp = getCurrentViewport()
vp.setValues(displayedObject=p)
a = myModel.rootAssembly
a.Instance(name='GBM', part=p, dependent=ON)
dat = a.instances['GBM'].datums
vp.assemblyDisplay.setValues(interactions=ON,constraints=ON, connectors=ON, engineeringFeatures=ON)
vp.assemblyDisplay.setValues(renderStyle=WIREFRAME)
# Assign the set
p=mdb.models['Model-1'].parts['GBM']
f=p.faces
for i in range(len(f)):
	p.Set(name='set-'+str(i+1),faces=f[i:i+1])
mdb.saveAs(pathName=jobname + '.cae')

# Export to the PFC version 6.0
def direct_scale(line):
  if line[1]>line[3]:
    line=[line[2],line[3],line[0],line[1]];
  direct=[line[2]-line[0],line[3]-line[1]];
  scale=sqrt(direct[0]**2+direct[1]**2);
  dip=(-direct[0])/scale;
  dip=acos(dip)*180/pi;
  return dip,scale

#[dip,scale]=direct_scale([0,1,0,2])

def GBM_to_PFC(partname,filename,BOUNDSIZE=[1,1]):
  f=open(filename,'w+')
  f.write('model new\n')
  f.write('model title \"GBM Simulation\"\n')
  f.write('; Set the domain extent\n')
  f.write('model domain extent %4.2f %4.2f %4.2f %4.3f\n' % (-0.25*BOUNDSIZE[0],1.25*BOUNDSIZE[0],-0.25*BOUNDSIZE[1],1.25*BOUNDSIZE[1]))
  f.write('contact cmat default model linearpbond property kn 5e7 dp_nratio 0.5 \n')
  f.write('; Generate walls \n')
  f.write('wall generate box %4.2f %4.2f %4.2f %4.2f\n' % (0,BOUNDSIZE[0],0,BOUNDSIZE[1]))
  f.write('; Distribute balls in the box.\n')
  f.write('model random 1001\n')
  f.write('ball distribute porosity 0.15 radius %f %f box  %4.2f %4.2f %4.2f %4.3f\n' % (BOUNDSIZE[0]/125.0,BOUNDSIZE[1]/100.0,0,BOUNDSIZE[0],0,BOUNDSIZE[1]))
  f.write('; Set ball attributes\n')
  f.write('ball attribute density 1000.0 damp 0.7\n')
  f.write('; Calm the system\n')
  f.write('model cycle 1000 calm 10\n')
  f.write('; Solve the system to a target limit (here the average force ratio)\n')
  f.write('; Use density scaling to quickly reach equilibrium\n')
  f.write('model mechanical timestep scale\n')
  f.write('model solve ratio-average 1e-3\n')
  f.write('model mechanical timestep auto\n')
  f.write('model calm\n')
  f.write('; delete walls\n')
  f.write('wall delete\n')
  f.write('; Install parallel bonds to particles in contact \n')
  f.write('; assign very high strength to prevent breakage\n')
  f.write('contact method bond gap 0.0\n')
  f.write('contact property pb_kn 1e8 pb_ks 1e8 pb_ten 1e12 pb_coh 1e12 pb_fa 30.0 \n')
  f.write('; Reset ball displacement \n')
  f.write('ball attribute displacement multiply 0.0 \n')
  f.write('; Set linear stiffness to 0.0 and force a reset of linear contact forces. \n')
  f.write('contact property kn 0.0 lin_force 0.0 0.0 \n')
  f.write('ball attribute force-contact multiply 0.0 moment-contact multiply 0.0 \n')
  f.write('model cycle 1 \n')
  f.write('model solve ratio-average 1e-5 \n')
  f.write('model save \'intact\'\n')
  f.write('program call \'GBM_Fracture\'\n')
  f.write('fracture property ...\n')
  f.write('         \'sj_kn\' 2e9 \'sj_ks\' 2e9 \'sj_fric\' 0.70 ...\n')
  f.write('         \'sj_coh\' 0.0 \'sj_ten\' 0.0 \'sj_large\' 1\n')
  f.write('; Apply smoothjoint contact model to contacts intercepted by fracture\n')
  f.write('fracture contact-model model \'smoothjoint\' install \n')
  f.write('; Ensure new contacts intersecting the fracture are set to the sj contact model \n')
  f.write('fracture contact-model model \'smoothjoint\' activate \n')
  f.write('model save \'GBM\'\n')
  f.close()
  f=open('GBM_Fracture.p2dat','w+')
  p=mdb.models['Model-1'].parts[partname]
  v=p.vertices
  e=p.edges
  e0=e[0:0];e1=e0;
  for i in range(len(e)):
    ix=0
    et=e[i]
    vt=et.pointOn
    if abs(vt[0][0])<1e-3 or abs(vt[0][1])<1e-3 or abs(vt[0][0]-BOUNDSIZE[0])<1e-3 or abs(vt[0][1]-BOUNDSIZE[1])<1e-3:
      e0=e0+e[i:i+1]
    else:
      e1=e1+e[i:i+1]
  for i in range(len(e1)):
    vt=e1[i].getVertices()
    vl=[v[vt[0]].pointOn[0][0],v[vt[0]].pointOn[0][1],v[vt[1]].pointOn[0][0],v[vt[1]].pointOn[0][1]]
    pos=[vl[0]/2.0+vl[2]/2.0,vl[1]/2.0+vl[3]/2.0]
    [dip,scale]=direct_scale(vl)
    f.write('fracture create  ...\n')
    f.write('         position (%f, %f)  ...\n' % (pos[0],pos[1]))
    f.write('         dip %f  ...\n' % dip)
    f.write('         size %f  ...\n' % scale)
    f.write('         dfn  \"GBM\" \n')
  f.close()

partname='IniGBM'
filename='GBM.p2dat'
domain=[BOUNDSIZE[0],BOUNDSIZE[1]];
GBM_to_PFC(partname,filename,domain)


# Export to UDEC
def GBM_to_UDEC(partname,filename,BOUNDSIZE=[1,1]):
  f=open(filename,'w+')
  f.write('new\n')
  f.write('round 0.001\n')  
  f.write('block  0 0 0 %4.2f %4.2f %4.2f %4.3f 0\n' % (BOUNDSIZE[1],BOUNDSIZE[0],BOUNDSIZE[1],BOUNDSIZE[1]))
  p=mdb.models['Model-1'].parts[partname]
  v=p.vertices
  e=p.edges
  e0=e[0:0];e1=e0;
  for i in range(len(e)):
    ix=0
    et=e[i]
    vt=et.pointOn
    if abs(vt[0][0])<1e-3 or abs(vt[0][1])<1e-3 or abs(vt[0][0]-BOUNDSIZE[0])<1e-3 or abs(vt[0][1]-BOUNDSIZE[1])<1e-3:
      e0=e0+e[i:i+1]
    else:
      e1=e1+e[i:i+1]
  for i in range(len(e1)):
    vt=e1[i].getVertices()
    vl=[v[vt[0]].pointOn[0][0],v[vt[0]].pointOn[0][1],v[vt[1]].pointOn[0][0],v[vt[1]].pointOn[0][1]]
    f.write('crack (%f, %f) (%f, %f)\n' % (vl[0],vl[1],vl[2],vl[3]))
  f.close()

partname='IniGBM'
filename='GBM.uddat'
domain=[BOUNDSIZE[0],BOUNDSIZE[1]];
GBM_to_UDEC(partname,filename,domain)

