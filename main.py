# -*- coding: utf-8 -*-  
#-----------------------------------------------------------------------
#     Authors: Qing-Xiang Meng 
#     Institute: Hohai University
#     Date: 2-July-2019
#     E-mail:mqx4088@gmail.com
#     Acknowledgements: The author would like to thank 
#     Stephane Lejeunes in LMA CNRS, France for the python script for Voronoi tessellation
#-----------------------------------------------------------------------
from abaqus import *
from abaqusConstants import *
from numpy import *
from math import *
from math import sqrt
import string
import random
import os
# Parameters for the Path
# Warning: you should change the path D:\Temp\VoroRock 
# to your own beofore running the code
dirpath='D:\Temp\VoroRock'

# Parameters for columnar jointed rock ; See detail in the paper
Jd=20.15;Jt=0.01;CV=0.40 #0.4418
# Prameters for the model domain
domain=(1,1)
# parameter for periodicity 1 is with periodic 0 is no periodity
isperiodic=1
#  Convergence error of CV
cv_error=0.01
# Small edge remove error using regulization
edge_error=0.05

os.chdir(dirpath)
from subrutine import *

'''
Part-1: Generate the seeds and voronoi diagram with the specified Jd Jt and CV using modified Lloyd's iteration
cv_error is judging the convergence of iteration
edge_error is judging the regulization
rcv is the cv during the iteration
'''
# Modified Lloyd's iteration to cv
#cv_error=0.01 #0.002
points=generate_random_seed(domain,Jd)
bound = [(0., 0.), (0., domain[1]), (domain[0], domain[1]), (domain[0], 0.)]
[points,rcv,iter_num]=iteration_lloyd_cv(bound,points,CV,cv_error);
print ('Convergenced with %d iterations' % (iter_num) )
print 'CV during iteration is:'
print rcv
# Generated periodic structure
if isperiodic==1:
	point2=generate_new_seed(domain,points)
	nbound2=generate_new_domain(domain)
else:
	point2=points
	nbound2=bound
voronoi2 = bounded_voronoi(nbound2, point2)
# Remove small edge using regulization
#edge_error=0.05
[node_origin,node_new,voronoi2]=regulization(voronoi2,nbound2, edge_error)
while len(node_origin)>0:
	[node_origin,node_new,voronoi2]=regulization(voronoi2,nbound2, edge_error)
# Add thickness based on Jt using polygon offset
offset=polygon_offset(voronoi2,point2,nbound2,Jt)

'''
Part-2: Generate the numerical model in Abaqus
Three sketches are generate s0 s1 sr
s0 for the initial diagram without offset
s1 for the diagram after offset
sd for the domain
'''
vpName = session.currentViewportName
modelName = session.sessionState[vpName]['modelName']
myModel = mdb.models[modelName]
s0 = myModel.ConstrainedSketch(name='Sketch Ini', sheetSize=200.0)
s1 = myModel.ConstrainedSketch(name='Sketch offset', sheetSize=200.0)
sr = myModel.ConstrainedSketch(name='Sketch Domain', sheetSize=200.0)
# Draw domain in the Abaqus
sr.rectangle(point1=(0.0, 0.0), point2=(domain[0], domain[1]))
# Draw initial digram in Abaqus
drawVoronoi(voronoi2,s0)
# Draw offset digram in Abaqus
drawVoronoi(offset,s1)
# Generate Part-IniRock
myModel.convertAllSketches()
p = myModel.Part(name="IniRock", dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
p.BaseShell(sketch=sr)
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s0)
# Generate Part-JointRock
myModel.convertAllSketches()
p = myModel.Part(name="JointRock", dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
p.BaseShell(sketch=sr)
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s1)
mdb.saveAs(pathName='column_jointed_rock.cae')
'''
Part-3: Generate the numerical model in PFC and UDEC
For PFC model, two files is generate,one for the generation of homogeneous model 
and another for generate the joint as fracture
'''
partname='IniRock'
filename='Joint_rock.p2dat'
Joint_to_PFC(mdb,partname,filename,domain)
print 'Write file to PFC is OK'
partname='IniRock'
filename='Joint_rock.uddat'
print 'Write file to UDEC is OK'
Joint_to_UDEC(mdb,partname,filename,domain)
print 'Finish, have a nice day!'