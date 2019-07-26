# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2017 replay file
# Internal Version: 2016_09_28-05.54.59 126836
# Run by MQX on Sat Mar 23 10:54:58 2019
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=356.101165771484, 
    height=253.399993896484)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
openMdb('GBM.cae')
#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
#* MdbError: incompatible release number, expected 2017, got 6.14-5
upgradeMdb("C:/Users/m-q-x/Desktop/OpenGBM/GBM-6.14-5.cae", 
    "C:/Users/m-q-x/Desktop/OpenGBM/GBM.cae", )
#: The model database "C:\Users\m-q-x\Desktop\OpenGBM\GBM_TEMP.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
#: The model database has been saved to "C:\Users\m-q-x\Desktop\OpenGBM\GBM.cae".
#: The model database "C:\Users\m-q-x\Desktop\OpenGBM\GBM-6.14-5.cae" has been converted.
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
a = mdb.models['Model-1'].rootAssembly
a.unlock()
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['GBM']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].Part(name='GBM-failed', 
    objectToCopy=mdb.models['Model-1'].parts['GBM'])
mdb.models['Model-1'].parts['GBM-failed'].Unlock(reportWarnings=False)
del mdb.models['Model-1'].parts['GBM']
mdb.models['Model-1'].parts.changeKey(fromName='GBM-failed', toName='GBM')
import assembly
mdb.models['Model-1'].rootAssembly.regenerate()
p1 = mdb.models['Model-1'].parts['GBM']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['GBM']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap=session.viewports['Viewport: 1'].colorMappings['Set']
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
mdb.models['Model-1'].parts['GBM'].deleteSets(setNames=('set-1', 'set-2', 
    'set-3', 'set-4', 'set-5', 'set-6', 'set-7', 'set-8', 'set-9', 'set-10', 
    'set-11', 'set-12', 'set-13', 'set-14', 'set-15', 'set-16', 'set-17', 
    'set-18', 'set-19', 'set-20', 'set-21', 'set-22', 'set-23', 'set-24', 
    'set-25', 'set-26', 'set-27', 'set-28', 'set-29', 'set-30', 'set-31', 
    'set-32', 'set-33', 'set-34', ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.67164, 
    farPlane=2.98522, width=1.48948, height=1.10175, viewOffsetX=0.0545431, 
    viewOffsetY=0.039957)
p = mdb.models['Model-1'].parts['GBM']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#0 #2 ]', ), )
p.Set(faces=faces, name='rock')
#: The set 'rock' has been created (1 face).
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.69819, 
    farPlane=2.95866, width=1.104, height=0.816617, viewOffsetX=-0.0646386, 
    viewOffsetY=-0.0744476)
p1 = mdb.models['Model-1'].parts['GBM']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.7425, 
    farPlane=2.91436, width=0.727674, height=0.538254, viewOffsetX=-0.196419, 
    viewOffsetY=-0.206999)
p = mdb.models['Model-1'].parts['GBM']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1000 #2 ]', ), )
p.Set(faces=faces, name='rock')
#: The set 'rock' has been edited (2 faces).
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.61947, 
    farPlane=3.03739, width=1.98989, height=1.47191, viewOffsetX=0.274302, 
    viewOffsetY=-0.100408)
mdb.models['Model-1'].parts['GBM'].sets.changeKey(fromName='rock', 
    toName='joint')
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.44695, 
    farPlane=3.2099, width=3.69181, height=2.7308, viewOffsetX=-0.341711, 
    viewOffsetY=-0.123747)
p = mdb.models['Model-1'].parts['GBM']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#ffffefff #1 ]', ), )
e = p.edges
edges = e.getSequenceFromMask(mask=('[#0 #50000000 #0:2 #ff000000 #fff ]', ), )
p.Set(edges=edges, faces=faces, name='rock')
#: The set 'rock' has been created (32 faces, 22 edges).
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.67993, 
    farPlane=2.97692, width=1.4159, height=1.04733, viewOffsetX=0.00531113, 
    viewOffsetY=-0.0587586)
del mdb.models['Model-1'].parts['GBM'].sets['rock']
del mdb.models['Model-1'].parts['GBM'].sets['joint']
session.viewports['Viewport: 1'].view.setValues(width=1.50123, height=1.11045, 
    viewOffsetX=0.200653, viewOffsetY=-0.161195)
p = mdb.models['Model-1'].parts['GBM']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#ffffefff #1 ]', ), )
p.Set(faces=faces, name='rock')
#: The set 'rock' has been created (32 faces).
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.78897, 
    farPlane=2.86788, width=0.33375, height=0.246872, viewOffsetX=-0.324952, 
    viewOffsetY=-0.363604)
p = mdb.models['Model-1'].parts['GBM']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1000 #2 ]', ), )
p.Set(faces=faces, name='joint')
#: The set 'joint' has been created (2 faces).
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.62426, 
    farPlane=3.0326, width=1.94508, height=1.43876, viewOffsetX=-0.128667, 
    viewOffsetY=-0.177638)
session.viewports['Viewport: 1'].view.fitView()
session.printOptions.setValues(reduceColors=False)
session.printToFile(fileName='hetegbm', format=PNG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)
s1.rectangle(point1=(-10.0, -10.0), point2=(10.0, 10.0))
p = mdb.models['Model-1'].Part(name='Part-3', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-3']
p.BaseShell(sketch=s1)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-3']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-3']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.Set(faces=faces, name='homo')
#: The set 'homo' has been created (1 face).
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap = session.viewports['Viewport: 1'].colorMappings['Set']
cmap.updateOverrides(overrides={'homo':(True, '#FFD700', 'Default', 
    '#FFD700')})
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap = session.viewports['Viewport: 1'].colorMappings['Set']
cmap.updateOverrides(overrides={'homo':(True, '#808000', 'Default', 
    '#808000')})
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap = session.viewports['Viewport: 1'].colorMappings['Set']
cmap.updateOverrides(overrides={'homo':(True, '#008080', 'Default', 
    '#008080')})
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
session.printToFile(fileName='homo', format=PNG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
