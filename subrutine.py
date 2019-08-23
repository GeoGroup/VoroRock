# -*- coding: utf-8 -*-  

from math import sqrt
import random
#lines are stored in (point, slope) form, where slope==None for vertical
#polygons are stored as a clockwise list of vertices
#all polygons are assumed to be convex
def bounded_voronoi(bounding_polygon, point_set):
    """Returns a [point]->[polygon] hash, where:
       -the point is from point_set,
       -polygon is the points of point's voronoi polygon
    All points in point_set should be contained in bounding_polygon.
    """
    polygons = {}
    for point in point_set:
        bound = bounding_polygon
        for other_point in point_set:
            if other_point == point:
                continue
            hp = halfplane(point, other_point)
            bound = polygon_intersect_halfplane(bound, hp)
        polygons[point] = bound
    return polygons

def update_diagram(diagram, new_point, bounding_polygon):
    point_set = diagram.keys()
    changed = []
    for point in point_set:
        bound = diagram[point]
        hp = halfplane(point, new_point)
        new_bound = polygon_intersect_halfplane(bound, hp)
        if bound != new_bound:
            diagram[point] = new_bound
            changed.append(point)
    bound = bounding_polygon
    for other_point in point_set:
        hp = halfplane(new_point, other_point)
        bound = polygon_intersect_halfplane(bound, hp)
    diagram[new_point] = bound

def halfplane(inner_point, outer_point):
    """Returns a pair (line, {'lt'|'gt'}), where 'lt' means that inner_point is 
    below the line (or to the left for a vertical line), and the line is the 
    perpendicular bisector of the segment defined by the two given points.
    """
    mp = midpoint(inner_point, outer_point)
    s = slope(inner_point, outer_point)
    if s == None:
        new_slope = 0
    elif s == 0:
        new_slope = None
    else:
        new_slope = -(1/s)
    perpendicular_bisector = (mp, new_slope) 
    if new_slope == None:
        if inner_point[0] < mp[0]:
            half = 'lt'
        else:
            half = 'gt'
        return (perpendicular_bisector, half)
    inner_b = inner_point[1]-inner_point[0]*new_slope
    real_b = mp[1]-mp[0]*new_slope
    if inner_b < real_b:
        half = 'lt'
    else:
        half = 'gt'
    return (perpendicular_bisector, half)

def polygon_intersect_halfplane(polygon, halfplane):
    """Returns a new polygon which is the intersection of the two arguments.
    This is not general use: it assumes that there is a nonempty intersection.
    """
    n = len(polygon)
    point_locations = []
    #print "polygon_intersect_halfplane: polygon", polygon
    #print "polygon_intersect_halfplane: halfplane", halfplane
    for point in polygon:
        point_locations.append(point_in_halfplane(point, halfplane))
    if len(filter(bool, point_locations)) == n:
        return polygon
    while not (point_locations[0] and not point_locations[-1]):
        #print point_locations
        #print polygon
        point_locations = lshift_list(point_locations)
        polygon = lshift_list(polygon)
    #print "polygon_intersect_halfplane: polygon", polygon
    #print "polygon_intersect_halfplane: point_locations", point_locations
    new_polygon = []
    intersects = 0
    for p in range(n):
        #print "P=", p
        if point_locations[p]:
            new_polygon.append(polygon[p])
        elif intersects == 0:
            intersects = 1
            if point_locations[p] == None:
                new_polygon.append(polygon[p])
            else:
                segment = (polygon[p-1], polygon[p])
                line = halfplane[0]
                new_point = segment_intersect_line(segment, line)
                new_polygon.append(new_point)
        if p == n-1:
            #print "On the last point of the polygon"
            if point_locations[p] == None:
                #print " And it's a None"
                new_point = polygon[p]
                if new_point != new_polygon[-1]:
                    new_polygon.append(new_point)
                    #print "  And it was appended"
            else:
                #print " And it isn't a None"
                segment = (polygon[p], polygon[0])
                line = halfplane[0]
                new_point = segment_intersect_line(segment, line)
                if new_point == None:
                    print "ERROR"
                    print "Found no intersection between segment and line"
                    print segment, line
                new_polygon.append(new_point)
    #print "----returning from polygon_intersect halfplane---"
    return new_polygon

def lshift_list(l):
    end = l[0]
    l = l[1:]
    l.append(end)
    return l

def point_in_halfplane(point, halfplane):
    #print "point_in_halfplane: point, halfplane", point, halfplane
    #special case: vertical line
    if halfplane[0][1] == None:
        x = halfplane[0][0][0]
        if halfplane[1] == 'lt':
            if point[0] < x:
                return True
            elif point[0] == x:
                return None
            else:
                return False
        else:
            if point[0] > x:
                return True
            elif point[0] == x:
                return None
            else:
                return False
    #otherwise: not vertical...
    real_b = halfplane[0][0][1]-halfplane[0][0][0]*halfplane[0][1]
    point_b = point[1]-point[0]*halfplane[0][1]
    if halfplane[1] == 'lt':
        if point_b < real_b:
            return True
    else:
        if point_b > real_b:
            return True
    if real_b == point_b:
        return None
    return False

def segment_intersect_line(segment, line):
    #returns either a point (x, y) or None; endpoints not included
    seg_line = (segment[0], slope(segment[0], segment[1]))
    intersection = line_intersect(seg_line, line)
    if intersection == None:
        return None
    if point_on_segment(intersection, segment):
        return intersection
    else:
        return None

def line_intersect(l1, l2):
    # test le coefficient directeur
    if l1[1] == l2[1]:  
        return None
    ((x1, y1), m1) = l1
    ((x2, y2), m2) = l2
    if m1 == None:
        x = x1
        y = m2*(x-x2)+y2
        return (x, y)
    if m2 == None:
        x = x2
        y = m1*(x-x1)+y1
        return (x, y)
    x = (m1*x1-y1-m2*x2+y2)/(m1-m2)
    y = m1*(x-x1)+y1
    return (x, y)

def point_on_segment(point, segment):
    #IMPORTANT: the point is assumed to be on the line defined by the segment
    d1 = dist(segment[0], segment[1])
    d2 = dist(point, segment[0])
    d3 = dist(point, segment[1])
    return ((d3 <= d1) and (d2 <= d1))

def dist(p1, p2):
    return sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

def midpoint(p1, p2):
    return ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2)

def normalize(v):
    length = sqrt(v[0]**2 + v[1]**2)
    v=(v[0]/length, v[1]/length)
    return v

def slope(p1, p2):
    if (p2[0] == p1[0]):
        #special case: vertical line
        return None
    return (p2[1]-p1[1])/(p2[0]-p1[0])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def point_on_bound(pt,bound):
   for i in range(0,len(bound)-1):
      seg1=(bound[i],bound[i+1])
      a1=slope(pt,bound[i+1])
      a2=slope(bound[i],pt)
      if(a1==a2 and point_on_segment(pt,seg1)):
         return True
   seg1=(bound[len(bound)-1],bound[0])
   a1=slope(pt,bound[len(bound)-1])
   a2=slope(bound[0],pt)
   if(a1==a2 and point_on_segment(pt,seg1)):
      return True
   return False

def drawVoronoi(voronoi,s):
   points=voronoi.keys()
   for i in range(0,len(points)):#len(points)
     #if (i%100==0):
     #  print 'Voro-num '+str(i)+'/'+str(len(points))
     poly = voronoi[points[i]] 
     if(len(poly)>0):
       for j in range(0,len(poly)-1):
         s.Line(poly[j+1],poly[j])
       s.Line(poly[len(poly)-1],poly[0])  


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pointsInBound(point,bound):
  if((point[0]>=bound[0][0]) and (point[0]<=bound[2][0]) and (point[1]>=bound[0][1]) and (point[1]<=bound[1][1])):
    return True
  else:
    return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def lineInBound(line,bound):
   point1=line[0]
   point2=line[1]  
   if(pointsInBound(point1,bound)): return True
   if(pointsInBound(point2,bound)): return True
   ml1=slope(point1,point2) 
   l1=(point1, ml1)
   for i in (0,3):
      if(i<3):
        seg=(bound[i],bound[i+1])
      else:
        seg=(bound[i],bound[0]) 
      pt3=segment_intersect_line(seg,l1)
      if(pt3!= None):
        if point_on_segment(pt3, line):
           return True
   return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def removePoints(voronoi,points,thepoint,eps):
  for i in range(0,len(points)):   
    poly = voronoi[points[i]]    
    npoly=len(poly)
    j=0
    while (j<npoly):
      d=dist(poly[j],thepoint)
      if(d<eps):
         voronoi[points[i]][j]=thepoint
      j=j+1
   # suppression des redondances
  for i in range(0,len(points)):   
    poly = voronoi[points[i]]    
    npoly=len(poly)
    j=0
    while (j<npoly):
      d=dist(poly[j],thepoint)
      if(j==len(poly)-1): 
        d=dist(poly[j],poly[0])
      else: 
         d=dist(poly[j+1],poly[j])
      if(d<1.e-10):  
         voronoi[points[i]].remove(poly[j])    
         npoly=npoly-1
      j=j+1
  return voronoi   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def verifVoronoi(voronoi,points,eps):
   for i in range(0,len(points)):   
     poly = voronoi[points[i]]    
     npoly=len(poly)
     j=0
     while (j<npoly):
       if(j==npoly-1): 
         d=dist(poly[j],poly[0])
       else: 
         d=dist(poly[j+1],poly[j])
       if(d<eps):
         thepoint=poly[j]
         voronoi[points[i]].remove(poly[j])   
         voronoi=removePoints(voronoi,points,thepoint,eps)   
         poly = voronoi[points[i]]      
         npoly=len(poly)
       j=j+1
   return voronoi 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pointsinlist(pt,liste):
  for i in range(0,len(liste)):
      if(dist(pt,liste[i])<1.e-4): return True
  return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# search the segment sharing the prev node in a list
def findnextseg(prev,seglist):
   for i in range(0,len(seglist)):
      seg=seglist[i]
      if(dist(prev,seg[0])<1.e-4):
        next=seg[1]   
        return True,next,i
   for i in range(0,len(seglist)):
      seg=seglist[i]
      if(dist(prev,seg[1])<1.e-4):
        next=seg[0]   
        return True,next,i
   return False,None,0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check the self intersection of polygon
def check_poly(newpointj):
 newpointj.append(newpointj[0])
 newpointj.append(newpointj[1])
 lg2=[]
 clock=[]
 newpointf=[]
 for j in range(len(newpointj)-2):
  pt1=midpoint(newpointj[j],newpointj[j+1])
  ml1=slope(newpointj[j],newpointj[j+1])
  lig=(pt1,ml1)
  lg2.append(lig)
  vect1=(newpointj[j+1][0]-newpointj[j][0], newpointj[j+1][1]-newpointj[j][1])
  vect2=(newpointj[j+2][0]-newpointj[j+1][0], newpointj[j+2][1]-newpointj[j+1][1])
  clock.append(vect1[0]*vect2[1]-vect1[1]*vect2[0])
 lg2.append(lg2[0])
 if clock[0]<0:
  newpointf.append(newpointj[1])
  for j in range(1,len(clock)):
   if clock[j]>0:
    ap=line_intersect(lg2[j+0], lg2[j+2])
    newpointf.append(ap)
    clock[j+1]=0
   elif clock[j]<0:
    newpointf.append(newpointj[j+1])
 if clock[0]>0:
  if clock[1]>0:
   ap=line_intersect(lg2[0], lg2[2])
   newpointf.append(ap)
   for j in range(2,len(clock)):
    if clock[j]>0:
     ap=line_intersect(lg2[j+0], lg2[j+2])
     newpointf.append(ap)
     clock[j+1]=0
    elif clock[j]<0:
     newpointf.append(newpointj[j+1])
  if clock[1]<0:
   ap=line_intersect(lg2[0], lg2[len(lg2)-2])
   newpointf.append(ap)
   for j in range(1,len(clock)):
    if clock[j]>0:
     ap=line_intersect(lg2[j+0], lg2[j+2])
     newpointf.append(ap)
     clock[j+1]=0
    elif clock[j]<0:
     newpointf.append(newpointj[j+1])
 del newpointj[len(newpointj)-1]
 return newpointf

# Remove the short edges less than error
def regulization(voronoi, bound, error):
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
            center=check_bound(center,bound,0.001)
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

def check_bound(point,bound,error):
	# check the point at the boundary
	if abs(point[0]-bound[0][0])<error:
		point=(bound[0][0],point[1])
	elif abs(point[0]-bound[2][0])<error:
		point=(bound[2][0],point[1])
	
	if abs(point[1]-bound[0][1])<error:
		point=(point[0],bound[0][1])
	elif abs(point[1]-bound[2][1])<error:
		point=(point[0],bound[2][1])
	
	return point

# generate n random points in domain
def generate_random_seed(domain,Jd):
	points = []
	num=int(Jd*domain[0]*domain[1])
	for i in range(num):
		newpoint = (random.uniform(0,domain[0]), random.uniform(0,domain[1]))
		points.append(newpoint)
	return points


# generate 9*n points in 3R*3R domain
def generate_new_seed(domain,points):
	point2=[]
	for i in range(0,len(points)):
		newpoint=points[i]
		for ix in (-1,1,0):
			for iy in (1,-1,0): 
				pointadd=(newpoint[0]+float(ix)*domain[0],newpoint[1]+float(iy)*domain[1])
				point2.append(pointadd)
	return point2

# generate new 3R*3R domain
def generate_new_domain(domain):
	bound = [(0., 0.), (0., domain[1]), (domain[0], domain[1]), (domain[0], 0.)]
	domain2 = [(-bound[2][0],-bound[2][1]), (-bound[2][0],2.*bound[2][1]), (2.*bound[2][0],2.*bound[2][1]), (2.*bound[2][0],-bound[2][1])]
	return domain2

# get the centroid and area of a polygon
def get_centroid_area(poly):
	poly.append(poly[0])
	s=0;gx=0;gy=0;tmp=0
	for i in range(len(poly)-1):
		tmp=0.5*(poly[i][0]*poly[i+1][1]-poly[i+1][0]*poly[i][1])
		#tmp=abs(tmp)
		gx=gx+tmp*(poly[i][0]+poly[i+1][0])/3
		gy=gy+tmp*(poly[i][1]+poly[i+1][1])/3
		s=s+tmp
	gx=gx/s
	gy=gy/s	
	del poly[len(poly)-1]
	return (gx,gy),abs(s)


# get the centroid and coefficient of varaiation of voronoi diagram
def get_cv_ct(voronoi,points):
	import numpy as np 
	Clist=[]
	Slist=[]
	for i in range(len(points)):
		poly = voronoi[points[i]] 
		[centre,s]=get_centroid_area(poly)
		Clist.append(centre)
		Slist.append(s)
	arr_std = np.std(Slist,ddof=1)
	arr_mean = np.mean(Slist)
	cv=arr_std/arr_mean
	return cv,Clist

# add for two lists
def list_add_mean(a,b):
	c = []
	for i in range(len(a)):
		c.append((a[i][0]/2+b[i][0]/2,a[i][1]/2+b[i][1]/2))
	return c

# modified Lloyd's iteration for a specified CV
def iteration_lloyd_cv(bound,points,CV,error):
	rcv=[]
	pt=[]
	voronoi = bounded_voronoi(bound, points)
	[tcv,ct]=get_cv_ct(voronoi,points)
	rcv.append(tcv)
	iter_num=1
	while abs(tcv-CV)>error:
		if tcv>CV:
			pt=points
			points=ct
		else:
			points=list_add_mean(points,pt)
		voronoi = bounded_voronoi(bound, points)
		[tcv,ct]=get_cv_ct(voronoi,points)
		rcv.append(tcv)
		iter_num=iter_num+1
	return points,rcv,iter_num

# Add joint thickness using polygon offset
def polygon_offset(voronoi2,point2,nbound2,Jt):
	offset={}
	for i in range(0,len(point2)):
		point=point2[i]
		#if (i%20==0):
		#	print 'Grain-num '+str(i)+'/'+str(len(point2))
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
			ptl=(pointmil[0]+scalar*(Jt/2.0)*vect[1],pointmil[1]-scalar*(Jt/2.0)*vect[0] )
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
		ptl=(pointmil[0]-scalar*(Jt/2.0)*vect[1], pointmil[1]+scalar*(Jt/2.0)*vect[0])
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
		offset[point]=newpointf
	return offset
# Add joint thickness using polygon offset
def polygon_offset2(voronoi2,point2,nbound2,Jt):
	offset={}
	for i in range(0,len(point2)):
		point=point2[i]
		if (abs(point[0]+nbound2[0][0]/2)<abs(nbound2[0][0]/2)*1.25) and (abs(point[1]+nbound2[0][1]/2)<abs(nbound2[0][1]/2)*1.25):
			#if (i%20==0):
			#	print 'Grain-num '+str(i)+'/'+str(len(point2))
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
				ptl=(pointmil[0]+scalar*(Jt/2.0)*vect[1],pointmil[1]-scalar*(Jt/2.0)*vect[0] )
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
			ptl=(pointmil[0]-scalar*(Jt/2.0)*vect[1], pointmil[1]+scalar*(Jt/2.0)*vect[0])
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
			offset[point]=newpointf
	return offset

# Export to the PFC version 6.0
# Get the direction and length of a line
def direct_scale(line):
  from math import acos
  from math import pi
  if line[1]>line[3]:
    line=[line[2],line[3],line[0],line[1]];
  direct=[line[2]-line[0],line[3]-line[1]];
  scale=sqrt(direct[0]**2+direct[1]**2);
  dip=(-direct[0])/scale;
  dip=acos(dip)*180/pi;
  return dip,scale

# Export the joint information to PFC model
def Joint_to_PFC(mdb,partname,filename,domain=[1,1]):
  f=open(filename,'w+')
  f.write('model new\n')
  f.write('model title \"Jointed Rock Simulation\"\n')
  f.write('; Set the domain extent\n')
  f.write('model domain extent %4.2f %4.2f %4.2f %4.3f\n' % (-0.25*domain[0],1.25*domain[0],-0.25*domain[1],1.25*domain[1]))
  f.write('contact cmat default model linearpbond property kn 5e7 dp_nratio 0.5 \n')
  f.write('; Generate walls \n')
  f.write('wall generate box %4.2f %4.2f %4.2f %4.2f\n' % (0,domain[0],0,domain[1]))
  f.write('; Distribute balls in the box.\n')
  f.write('model random 1001\n')
  f.write('ball distribute porosity 0.15 radius %f %f box  %4.2f %4.2f %4.2f %4.3f\n' % (domain[0]/125.0,domain[1]/100.0,0,domain[0],0,domain[1]))
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
  p=mdb.models['Model-1'].parts[partname]
  v=p.vertices
  e=p.edges
  e0=e[0:0];e1=e0;
  for i in range(len(e)):
    ix=0
    et=e[i]
    vt=et.pointOn
    if abs(vt[0][0])<1e-3 or abs(vt[0][1])<1e-3 or abs(vt[0][0]-domain[0])<1e-3 or abs(vt[0][1]-domain[1])<1e-3:
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
    f.write('         dfn  \"Joint\" \n')
  #f.write('program call \'Fracture\'\n')
  f.write('fracture property ...\n')
  f.write('         \'sj_kn\' 2e9 \'sj_ks\' 2e9 \'sj_fric\' 0.70 ...\n')
  f.write('         \'sj_coh\' 0.0 \'sj_ten\' 0.0 \'sj_large\' 1\n')
  f.write('; Apply smoothjoint contact model to contacts intercepted by fracture\n')
  f.write('fracture contact-model model \'smoothjoint\' install \n')
  f.write('; Ensure new contacts intersecting the fracture are set to the sj contact model \n')
  f.write('fracture contact-model model \'smoothjoint\' activate \n')
  f.write('model save \'Joint_Rock\'\n')
  f.close()


# Export joint information to UDEC
def Joint_to_UDEC(mdb,partname,filename,domain=[1,1]):
  f=open(filename,'w+')
  f.write('new\n')
  f.write('round 0.001\n')  
  f.write('block  0 0 0 %4.2f %4.2f %4.2f %4.3f 0\n' % (domain[1],domain[0],domain[1],domain[1]))
  p=mdb.models['Model-1'].parts[partname]
  v=p.vertices
  e=p.edges
  e0=e[0:0];e1=e0;
  for i in range(len(e)):
    ix=0
    et=e[i]
    vt=et.pointOn
    if abs(vt[0][0])<1e-3 or abs(vt[0][1])<1e-3 or abs(vt[0][0]-domain[0])<1e-3 or abs(vt[0][1]-domain[1])<1e-3:
      e0=e0+e[i:i+1]
    else:
      e1=e1+e[i:i+1]
  for i in range(len(e1)):
    vt=e1[i].getVertices()
    vl=[v[vt[0]].pointOn[0][0],v[vt[0]].pointOn[0][1],v[vt[1]].pointOn[0][0],v[vt[1]].pointOn[0][1]]
    f.write('crack (%f, %f) (%f, %f)\n' % (vl[0],vl[1],vl[2],vl[3]))
  f.close()
