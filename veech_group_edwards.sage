############################################################
############################################################
##  Classes and functions used within veech_group() below ##
############################################################
############################################################

# Will print progress of various for loops
import sys
def print_percent_complete(msg, i, n_i):
    i, n_i = int(i)+1, int(n_i)
    sys.stdout.write('\r')
    sys.stdout.write("{} {}/{}" .format(msg, i, n_i)) #, (100/(n_i)*i).n()))
    sys.stdout.flush()
    if i/n_i == 1:
        print()


# A geodesic in the upper-half plane
class geodesic:
    # We can define a geodesic in one of three ways:
    # (i) Provide a 2x2 matrix M, in which case self is the perpendicular bisector of the geodesic connecting I and M*I in the upper-half plane
    # (ii) Provide one real number, foot, in which case self is the vertical geodesic with real part foot
    # (iii) Provide two real numbers, foot_left<foot_right, in which case self is the semicircular geodesic with these two feet
    def __init__(self,matrix=None,foot=None,foot_left=None,foot_right=None):
        self.is_polygon_side=False
        self.vertex_counterclockwise=None
        self.vertex_clockwise=None
        # Case (i)
        if matrix!=None:
            self.matrix=matrix
            # Here self can be one of three types:
            # (a) Vertical: A vertical line
            # (b) Encloses I: A semicircle centered on the real axis that encloses I
            # (c) Exposes I: A semicircle centered on the real axis that does not enclose I
            # We find what type self is and its feet (note that a vertical geodesic has only one finite foot while the other types have two distinct feet)
            # We also, by default, determine the counterclockwise and clockwise vertices corresponding to self depending on its feet; if self is found to be an edge of a hyperbolic polygon, these vertices may change to represent the vertices on the edge of the polygon corresponding to self
            M=matrix
            M_act_on_I=(QQbar(M[0][0])*QQbar(I)+QQbar(M[0][1]))/(QQbar(M[1][0])*QQbar(I)+QQbar(M[1][1]))
            x=M_act_on_I.real()
            y=M_act_on_I.imag()
            # Case (a)
            if y==1:
                self.type='Vertical'
                self.foot=x/2
            else:
                # If the real part of M*I is 0, then self is a semi-circle with feet at -y^(1/2) and y^(1/2)
                if x==0:
                    self.foot_left=-y^(1/2)
                    self.foot_right=y^(1/2)
                else: 
                    # The center and radius, respectively, of the semicircular geodesic connecting I and M_act_on_I
                    c=(x^2+y^2-1)/(2*x)
                    R=(1+c^2)^(1/2)
                    # Left foot of h-line through I and M_act_on_I
                    l=c-R
                    # Imaginary part of hyperbolic midpoint between $I$ and $phi(M_act_on_I)$, where $phi$ is the Mobius transformation taking the semicircular geodesic connecting I and M_act_on_I with left and right feet l and r, resp., to imaginary axis where l -> 0, i -> i, and r -> infinity
                    # The map $phi$ is given by $phi(z)=(z-l)/(l*z+1)$, and $phi^{-1}(z)=(z+l)/(-l*z+1)$
                    # imag_mdpt is the squareroot of the imaginary part of $phi(M_act_on_I)$
                    imag_mdpt=((y*(1+l^2))/(l^2*(x^2+y^2)+2*l*x+1))^(1/2)
                    # The left and right feet of the perpendicular bisector of the geodesic between $I$ and $phi(M_act_on_I)$ are -imag_mdpt and imag_mdpt
                    # The feet of the perpendicular bisector of the geodesic between I and M_act_on_I are given by $phi^{-1}(-imag_mdpt)$ and $phi^{-1}(imag_mdpt)$
                    foot0=(-imag_mdpt+l)/(l*imag_mdpt+1)
                    foot1=(imag_mdpt+l)/(-l*imag_mdpt+1)
                    if foot0<foot1:
                        self.foot_left=foot0
                        self.foot_right=foot1
                    else:
                        self.foot_left=foot1
                        self.foot_right=foot0
                # Case (b)
                if self.foot_left*self.foot_right<-1:
                    self.type='Encloses I'
                # Case (c)  
                else:
                    self.type='Exposes I'
                # Center and radius of the semicircle describing self
                self.center=(self.foot_right+self.foot_left)/2
                self.radius=(self.foot_right-self.foot_left)/2

            # Hyperbolic distance from I to the perpendicular bisector (as the infimum of distances from I to any point of the perpendicular bisector)
            frobenius_norm=sqrt(M[0][0]^2+M[0][1]^2+M[1][0]^2+M[1][1]^2)
            self.distance_from_I=log(nu(frobenius_norm))
        # Case (ii)
        elif foot!=None:
            self.type='Vertical'
            self.foot=foot
        # Case (iii)        
        elif foot_left!=None and foot_right!=None:
            self.foot_left=foot_left
            self.foot_right=foot_right
            if self.foot_left*self.foot_right<-1:
                self.type='Encloses I'
            else:
                self.type='Exposes I'
            # Center and radius of semicirlce describing self
            self.center=(self.foot_right+self.foot_left)/2
            self.radius=(self.foot_right-self.foot_left)/2

    def __repr__(self):
        if self.type=='Vertical':
            return 'Vertical hyperbolic geodesic with foot at {}'.format(self.foot)
        elif self.type=='Encloses I':
            return 'Hyperbolic geodesic enclosing I with feet at {} and {}'.format(self.foot_left,self.foot_right)
        elif self.type=='Exposes I':
            return 'Hyperbolic geodesic exposing I with feet at {} and {}'.format(self.foot_left,self.foot_right)

    # Option to plot self
    def plot(self,x_min=-5,x_max=5,y_min=0,y_max=5):
        if self.type=='Vertical':
            x=var('x')
            y=var('y')
            return implicit_plot(x-self.foot,(x,x_min,x_max),(y,y_min,y_max),color='black')
        else:
            return arc((self.center,0),self.radius,sector=(0,pi),color='black')

    # Here side is a geodesic (different from self) corresponding to an edge of a hyperbolic polygon (side.is_polygon_side=True) with vertices side.vertex_counterclockwise and side.vertex_clockwise
    # If self and side intersect, then we test if the intersection is in the interior of the segment of side between side.vertex_counterclockwise and side.vertex_clockwise; if so, we update the vertices of self and side and set self.is_polygon_side=True
    def intersect_polygon_side(self,side,timeout=1,prec=53):
        # If both self and side are vertical, then they cannot intersect
        if self.type=='Vertical':
            if side.type!='Vertical':
                # Check if self and side intersect at all
                if strictly_less_than(self.foot,side.foot_right,timeout=timeout,prec=prec) and strictly_greater_than(self.foot,side.foot_left,timeout=timeout,prec=prec):
                    # Real and imaginary part of intersection
                    x=self.foot
                    y=sqrt(side.radius^2-(x-side.center)^2)
                    # Check if the intersection is between the vertices of the polygon edge that side corresponds to
                    if strictly_greater_than(x,min(side.vertex_counterclockwise[0],side.vertex_clockwise[0]),timeout=timeout,prec=prec) and strictly_less_than(x,max(side.vertex_counterclockwise[0],side.vertex_clockwise[0]),timeout=timeout,prec=prec):
                        self.is_polygon_side=True
                        # Update the vertices of both self and side corresponding to this intersection
                        if (side.type=='Encloses I' and self.foot<0) or (side.type=='Exposes I' and self.foot>0):
                            self.vertex_clockwise=[x,y]
                            side.vertex_counterclockwise=[x,y]
                        elif (side.type=='Encloses I' and self.foot>0) or (side.type=='Exposes I' and self.foot<0):
                            self.vertex_counterclockwise=[x,y]
                            side.vertex_clockwise=[x,y]
        # self is not vertical
        else:
            # side is vertical
            if side.type=='Vertical':
                # Check if self and side intersect at all
                if strictly_greater_than(self.foot_right,side.foot,timeout=timeout,prec=prec) and strictly_less_than(self.foot_left,side.foot,timeout=timeout,prec=prec):
                    # Real and imaginary part of intersection
                    x=side.foot
                    y=sqrt(self.radius^2-(x-self.center)^2)
                    # Check if the intersection is between the vertices of the polygon edge that side corresponds to
                    if strictly_greater_than(y,min(side.vertex_counterclockwise[1],side.vertex_clockwise[1]),timeout=timeout,prec=prec) and strictly_less_than(y,max(side.vertex_counterclockwise[1],side.vertex_clockwise[1]),timeout=timeout,prec=prec):
                        self.is_polygon_side=True
                        # Update the vertices of both self and side corresponding to this intersection
                        if (self.type=='Enloses I' and side.foot<0) or (self.type=='Exposes I' and side.foot>0):
                                self.vertex_counterclockwise=[x,y]
                                side.vertex_clockwise=[x,y]
                        elif (self.type=='Enloses I' and side.foot>0) or (self.type=='Exposes I' and side.foot<0):
                                self.vertex_clockwise=[x,y]
                                side.vertex_counterclockwise=[x,y]  
            # side is not vertical
            else:
                # Check if self and side intersect at all (this can happen in one of two ways):
                # (i) Check if the feet of self straddle the left foot of side
                if strictly_less_than(self.foot_right,side.foot_right,timeout=timeout,prec=prec) and strictly_greater_than(self.foot_right,side.foot_left,timeout=timeout,prec=prec) and strictly_less_than(self.foot_left,side.foot_left,timeout=timeout,prec=prec):
                    # Real and imaginary part of intersection
                    x=(self.radius^2-self.center^2-(side.radius^2-side.center^2))/(2*(side.center-self.center))
                    y=sqrt(self.radius^2-(x-self.center)^2)
                    # Check if the intersection is between the vertices of the polygon edge that side corresponds to
                    if strictly_greater_than(x,min(side.vertex_counterclockwise[0],side.vertex_clockwise[0]),timeout=timeout,prec=prec) and strictly_less_than(x,max(side.vertex_counterclockwise[0],side.vertex_clockwise[0]),timeout=timeout,prec=prec):
                        self.is_polygon_side=True
                        # Update the vertices of both self and side corresponding to this intersection
                        if self.type==side.type:
                            self.vertex_counterclockwise=[x,y]
                            side.vertex_clockwise=[x,y]
                        else:
                            self.vertex_clockwise=[x,y]
                            side.vertex_counterclockwise=[x,y]
                # (ii) Otherwise, check if the feet of self straddle the right foot of side                 
                elif strictly_greater_than(self.foot_left,side.foot_left,timeout=timeout,prec=prec) and strictly_less_than(self.foot_left,side.foot_right,timeout=timeout,prec=prec) and strictly_greater_than(self.foot_right,side.foot_right,timeout=timeout,prec=prec):
                    # Real and imaginary part of intersection
                    x=(self.radius^2-self.center^2-(side.radius^2-side.center^2))/(2*(side.center-self.center))
                    y=sqrt(self.radius^2-(x-self.center)^2)
                    # Check if the intersection is between the vertices of the polygon edge that side corresponds to
                    if strictly_greater_than(x,min(side.vertex_counterclockwise[0],side.vertex_clockwise[0]),timeout=timeout,prec=prec) and strictly_less_than(x,max(side.vertex_counterclockwise[0],side.vertex_clockwise[0]),timeout=timeout,prec=prec):
                        self.is_polygon_side=True
                        # Update the vertices of both self and side corresponding to this intersection
                        if self.type==side.type:
                            self.vertex_clockwise=[x,y]
                            side.vertex_counterclockwise=[x,y]
                        else:
                            self.vertex_counterclockwise=[x,y]
                            side.vertex_clockwise=[x,y]                         


# Various comparisons throughout the code get hung up on exact or approximate computations depending on the surface; these allow for an exact comparison to run for at most a specified time (timeout), then switch to a specified approximate comparison (prec bits of precision) if the time limit is reached
def is_strictly_less_than(a,b,timeout=1):
    @fork(timeout=timeout)
    def less_than(c,d):
        return bool(c<d)
    return less_than(a,b)
def is_strictly_greater_than(a,b,timeout=1):
    @fork(timeout=timeout)
    def gr_than(c,d):
        return bool(c>d)
    return gr_than(a,b)    
def is_less_than_or_equal(a,b,timeout=1):
    @fork(timeout=timeout)
    def less_or_eq(c,d):
        return bool(c<=d)
    return less_or_eq(a,b)
def is_greater_than_or_equal(a,b,timeout=1):
    @fork(timeout=timeout)
    def gr_or_eq(c,d):
        return bool(c>=d)
    return gr_or_eq(a,b)


def strictly_less_than(a,b,timeout=1,prec=53):
    is_it=is_strictly_less_than(a,b,timeout=timeout)
    if is_it=='NO DATA (timed out)':
        is_it=bool(a.n(prec)<b.n(prec))
    return is_it
def strictly_greater_than(a,b,timeout=1,prec=53):
    is_it=is_strictly_greater_than(a,b,timeout=timeout)
    if is_it=='NO DATA (timed out)':
        is_it=bool(a.n(prec)>b.n(prec))
    return is_it    
def less_than_or_equal(a,b,timeout=1,prec=53):
    is_it=is_less_than_or_equal(a,b,timeout=timeout)
    if is_it=='NO DATA (timed out)':
        is_it=bool(a.n(prec)<=b.n(prec))
    return is_it        
def greater_than_or_equal(a,b,timeout=1,prec=53):
    is_it=is_greater_than_or_equal(a,b,timeout=timeout)
    if is_it=='NO DATA (timed out)':
        is_it=bool(a.n(prec)>=b.n(prec))
    return is_it            

# To be used below so that holonomies may be sorted counter-clockwise
def ang(holonomy,prec=53):
    x=holonomy[0]
    y=holonomy[1]
    if x>0:
        if y>=0:
            theta=arctan(y/x)/(2*pi)
        else:
            theta=arctan(y/x)/(2*pi)+1
    elif x<0:
            theta=arctan(y/x)/(2*pi)+1/2
    else:
        if y>0:
            theta=1/4
        else: theta=3/4
    return theta.n(prec)

# Returns the maximum singular value of a matrix in SL_2(R) with Frobenius norm x
def nu(x):
    return sqrt((x^2+sqrt(x^4-4))/2)

# Determines the Voronoi staples of Surface by taking intersections of half-planes determined by a finite collection of marked segments (CandidateSegments)
def Staples(Surface,graphic=False):
    # Here we find the max length of an edge in the Delaunay triangulation of Surface
    Y=Surface.delaunay_triangulation()
    max_distances_in_polygons=[]
    for poly in range(Y.num_polygons()):
        distances_in_poly=[]
        vertex_pairs=[[i,j] for i in range(Y.polygon(poly).num_edges()) for j in range(i+1,Y.polygon(poly).num_edges())]
        for pair in vertex_pairs:
            distances_in_poly.append(sqrt((Y.polygon(poly).vertex(pair[1])[0]-Y.polygon(poly).vertex(pair[0])[0])^2+(Y.polygon(poly).vertex(pair[1])[1]-Y.polygon(poly).vertex(pair[0])[1])^2))
        max_distances_in_polygons.append(max(distances_in_poly))
    # Gives an upper bound on max distance between distinct vertices of polygons of Delaunay triangulation of Surface
    diameter_upper_bound=AA(max(max_distances_in_polygons))

    # Each element of CandidateSegments will be a list corresponding to a singularity P of Surface that contains info about all saddle connections leaving P with length bounded by diameter_upper_bound; this info is recoreded as [sc.holonomy(),sc,cumulative_angle] where sc.holonomy() is the holonomy vector in R^2 corresponding to the saddle connection, sc is the saddle connection as recored in flatsurf, and cumulative_angle is the cumulative angle of this saddle connection from the ray in direction (1,0) begining at the first vertex corresponding to singularity P in Surface.vertices_in_disjoint_planes
    CandidateSegments=[[] for P in range(Surface.num_singularities())]
    for P in range(Surface.num_singularities()):
        for copy_of_plane in range(Surface.cone_angles[P]):
            for vertex_data in Surface.vertices_in_disjoint_planes[P][copy_of_plane]:
                # Compute saddle connections of length bounded by ceiling of diameter_upper_bound beginning at vertex_data
                saddle_connections=Surface.saddle_connections(ceil(diameter_upper_bound^2),vertex_data[0],vertex_data[1])
                # First see if marked segments from the current vertex could develop in distinct copies of the plane; if so, we treat the cumulative angle accordingly
                if vertex_data==Surface.vertices_in_disjoint_planes[P][copy_of_plane][0]:
                    right_most_edge=Surface.polygon(vertex_data[0]).edge(vertex_data[1])
                    left_most_edge=-Surface.polygon(vertex_data[0]).edge((vertex_data[1]-1)%Surface.polygon(vertex_data[0]).num_edges())
                    right_most_edge_angle=ang(right_most_edge)
                    left_most_edge_angle=ang(left_most_edge)
                    for sc in saddle_connections:
                        holonomy_angle=ang(sc.holonomy())       
                        # If sc leaves its singularity clockwise of the vector (1,0), then the marked segment is in the previous copy of the plane and we assign the cumulative angle accordingly
                        if right_most_edge_angle>left_most_edge_angle and holonomy_angle>=right_most_edge_angle:
                            cumulative_angle=(copy_of_plane-1)%Surface.cone_angles[P]+holonomy_angle
                        # Otherwise its holonomy vector is in the current copy of the plane
                        else:
                            cumulative_angle=copy_of_plane+holonomy_angle
                        CandidateSegments[P].append([sc.holonomy(),sc,cumulative_angle])
                # These are the unambiguous vertices
                else:   
                    for sc in saddle_connections:
                        holonomy_angle=ang(sc.holonomy())
                        cumulative_angle=copy_of_plane+holonomy_angle
                        CandidateSegments[P].append([sc.holonomy(),sc,cumulative_angle])
        # Sort by cumulative angle
        CandidateSegments[P].sort(key=lambda x:x[2])

    # Will record perpendicular bisectors determining edges of Voronoi 2-cells; we record each of these as [[vertex_clockwise,vertex_counterclockwise],voronoi_staple] where vertex_clockwise (vertex_counterclockwise) is a pair [x,y] of x- and y-coordinates of the clockwise (counterclockwise) vertex that the edge determines, and voronoi_staple is the Voronoi staple (from CandidateSegments) determining this edge
    Voronoi2Cells=[[] for P in range(Surface.num_singularities())]
    for P in range(Surface.num_singularities()):
        # List of cumulative angles of marked segments corresponding to singularity P from CandidateSegments
        CandidateAngles=[candidate[2] for candidate in CandidateSegments[P]]            
        # We send out 'feeler' rays from the origin in various directions theta and include the closest element of CandidateSegments in direction theta as a side of the Voronoi 2 cell.  We continue until a compact region is found
        compact_cell=False
        feeler_directions=[k/5 for k in range(5)]
        Cell_P_indices=[]
        while compact_cell==False:
            for copy_of_plane in range(Surface.cone_angles[P]):
                for theta in feeler_directions:
                    tan_theta=tan(2*pi*theta)
                    # Find indices of all candidates in open pi/2-sector on either side of feeler ray in direction copy_of_plane+theta, as these are the only candidates whose perpendicular bisectors can intersect this ray
                    sector_indices=[]
                    for j in range(len(CandidateAngles)):
                        if abs(copy_of_plane+theta-CandidateAngles[j])<1/4 or abs(copy_of_plane+theta-CandidateAngles[j])>Surface.cone_angles[P]-1/4:
                            sector_indices.append(j)
                    # Of candidates in the sector, find which ones have perpendicular bisector intersecting the feeler ray nearest the origin
                    index_and_intersection_distance=[]
                    for j in sector_indices:
                        # Coordinates of holonomy of j^th candidate
                        x0=CandidateSegments[P][j][0][0]
                        y0=CandidateSegments[P][j][0][1]
                        # Algebra gives distance of intersection of perpendicular bisector of holonomy and feeler in direction theta (some segments may be pi/2 from theta; we weed those out here)
                        if (y0*tan_theta+x0).n()!=0:
                            intersection_distance=((x0^2+y0^2)/(2*(y0*tan_theta+x0)))^2*(1+tan_theta^2)
                            index_and_intersection_distance.append([j,intersection_distance])
                    # Sort by distance and find the closest intersections to the origin
                    index_and_intersection_distance.sort(key=lambda x:x[1])
                    closest=[]
                    for candidate in index_and_intersection_distance:
                        if candidate[1]==index_and_intersection_distance[0][1]:
                            closest.append(candidate)
                        else:
                            break
                    # If only one or two perpendicular bisectors are the closest and these have not already been included as edges of the 2-cell, then we include them
                    if len(closest)<3:
                        for candidate in closest:
                            if candidate[0] not in Cell_P_indices:
                                Cell_P_indices.append(candidate[0])
                    # Otherwise we must decide which two perpendicular bisectors 'cut into' the 2-cell the most and include these if they haven't already been included; these correspond to the clockwise- and clounterclockwise-most holonomies in question
                    else:
                        index_and_angle=[]
                        for candidate in closest:
                            index_and_angle.append([candidate[0],CandidateAngles[candidate[0]]])
                        index_and_angle.sort(key=lambda x:x[1])
                        # If the corresponding holonomies are not in both the first and final quadrant, then we take the first and last holonomies as sorted.  Otherwise we take the last in the first quadrant and the first in the final quadrant
                        if index_and_angle[-1][1]-index_and_angle[0][1]<1/2:
                            if index_and_angle[0][0] not in Cell_P_indices:
                                Cell_P_indices.append(index_and_angle[0][0])
                            if index_and_angle[-1][0] not in Cell_P_indices:
                                Cell_P_indices.append(index_and_angle[-1][0])
                        else:
                            j0=0
                            while index_and_angle[j0][1]<1/4:
                                j0+=1
                            if index_and_angle[j0][0] not in Cell_P_indices:
                                Cell_P_indices.append(index_and_angle[j0][0])
                            if index_and_angle[j0-1][0] not in Cell_P_indices:
                                Cell_P_indices.append(index_and_angle[j0-1][0])

            # We've found the perpendicular bisectors closest to the origin in the direction of each feeler; now we check if these determine a compact region.  In particular, we must have for each corresponding marked segment that the subsequent (counterclockwise) marked segment is within an angle of pi of the current one
            Cell_P_indices.sort()
            # A flag that will be made False if our region is determined to be non-compact
            is_compact=True
            if CandidateAngles[Cell_P_indices[-1]]-CandidateAngles[Cell_P_indices[0]]<=Surface.cone_angles[P]-1/2:
                is_compact=False
                # In the next iteration, we'll send a feeler ray between the first and last found Voronoi edges since these don't intersect
                direction=CandidateAngles[Cell_P_indices[-1]]+(CandidateAngles[Cell_P_indices[0]]+(Surface.cone_angles[P]-CandidateAngles[Cell_P_indices[-1]]))/2
                if direction<Surface.cone_angles[P]:
                    feeler_directions.append(direction-floor(direction))
                else:
                    direction=CandidateAngles[Cell_P_indices[0]]-(CandidateAngles[Cell_P_indices[0]]+(Surface.cone_angles[P]-CandidateAngles[Cell_P_indices[-1]]))/2
                    feeler_directions.append(direction-floor(direction))
            for j in range(len(Cell_P_indices)-1):
                if CandidateAngles[Cell_P_indices[j+1]]-CandidateAngles[Cell_P_indices[j]]>=1/2:
                    is_compact=False            
                    # In the next iteration, we'll send a feeler ray between the jth and (j+1)st found Voronoi edges since these don't intersect
                    direction=CandidateAngles[Cell_P_indices[j]]+(CandidateAngles[Cell_P_indices[j+1]]+-CandidateAngles[Cell_P_indices[j]])/2
                    feeler_directions.append(direction-floor(direction))
            if is_compact==True:
                compact_cell=True
            else:
                feeler_directions.sort()

        # Include the candidates found above as edges of Voronoi2Cells[P] and update edges based on intersections
        for j in Cell_P_indices:
            Voronoi2Cells[P].append([[[None,None],[None,None]],CandidateSegments[P][j]])
        for j in range(len(Voronoi2Cells[P])):
            # Find where subsequent edges intersect and update vertices
            staple0=Voronoi2Cells[P][j]
            x0=staple0[1][0][0]
            y0=staple0[1][0][1]
            staple1=Voronoi2Cells[P][(j+1)%len(Voronoi2Cells[P])]
            x1=staple1[1][0][0]
            y1=staple1[1][0][1]
            # Algebra gives where the corresponding edges determined by these staples intersect
            x_intersection=((x1^2+y1^2)*y0-(x0^2+y0^2)*y1)/(2*(x1*y0-x0*y1))
            if y0!=0:
                y_intersection=-(x0/y0)*(x_intersection-x0/2)+y0/2
            else:
                y_intersection=-(x1/y1)*(x_intersection-x1/2)+y1/2
            # Update counterclockwise vertex of staple0 and clockwise vertex of staple1
            staple0[0][1]=[x_intersection,y_intersection]
            staple1[0][0]=[x_intersection,y_intersection]

        # For each 2-cell edge, we search the remaining candidates in CandidateSegments to see which, if any, has a perpendicular bisector intersecting this edge.  We continue this search until no new edges are found
        previous_num_edges=len(Voronoi2Cells[P])-1
        while len(Voronoi2Cells[P])>previous_num_edges:
            previous_num_edges=len(Voronoi2Cells[P])
            previous_Cell_P_indices=deepcopy(Cell_P_indices)
            for j in range(len(previous_Cell_P_indices)):
                staple=Voronoi2Cells[P][j]
                x0=staple[1][0][0]
                y0=staple[1][0][1]
                # Indices of candidates between subsequent staples
                if j==len(previous_Cell_P_indices)-1:
                    sector_indices=list(range(previous_Cell_P_indices[0]))+list(range(previous_Cell_P_indices[j]+1,len(CandidateSegments[P])))
                else:
                    sector_indices=range(previous_Cell_P_indices[j]+1,previous_Cell_P_indices[j+1])
                indices_and_intersection_info=[]
                for k in sector_indices:
                    candidate=CandidateSegments[P][k]
                    x1=candidate[0][0]
                    y1=candidate[0][1]
                    # Algebra gives where the corresponding edges determined by these staples intersect
                    x_intersection=((x1^2+y1^2)*y0-(x0^2+y0^2)*y1)/(2*(x1*y0-x0*y1))
                    if y0!=0:
                        y_intersection=-(x0/y0)*(x_intersection-x0/2)+y0/2
                    else:
                        y_intersection=-(x1/y1)*(x_intersection-x1/2)+y1/2
                    # We want to consider only those perpendicular bisectors that intersect the staple's corresponding edge
                    if y0!=0:
                        if x_intersection>min(staple[0][0][0],staple[0][1][0]) and x_intersection<max(staple[0][0][0],staple[0][1][0]):
                            indices_and_intersection_info.append([k,x_intersection])
                    else:
                        if y_intersection>min(staple[0][0][1],staple[0][1][1]) and y_intersection<max(staple[0][0][1],staple[0][1][1]):
                            indices_and_intersection_info.append([k,y_intersection])
                if len(indices_and_intersection_info)>0:
                    # Find the clockwise-most intersection
                    indices_and_intersection_info.sort(key=lambda x:x[1])
                    clockwise_most=[]
                    if y0!=0:
                        # If y0>0, we want the right-most x_intersection
                        if y0>0:
                            indices_and_intersection_info.reverse()
                    else:
                        # If x0<0, we want the upper-most y_intersection
                        if x0<0:
                            indices_and_intersection_info.reverse()
                    for candidate in indices_and_intersection_info:
                        if candidate[1]==indices_and_intersection_info[0][1]:
                            clockwise_most.append(candidate)
                        else:
                            break
                    # If exactly one candidate's perpendicular bisector intersects the edge clockwise-most, then it is also an edge of the 2-cell
                    if len(clockwise_most)==1:
                        Cell_P_indices.append(clockwise_most[0][0])
                    # Otherwise, we need to find the one that 'cuts into' the cell the most; this is the candidate whose holonomy is counterclockwise-most 
                    else:
                        indices_and_holonomy_angles=[]
                        for candidate in clockwise_most:
                            indices_and_holonomy_angles.append([candidate[0],CandidateSegments[P][candidate[0]][2]])
                        indices_and_holonomy_angles.sort(key=lambda x:x[1])
                        # If holonomies are in both the first and final quadrants, we need to find the last one in the first quadrant.  Otherwise we choose the last one
                        if indices_and_holonomy_angles[0][1]<1/4 and indices_and_holonomy_angles[-1][1]>Surface.cone_angles[P]-1/4:
                            j0=0
                            angle_j0=indices_and_holonomy_angles[j0][1]
                            while angle_j0<1/4:
                                j0+=1
                                angle_j0=indices_and_holonomy_angles[j0][1]
                            Cell_P_indices.append(indices_and_holonomy_angles[j0-1][0])
                        else:   
                            Cell_P_indices.append(indices_and_holonomy_angles[-1][0])

            # Update edges and vertices
            Cell_P_indices.sort()
            Voronoi2Cells[P]=[]
            for j in Cell_P_indices:
                Voronoi2Cells[P].append([[[None,None],[None,None]],CandidateSegments[P][j]])
            for j in range(len(Voronoi2Cells[P])):
                # Find where subsequent edges intersect and update vertices
                staple0=Voronoi2Cells[P][j]
                x0=staple0[1][0][0]
                y0=staple0[1][0][1]
                staple1=Voronoi2Cells[P][(j+1)%len(Voronoi2Cells[P])]
                x1=staple1[1][0][0]
                y1=staple1[1][0][1]
                # Algebra gives where the corresponding edges determined by these staples intersect
                x_intersection=((x1^2+y1^2)*y0-(x0^2+y0^2)*y1)/(2*(x1*y0-x0*y1))
                if y0!=0:
                    y_intersection=-(x0/y0)*(x_intersection-x0/2)+y0/2
                else:
                    y_intersection=-(x1/y1)*(x_intersection-x1/2)+y1/2
                # Update counterclockwise vertex of staple0 and clockwise vertex of staple1
                staple0[0][1]=[x_intersection,y_intersection]
                staple1[0][0]=[x_intersection,y_intersection]

    # We've found all the Voronoi staples as segments determining edges of Voronoi2Cells; now we represent these in a list separated by singularity and copy of plane
    voronoi_staples=[[[] for copy_of_plane in range(Surface.cone_angles[P])] for P in range(Surface.num_singularities())]
    for P in range(Surface.num_singularities()):
        j=0
        angle=Voronoi2Cells[P][j][1][2]
        for copy_of_plane in range(Surface.cone_angles[P]):
            while angle<copy_of_plane+1:
                voronoi_staples[P][copy_of_plane].append(Voronoi2Cells[P][j][1])
                j+=1
                if j<len(Voronoi2Cells[P]):
                    angle=Voronoi2Cells[P][j][1][2]
                else:
                    break

    if graphic==False:
        return voronoi_staples
    else:
        VoronoiGraphics=[[line([(0,0),(0,0)],aspect_ratio=1) for copy_of_plane in range(Surface.cone_angles[P])] for P in range(Surface.num_singularities())]
        polygons=[[[[],[]] for copy_of_plane in range(Surface.cone_angles[P])] for P in range(Surface.num_singularities())]
        for P in range(Surface.num_singularities()):
            j=0
            angle=Voronoi2Cells[P][j][1][2]
            for copy_of_plane in range(Surface.cone_angles[P]):
                while angle<copy_of_plane+1:
                    straight_line=Voronoi2Cells[P][j][0]
                    x0=straight_line[0][0]
                    y0=straight_line[0][1]
                    x1=straight_line[1][0]
                    y1=straight_line[1][1]              
                    # If straight_line is vertical
                    if x0==x1:
                        # If the straight_line crosses the positive x-axis, it is potentially represented in multiple copies of the plane
                        if x0>0:
                            if y0<0 and y1>0:
                                VoronoiGraphics[P][copy_of_plane]+=line([(x0,0),(x1,y1)],aspect_ratio=1)
                                VoronoiGraphics[P][(copy_of_plane-1)%Surface.cone_angles[P]]+=line([(x0,y0),(x1,0)],aspect_ratio=1)                             
                            elif y0>=0:
                                VoronoiGraphics[P][copy_of_plane]+=line(straight_line,aspect_ratio=1)
                            else:
                                VoronoiGraphics[P][(copy_of_plane-1)%Surface.cone_angles[P]]+=line(straight_line,aspect_ratio=1)
                        # Otherwise it's in one copy of the plane   
                        else:
                            VoronoiGraphics[P][copy_of_plane]+=line(straight_line,aspect_ratio=1)
                    # If straight_line is horizontal
                    elif y0==y1:
                        VoronoiGraphics[P][copy_of_plane]+=line(straight_line,aspect_ratio=1)
                    # Otherwise
                    else:
                        x_intercept=x0-((x1-x0)/(y1-y0))*y0
                        # If extension of straight_line crosses the positive x-axis, it is potentially represented in multiple copies of the plane
                        if x_intercept>0:
                            if y0<0 and y1>0:
                                if angle-copy_of_plane<1/2:
                                    VoronoiGraphics[P][copy_of_plane]+=line([(x_intercept,0),(x1,y1)],aspect_ratio=1)
                                    VoronoiGraphics[P][(copy_of_plane-1)%Surface.cone_angles[P]]+=line([(x0,y0),(x_intercept,0)],aspect_ratio=1)
                                else:
                                    VoronoiGraphics[P][(copy_of_plane+1)%Surface.cone_angles[P]]+=line([(x_intercept,0),(x1,y1)],aspect_ratio=1)
                                    VoronoiGraphics[P][copy_of_plane]+=line([(x0,y0),(x_intercept,0)],aspect_ratio=1)                                   
                            elif y0>=0:
                                if angle-copy_of_plane<1/2:
                                    VoronoiGraphics[P][copy_of_plane]+=line(straight_line,aspect_ratio=1)
                                else:
                                    VoronoiGraphics[P][(copy_of_plane+1)%Surface.cone_angles[P]]+=line(straight_line,aspect_ratio=1)
                            else:
                                if angle-copy_of_plane<1/2:
                                    VoronoiGraphics[P][(copy_of_plane-1)%Surface.cone_angles[P]]+=line(straight_line,aspect_ratio=1)
                                else:
                                    VoronoiGraphics[P][copy_of_plane]+=line(straight_line,aspect_ratio=1)
                        # Otherwise it's in one copy of the plane
                        else:
                            VoronoiGraphics[P][copy_of_plane]+=line(straight_line,aspect_ratio=1)   
                    if j+1<len(Voronoi2Cells[P]):
                        j+=1
                        angle=Voronoi2Cells[P][j][1][2]
                    else:
                        break
        return [voronoi_staples,VoronoiGraphics]    

# Returns all marked segments corresponding to singularity P of radius less than or equal to Radii[P]
def marked_segments(Surface,Radii):
    MSr=[]
    for P in range(Surface.num_singularities()):
        MSr_P=[[] for i in range(Surface.cone_angles[P])]
        for copy_of_plane in range(len(Surface.vertices_in_disjoint_planes[P])):
            for vertex_data in Surface.vertices_in_disjoint_planes[P][copy_of_plane]:           
                polygon_start=vertex_data[0]
                vertex_start=vertex_data[1]
                # saddle_connections=load('sc{}{}.sobj'.format(polygon_start,vertex_start))
                saddle_connections=Surface.saddle_connections(ceil(Radii[P]^2),polygon_start,vertex_start)
                # First see if the current vertex is ambiguous (in the sense that saddle connections from it may be assigned to the current copy of the plane or the previous copy)
                if vertex_data==Surface.vertices_in_disjoint_planes[P][copy_of_plane][0]:
                    right_most_edge=Surface.polygon(polygon_start).edge(vertex_start)
                    left_most_edge=-Surface.polygon(polygon_start).edge((vertex_start-1)%Surface.polygon(polygon_start).num_edges())
                    right_most_edge_angle=ang(right_most_edge)
                    left_most_edge_angle=ang(left_most_edge)
                    for i in range(len(saddle_connections)):
                        sc=saddle_connections[i]
                        holonomy=sc.holonomy()
                        holonomy_angle=ang(holonomy)
                        # If sc leaves its singularity clockwise of the vector (1,0), then its holonomy vector is assigned to the previous copy of the plane
                        if right_most_edge_angle>left_most_edge_angle and holonomy_angle>=right_most_edge_angle:
                            cumulative_angle=(copy_of_plane-1)%Surface.cone_angles[P]+holonomy_angle
                            MSr_P[copy_of_plane-1].append([holonomy,sc,cumulative_angle])
                        # Otherwise its holonomy vector is assigned to the current copy of the plane
                        else:
                            cumulative_angle=copy_of_plane+holonomy_angle
                            MSr_P[copy_of_plane].append([holonomy,sc,cumulative_angle])
                # These are the unambiguous vertices
                else:   
                    for i in range(len(saddle_connections)):
                        sc=saddle_connections[i]
                        holonomy=sc.holonomy()
                        holonomy_angle=ang(holonomy)
                        cumulative_angle=copy_of_plane+holonomy_angle
                        MSr_P[copy_of_plane].append([holonomy,sc,cumulative_angle])
                # print(vertex_data)
        MSr.append(MSr_P)
        # Sort holonomies in each copy of the plane by cumulative counterclockwise angle from (1,0) on 0th copy of plane
        for copy_of_plane in range(len(MSr[P])):
            MSr[P][copy_of_plane].sort(key=lambda x:x[2])
    return MSr          


# Returns all elements of the Veech group whose Frobenius norm is bounded above by NormBound
# PreviousSLa is an optional list of matrices already found to be in the Veech group with a smaller NormBound
# MatricesChecked is a list of matrices already tested for a smaller NormBound
# MatricesToCheck is a list of matrices to check that were found in a previous iteration but did not satisfy the norm bound at the time
def SLa(Surface,NormBound,PreviousSLa=[matrix.identity(2)],MatricesChecked=[matrix.identity(2)],MatricesToCheck=[],PreviousNormBound=None):
    # MatricesToCheck_Now will be a list comprised of matrices in $\text{SL}_2\mathbb{R}$ whose Frobenius norms are bounded above by NormBound and whose inverses take a basis of Voronoi staples to another basis of vectors in marked_periods(radius) for the appropriate radius
    MatricesToCheck_Now=[]
    MatricesToCheck_Later=MatricesToCheck
    remove=[]
    for j in range(len(MatricesToCheck_Later)):
        M=MatricesToCheck_Later[j]
        if M[0][0]^2+M[0][1]^2+M[1][0]^2+M[1][1]^2<=NormBound^2:
            MatricesToCheck_Now.append(M)
            remove.append(j)
    remove.reverse()
    for j in remove:
        MatricesToCheck_Later.pop(j)
    RadiiUpper=[r*nu(NormBound) for r in Surface.radii]

    # Marked segments of length bounded by RadiiUpper[P] for each singularity P
    MS_Radii=marked_segments(Surface,RadiiUpper)
    # List of all holonomies of marked segments bounded by Radii
    MS_Radii_holonomies=[MS_Radii[P][copy_of_plane][j][0] for P in range(len(MS_Radii)) for copy_of_plane in range(len(MS_Radii[P])) for j in range(len(MS_Radii[P][copy_of_plane]))]   

    # We find which cone angle produces the least number of marked segments and find candidate matrices based on singularities of this cone angle
    num_marked_segments_of_same_angle=[]
    num_marked_segments_of_current_angle=0
    previous_cone_angle=Surface.cone_angles[0]
    for P in range(Surface.num_singularities()):
        if Surface.cone_angles[P]==previous_cone_angle:
            num_marked_segments_of_current_angle+=sum([len(MS_Radii[P][copy_of_plane]) for copy_of_plane in range(len(MS_Radii[P]))])
        else:
            num_marked_segments_of_same_angle.append(num_marked_segments_of_current_angle)
            num_marked_segments_of_current_angle=sum([len(MS_Radii[P][copy_of_plane]) for copy_of_plane in range(len(MS_Radii[P]))])
        if P==Surface.num_singularities()-1:
            num_marked_segments_of_same_angle.append(num_marked_segments_of_current_angle)
        previous_cone_angle=Surface.cone_angles[P]
    index=num_marked_segments_of_same_angle.index(min(num_marked_segments_of_same_angle))
    j=0
    previous_cone_angle=Surface.cone_angles[0]
    for P in range(Surface.num_singularities()):
        if Surface.cone_angles[P]!=previous_cone_angle:
            j+=1
        if j==index:
            # These are the singularities of cone angle producing the minimal number of marked segments
            singularity_range=range(P,P+Surface.cone_angles.count(Surface.cone_angles[P]))
            break
        previous_cone_angle=Surface.cone_angles[P]

    # Basis from projection of Surface.VoronoiStaples onto complex plane; note that since the Voronoi staples are sorted by cumulative angle, v0,v1 have standard orientation
    v0=Surface.VoronoiStaples[singularity_range[0]][0][0][0]
    v1=Surface.VoronoiStaples[singularity_range[0]][0][1][0]
    # Lengths (squared) of v0,v1
    max_basis_len=max([v0[0]^2+v0[1]^2,v1[0]^2+v1[1]^2])
    min_basis_len=min([v0[0]^2+v0[1]^2,v1[0]^2+v1[1]^2])
    # lower_radius is the most the shortest of v0,v1 could have expanded under a matrix with PreviousNormBound
    if PreviousNormBound==None:
        lower_radius=0
    else:
        lower_radius=min_basis_len*(nu(PreviousNormBound))^2
    T=matrix([[v0[0],v1[0]],[v0[1],v1[1]]])
    sing=0
    for P in singularity_range:
        print_percent_complete('Singularity being checked in search for candiate matrices of norm bound {}:'.format(NormBound),sing,len(singularity_range))
        sing+=1
        # Single lists of marked periods in MS_Radii[P], sorted by cumulative angle
        projections_to_plane=[]
        for copy_of_plane in range(Surface.cone_angles[P]):
            for ms in MS_Radii[P][copy_of_plane]:
                # NormBound restricts how much v0,v1 can stretch, so we only search among marked segments satisfying the following constraint
                if ms[0][0]^2+ms[0][1]^2<=max_basis_len*(nu(NormBound))^2:
                    projections_to_plane.append(ms)
        # Construct new basis w0,w1 where the length of w0 is greater than lower_radius (to avoid repitition from previous runs, since if the lengths of both w0,w1 are less than lower_radius, the resulting matrix would have been detected in a previous iteration)
        for j in range(len(projections_to_plane)):
            w0=projections_to_plane[j]
            w0_hol=w0[0]
            if w0_hol[0]^2+w0_hol[1]^2>lower_radius:
                # First we search for w1 within pi-sector counterclockwise of w0
                k=j+1
                w1=projections_to_plane[k%len(projections_to_plane)]
                w1_hol=w1[0]
                S=matrix([[w0_hol[0],w1_hol[0]],[w0_hol[1],w1_hol[1]]])
                # Since projections_to_plane are sorted by cumulative angle, we know that once det(S) is non-postive, we're no longer in the necessary pi-sector
                while det(S)>0: 
                    # M sends v0,v1 to w0,w1, respectively 
                    M=S*T.inverse()
                    if det(M)==1 and M not in MatricesChecked and M not in MatricesToCheck_Now and M not in MatricesToCheck_Later:
                        if M[0][0]^2+M[0][1]^2+M[1][0]^2+M[1][1]^2<=NormBound^2:
                            MatricesToCheck_Now.append(M)
                        else:
                            MatricesToCheck_Later.append(M)
                    k+=1
                    w1=projections_to_plane[k%len(projections_to_plane)]
                    w1_hol=w1[0]
                    S=matrix([[w0_hol[0],w1_hol[0]],[w0_hol[1],w1_hol[1]]])             
                # Next we search for w1 within pi-sector clockwise of w0
                k=j-1
                w1=projections_to_plane[k%len(projections_to_plane)]
                w1_hol=w1[0]
                S=matrix([[w1_hol[0],w0_hol[0]],[w1_hol[1],w0_hol[1]]])
                # Since projections_to_plane are sorted by cumulative angle, we know that once det(S) is non-postive, we're no longer in the necessary pi-sector
                while det(S)>0: 
                    # M sends v0,v1 to w1,w0, respectively 
                    M=S*T.inverse()                     
                    if det(M)==1 and M not in MatricesChecked and M not in MatricesToCheck_Now and M not in MatricesToCheck_Later:
                        if M[0][0]^2+M[0][1]^2+M[1][0]^2+M[1][1]^2<=NormBound^2:
                            MatricesToCheck_Now.append(M)
                        else:
                            MatricesToCheck_Later.append(M)
                    k-=1
                    w1=projections_to_plane[k%len(projections_to_plane)]
                    w1_hol=w1[0]
                    S=matrix([[w1_hol[0],w0_hol[0]],[w1_hol[1],w0_hol[1]]])
            print_percent_complete('    Holonomies checked emenating from this singularity:',j,len(projections_to_plane))                              

    # SL will keep the matrices with Frobenius norm bounded by NormBound that are in the Veech group
    SL=PreviousSLa                      

    # We find products of matrices already found to be in the Veech group whose Frobenius norms are bounded by NormBound
    len_PreviousSLa=len(PreviousSLa)
    for j in range(len_PreviousSLa):
        M=SL[j]
        for k in range(j,len_PreviousSLa):
            N=SL[k]
            if M*N in MatricesToCheck_Now:
                MatricesToCheck_Now.remove(M*N) 
                MatricesChecked.append(M*N)
                SL.append(M*N)
            if N*M in MatricesToCheck_Now:
                MatricesToCheck_Now.remove(N*M) 
                MatricesChecked.append(N*M)
                SL.append(N*M)

    # Since MatricesToCheck_Now contains all (unchecked) matrices in the Veech group with Frobenius norm bounded by NormBound, if M and N belong to Veech group, then if M*N (or N*M) has Frobenius norm bounded by NormBound, this product must belong to Veech group and must be represented in MatricesToCheck_Now or SL; we get rid of matrices M that do not satisfy this criterion
    remove=[]
    for j in range(len(MatricesToCheck_Now)):
        M=MatricesToCheck_Now[j]
        for N in SL:
            K=M*N
            L=N*M
            if K[0][0]^2+K[0][1]^2+K[1][0]^2+K[1][1]^2<=NormBound^2 and K not in MatricesToCheck_Now and K not in SL:
                remove.append(j)
                MatricesChecked.append(M)
                break
            elif L[0][0]^2+L[0][1]^2+L[1][0]^2+L[1][1]^2<=NormBound^2 and L not in MatricesToCheck_Now and L not in SL:
                remove.append(j)
                MatricesChecked.append(M)
                break
        # print_percent_complete('matrix',j,len(MatricesToCheck_Now))
    remove.reverse()
    for j in remove:
        MatricesToCheck_Now.pop(j)                                              
                                
    # We test Veech group membership of each M in MatricesToCheck_Now
    MatricesToCheck_Now_copy=deepcopy(MatricesToCheck_Now)
    num_matrices=len(MatricesToCheck_Now_copy)
    current_iteration=0
    for M in MatricesToCheck_Now_copy[:6000]:
        MatricesToCheck_Now.remove(M)
        if M not in MatricesChecked:
            M_in_SL_flag=False
            holonomies_agree=True
            MatricesChecked.append(M)

            # Apply M to VoronoiStaples
            M_VoronoiStaples=[]
            angle_of_image_of_e1=ang(M*vector((1,0)))
            for P in range(Surface.num_singularities()):
                M_VoronoiStaples_P=[[] for c in range(Surface.cone_angles[P])]
                for copy_of_plane in range(Surface.cone_angles[P]):
                    for j in range(len(Surface.VoronoiStaples[P][copy_of_plane])):
                        # Holonomy of a staple and its image holonomy under M
                        hol=Surface.VoronoiStaples[P][copy_of_plane][j][0]
                        M_hol=M*vector([hol[0],hol[1]])
                        # M_hol must be in MS_Radii_holonomies for M to possibly be in the Veech group; we use this as a quick screen while applying M to each staple
                        if M_hol in MS_Radii_holonomies:
                            angle_of_M_hol=ang(M_hol)
                            # If the angle of M_hol is greater than or equal to angle_of_image_of_e1, then we assign this marked segment to this copy_of_plane; otherwise,it gets assigned to the next copy_of_plane
                            # We store this data as a list: [M_hol, hol's corresponding saddle connection on Surface (as recorded in flatsurf), cumulative angle of M_hol from (1,0) in the 0th copy of the plane, cumulative angle of hol from (1,0) in 0th copy of plane]
                            if angle_of_M_hol>=angle_of_image_of_e1:
                                M_VoronoiStaples_P[copy_of_plane].append([M_hol,Surface.VoronoiStaples[P][copy_of_plane][j][1],copy_of_plane+angle_of_M_hol,Surface.VoronoiStaples[P][copy_of_plane][j][2]])
                            else:
                                M_VoronoiStaples_P[(copy_of_plane+1)%Surface.cone_angles[P]].append([M_hol,Surface.VoronoiStaples[P][copy_of_plane][j][1],(copy_of_plane+1)%Surface.cone_angles[P]+angle_of_M_hol,Surface.VoronoiStaples[P][copy_of_plane][j][2]])
                        else:
                            holonomies_agree=False
                            break   
                    if holonomies_agree==False:
                        break
                if holonomies_agree==False:
                    break
                # Sort holonomies in each copy of the plane by counterclockwise angle
                for copy_of_plane in range(Surface.cone_angles[P]):
                    M_VoronoiStaples_P[copy_of_plane].sort(key=lambda x:x[2])
                M_VoronoiStaples.append(M_VoronoiStaples_P)

            # We only continue testing M if the image of Voronoi staple holonomies under M is contained in holonomies of marked segments bounded by Radii
            if holonomies_agree==True:
                # Now we test various permutations of singularities of the same cone angle within M_VoronoiStaples as well as cyclic permutations of copies of the plane corresponding to each singularity
                # We first check if the holonomies of a permutation of M_VoronoiStaples match those of MS_Radii; if so, we then test to see if the Z2-action is the same
                for permute_singularities in Surface.permutations_of_singularities.list():
                    M_VoronoiStaples_copy=deepcopy(M_VoronoiStaples)
                    place=0
                    for n in range(len(Surface.permutation_group_orders)):
                        num_sings_of_same_angle=Surface.permutation_group_orders[n]
                        M_VoronoiStaples_copy[place:place+num_sings_of_same_angle]=permute_singularities[n](M_VoronoiStaples_copy[place:place+num_sings_of_same_angle])
                        place+=num_sings_of_same_angle
                    # For each permutation of singularities, we cyclically permute the copies of the plane associated to a particular singularity
                    for permute_planes in Surface.permutations_of_copies_of_plane.list():
                        end_permutation_check=False
                        M_VoronoiStaples_TransO=deepcopy(M_VoronoiStaples_copy)
                        for P in range(Surface.num_singularities()):
                            M_VoronoiStaples_TransO[P]=permute_planes[P](M_VoronoiStaples_TransO[P])                        
                        # Find only the marked segments coinciding with M_VoronoiStaples_TransO (the image of M_VoronoiStaples under the current permutations)
                        MS_Radii_agree=[[[] for copy_of_plane in range(Surface.cone_angles[P])] for P in range(Surface.num_singularities())]
                        for P in range(Surface.num_singularities()):
                            for copy_of_plane in range(Surface.cone_angles[P]):
                                for marked_segment in MS_Radii[P][copy_of_plane]:
                                    if marked_segment[0] in [vs[0] for vs in M_VoronoiStaples_TransO[P][copy_of_plane]]:
                                        MS_Radii_agree[P][copy_of_plane].append(marked_segment)
                                if len(MS_Radii_agree[P][copy_of_plane])!=len(M_VoronoiStaples_TransO[P][copy_of_plane]):
                                    end_permutation_check=True
                                    break
                            if end_permutation_check==True:
                                break
                        # If we've made it this far, the holonomies agree; we need to check the Z2-action
                        if end_permutation_check==False:
                            # List of saddle connection data; this is not associated to a particular singularity and copy of the plane
                            M_VoronoiStaples_TransO_saddle_connections=[M_VoronoiStaples_TransO[P][copy_of_plane][j][1] for P in range(Surface.num_singularities()) for copy_of_plane in range(Surface.cone_angles[P]) for j in range(len(M_VoronoiStaples_TransO[P][copy_of_plane]))]
                            MS_Radii_agree_saddle_connections=[MS_Radii_agree[P][copy_of_plane][j][1] for P in range(Surface.num_singularities()) for copy_of_plane in range(Surface.cone_angles[P]) for j in range(len(MS_Radii_agree[P][copy_of_plane]))]
                            Z2_action_agrees_flag=True
                            for j in range(len(M_VoronoiStaples_TransO_saddle_connections)):
                                # The saddle connections on Surface corresponding to the same marked semgment
                                M_VoronoiStaples_TransO_sc=M_VoronoiStaples_TransO_saddle_connections[j]
                                MS_Radii_agree_sc=MS_Radii_agree_saddle_connections[j]
                                # These are the identical but oppositely oriented saddle connections of those found above
                                M_VoronoiStaples_TransO_sc_inverse=M_VoronoiStaples_TransO_sc.invert()
                                MS_Radii_agree_sc_inverse=MS_Radii_agree_sc.invert()
                                # For the Z2-action to match, these inverse saddle connections must have the same indices; otherwise the Z2-action is different and M_VoronoiStaples_TransO does not equal MS_Radii_agree
                                if M_VoronoiStaples_TransO_saddle_connections.index(M_VoronoiStaples_TransO_sc_inverse)!=MS_Radii_agree_saddle_connections.index(MS_Radii_agree_sc_inverse):
                                    Z2_action_agrees_flag=False
                                    break
                            if Z2_action_agrees_flag==True:
                                M_in_SL_flag=True
                                SL.append(M)
                                # print(M)
                                # For N in SL, if M*N is in MatricesToCheck_Now, then M is in SL if and only if M*N is in SL (and sim. for N*M and replacing M with M.inverse())
                                previous_len_SL=len(SL)-1
                                first_run=True
                                while len(SL)>previous_len_SL:
                                    if first_run==True:
                                        index=0
                                        first_run=False
                                    else:
                                        index=previous_len_SL
                                    previous_len_SL=len(SL)
                                    SL_copy=deepcopy(SL)    
                                    for N in SL_copy[index:]:
                                # j=0
                                # while j<len(SL):
                                #     j+=1
                                #     for N in SL[j:]
                                # previous_len_SL=len(SL)-1
                                # while len(SL)>previous_len_SL:
                                    # previous_len_SL=len(SL)
                                    # SL_copy=deepcopy(SL)
                                    # for N in SL_copy:
                                        if M*N in MatricesToCheck_Now and M*N not in MatricesChecked:
                                            MatricesChecked.append(M*N)
                                            SL.append(M*N)
                                        if N*M in MatricesToCheck_Now and N*M not in MatricesChecked:
                                            MatricesChecked.append(N*M)
                                            SL.append(N*M)                                                  
                                        if M.inverse()*N in MatricesToCheck_Now and M.inverse()*N not in MatricesChecked:
                                            MatricesChecked.append(M.inverse()*N)
                                            SL.append(M.inverse()*N)
                                        if N*M.inverse() in MatricesToCheck_Now and N*M.inverse() not in MatricesChecked:
                                            MatricesChecked.append(N*M.inverse())
                                            SL.append(N*M.inverse())
                                break
                    if M_in_SL_flag==True:
                        break
            # If M not in SL and N in SL, then M*N in MatricesToCheck_Now not in SL, for otherwise M=(M*N)*N.inverse() is in SL
            if M_in_SL_flag==False:
                for N in SL:
                    if M*N in MatricesToCheck_Now and M*N not in MatricesChecked:
                        MatricesChecked.append(M*N)
                    if N*M in MatricesToCheck_Now and N*M not in MatricesChecked:
                        MatricesChecked.append(N*M)
                    if M.inverse()*N in MatricesToCheck_Now and M.inverse()*N not in MatricesChecked:
                        MatricesChecked.append(M.inverse()*N)
                    if N*M.inverse() in MatricesToCheck_Now and N*M.inverse() not in MatricesChecked:
                        MatricesChecked.append(N*M.inverse())
        print_percent_complete('Matrices checked:'.format(NormBound),current_iteration,num_matrices)
        current_iteration+=1                        
    #print([SL,RadiiUpper,NormBound])
    return([SL,MatricesChecked,MatricesToCheck_Later])          

# Takes as input list of (closed) intervals called 'real_subset' and another (closed) interval called 'interval', and outputs the (closure) of the difference as a list of intervals
def difference(real_subset,interval,timeout=1,prec=53):
    if real_subset!=[]:
        c=interval[0]
        d=interval[1]
        n=len(real_subset)-1
        # List of left endpoints of intervals in real_subset
        a=[J[0] for J in real_subset]
        # List of right endpoints of intervals in real_subset
        b=[J[1] for J in real_subset]

        if less_than_or_equal(c,a[0],timeout=timeout,prec=prec):
            left=[]
        elif strictly_greater_than(c,a[0],timeout=timeout,prec=prec) and less_than_or_equal(c,b[0],timeout=timeout,prec=prec):
            left=[[a[0],c]]
        elif greater_than_or_equal(c,b[n],timeout=timeout,prec=prec):
            left=real_subset
        else:
            for j in range(1,n+1):
                if greater_than_or_equal(c,b[j-1],timeout=timeout,prec=prec) and less_than_or_equal(c,a[j],timeout=timeout,prec=prec):
                    left=real_subset[:j]
                    break
                elif strictly_greater_than(c,a[j],timeout=timeout,prec=prec) and less_than_or_equal(c,b[j],timeout=timeout,prec=prec):
                    left=real_subset[:j]+[[a[j],c]]
                    break

        if greater_than_or_equal(d,b[n],timeout=timeout,prec=prec):
            right=[]
        elif greater_than_or_equal(d,a[n],timeout=timeout,prec=prec) and strictly_less_than(d,b[n],timeout=timeout,prec=prec):
            right=[[d,b[n]]]
        elif less_than_or_equal(d,a[0],timeout=timeout,prec=prec):
            right=real_subset
        else:
            for j in range(n):
                if greater_than_or_equal(d,b[j],timeout=timeout,prec=prec) and less_than_or_equal(d,a[j+1],timeout=timeout,prec=prec):
                    right=real_subset[j+1:]
                    break
                elif greater_than_or_equal(d,a[j],timeout=timeout,prec=prec) and strictly_less_than(d,b[j],timeout=timeout,prec=prec):
                    right=[[d,b[j]]]+real_subset[j+1:]
                    break
        return left+right

    else:
        return real_subset

# This function constructs the hyperbolic polygon formed as the intersection of half-planes containing I that are defined by the inputed geodesics
def Omega(geodesics,prec=53,timeout=1):
    geodesics_copy=deepcopy(geodesics)

    # This will be needed below to find all geodesics exposing I that form sides of polygon
    FreeSides=[[-infinity,infinity]]

    # We sort the given geodesics by type
    geodesics_vertical=[]
    geodesics_enclosing_I=[]
    geodesics_exposing_I=[]

    for gamma in geodesics_copy:
        # We will flip this to be True when various gammas in geodesics_copy are shown to correspond to sides of the polygon
        gamma.is_polygon_side=False
        # We by default set the vertices of each gamma to be where it intersects the boundary of the upper-half plane; these will be updated as geodesics are shown to intersect as sides of the polygon
        if gamma.type=='Vertical':
            if gamma.foot<0:
                gamma.vertex_counterclockwise=[gamma.foot,0]
                gamma.vertex_clockwise=[gamma.foot,infinity]
            else:
                gamma.vertex_counterclockwise=[gamma.foot,infinity]
                gamma.vertex_clockwise=[gamma.foot,0]
            geodesics_vertical.append(gamma)
        elif gamma.type=='Encloses I':
            gamma.vertex_counterclockwise=[gamma.foot_left,0]
            gamma.vertex_clockwise=[gamma.foot_right,0]
            geodesics_enclosing_I.append(gamma)
        elif gamma.type=='Exposes I':
            gamma.vertex_counterclockwise=[gamma.foot_right,0]
            gamma.vertex_clockwise=[gamma.foot_left,0]
            geodesics_exposing_I.append(gamma)

    # Sort the geodesics by distance from I             
    geodesics_vertical.sort(key=lambda x:x.distance_from_I.n(prec))
    geodesics_enclosing_I.sort(key=lambda x:x.distance_from_I.n(prec))
    geodesics_exposing_I.sort(key=lambda x:x.distance_from_I.n(prec))

    # Will record the geodesics contributing to sides of the polygon
    polygon=[]

    # Find the geodesics enclosing I (if any) that are closest to I, and include them as sides of polygon
    num_close_geods_enclosing_I=0
    if len(geodesics_enclosing_I)>0:
        enclosing_min_distance=geodesics_enclosing_I[0].distance_from_I.n(prec)
        j=0
        gamma=geodesics_enclosing_I[0]
        while gamma.distance_from_I.n(prec)==enclosing_min_distance:
            num_close_geods_enclosing_I+=1
            for side in polygon:
                # Update the vertices of gamma and side if they intersect between the already existing vertices of side
                gamma.intersect_polygon_side(side,timeout=timeout,prec=prec)
            FreeSides=difference(FreeSides,[-infinity,gamma.foot_left],timeout=timeout,prec=prec)
            FreeSides=difference(FreeSides,[gamma.foot_right,infinity],timeout=timeout,prec=prec)
            # Since gamma is one of the closest geodesics to I, it is a polygon side regardless of whether or not it intersects already found polygon sides
            gamma.is_polygon_side=True
            polygon.append(gamma)           
            j+=1
            if j==len(geodesics_enclosing_I):
                break
            else:
                gamma=geodesics_enclosing_I[j]
        geodesics_enclosing_I=geodesics_enclosing_I[num_close_geods_enclosing_I:]

    # Find the geodesics exposing I that are closest to I, and include them as sides of polygon
    exposing_min_distance=geodesics_exposing_I[0].distance_from_I.n(prec)
    num_close_geods_exposing_I=0
    j=0
    gamma=geodesics_exposing_I[0]
    while gamma.distance_from_I.n(prec)==exposing_min_distance:
        num_close_geods_exposing_I+=1
        for side in polygon:
            # Update the vertices of gamma and side if they intersect between the already existing vertices of side
            gamma.intersect_polygon_side(side,timeout=timeout,prec=prec)
        FreeSides=difference(FreeSides,[gamma.foot_left,gamma.foot_right],timeout=timeout,prec=prec)
        # Since gamma is one of the closest geodesics to I, it is a polygon side regardless of whether or not it intersects already found polygon sides
        gamma.is_polygon_side=True
        polygon.append(gamma)
        j+=1
        if j==len(geodesics_exposing_I):
            break
        else:
            gamma=geodesics_exposing_I[j]
    geodesics_exposing_I=geodesics_exposing_I[num_close_geods_exposing_I:]

    # Find the vertical sides contributing to polygon, if any
    if len(geodesics_vertical)>=2:
        # We only need the closest two vertical geodesics
        for gamma in geodesics_vertical[:2]:
            # If there are no geodesics enclosing I, then we must include the two closest vertical geodesics as sides of the polygon for the polygon to have finite area
            if num_close_geods_enclosing_I==0:
                gamma.is_polygon_side=True
            # Otherwise, the vertical geodesics must be shown to intersect a segment of some polygon side contributing to the polygon (even if there are no geodesics enclosing I, we continue with this process to see if gamma intersects any sides and if we need to update vertices)
            for side in polygon:
                # Check if gamma and side intersect, and update the vertices of both geodesics if so
                gamma.intersect_polygon_side(side,timeout=timeout,prec=prec)

            if gamma.is_polygon_side==True:
                if gamma.foot<0:
                    FreeSides=difference(FreeSides,[-infinity,gamma.foot],timeout=timeout,prec=prec)
                else:
                    FreeSides=difference(FreeSides,[gamma.foot,infinity],timeout=timeout,prec=prec)
                polygon.append(gamma)

    # Now we see which other geodesics might contribute to the sides of polygon; these sides are not the closest sides of polygon to I, as those were all found above
    first_run=True
    while len(FreeSides)>0 or first_run==True:
        first_run=False
        skip=[]
        for j in range(len(geodesics_enclosing_I)):
            gamma=geodesics_enclosing_I[j]
            for side in polygon:
                # Check if gamma and side intersect, and update the vertices of both geodesics if so
                gamma.intersect_polygon_side(side,timeout=timeout,prec=prec)
            if gamma.is_polygon_side==True:
                skip.append(j)
                FreeSides=difference(FreeSides,[-infinity,gamma.foot_left],timeout=timeout,prec=prec)
                FreeSides=difference(FreeSides,[gamma.foot_right,infinity],timeout=timeout,prec=prec)
                polygon.append(gamma)
        skip.reverse()
        for j in skip:
            geodesics_enclosing_I.pop(j)

        # Here we find geodesics gamma exposing I that do not intersect an already existing side of polygon, but such that the interval between the feet of gamma is contained in a free side of the polygon, and this interval is maximal in the sense that for all other gamma' satisfying this condition, the interval beween the feet of gamma' does not contain the interval between the feet of gamma
        FreeSides_copy=deepcopy(FreeSides)
        for J in FreeSides_copy:
            lower_bound=J[0].n(prec)
            # Find all geodesics gamma whose left foot is at the lower bound of the interval J of FreeSides_copy
            gammas_with_left_foot_at_lower_bound=[]
            for gamma in geodesics_exposing_I:
                if gamma.foot_left.n(prec)==lower_bound:
                    gammas_with_left_foot_at_lower_bound.append(gamma)
            if len(gammas_with_left_foot_at_lower_bound)>0:
                # Sort these by right foot
                gammas_with_left_foot_at_lower_bound.sort(key=lambda x:x.foot_right.n(prec))
                # The geodesic with largest right foot is maximal in the sense described above and will contribute to the side of polygon
                gamma=gammas_with_left_foot_at_lower_bound[-1]
                # We make sure this geodesic doesn't intersect a side of polygon; if it does, it will be added to the polygon below, as we'll need to update vertices of geodesics contributing to the polgyon
                if gamma.foot_right.n(prec)<=J[1].n(prec):
                    FreeSides=difference(FreeSides,[gamma.foot_left,gamma.foot_right],timeout=timeout,prec=prec)
                    gamma.is_polygon_side=True
                    polygon.append(gamma)
                    geodesics_exposing_I.remove(gamma)                  

        skip=[]
        for j in range(len(geodesics_exposing_I)):
            gamma=geodesics_exposing_I[j]
            for side in polygon:
                # Check if gamma and side intersect, and update the vertices of both geodesics if so
                gamma.intersect_polygon_side(side,timeout=timeout,prec=prec)
            if gamma.is_polygon_side==True:
                skip.append(j)
                FreeSides=difference(FreeSides,[gamma.foot_left,gamma.foot_right],timeout=timeout,prec=prec)
                polygon.append(gamma)
        skip.reverse()
        for j in skip:
            geodesics_exposing_I.pop(j)
    return polygon


# Calculates the volume of a hyperbolic polygon, where the polygon is a list of geodesics comprising the sides of the polygon
def PolygonVolume(polygon,prec=53):
    # For each geodesic gamma contributing to a side of the polygon, we find the unit tangent vectors [[t0x,t0y],[t1x,t1y]] (directed toward the interior of the side of the polygon) at the counterclockwise- and clockwise-most vertices, respectively, of gamma and append this list to the list describing gamma; if the vertex is at infinity we let the tangent vector be [0,-1]
    for gamma in polygon:
        # We first consider vertical geodesics
        if gamma.type=='Vertical':
            # If gamma lies left of the imaginary axis, then the counterclockwise tangent vector will point straight up and the clockwise will point straight down
            if gamma.foot<0:
                gamma.tangent_vector_counterclockwise=[0,1]
                gamma.tangent_vector_clockwise=[0,-1]
            # If gamma lies right of the imaginary axis, then the counterclockwise vector will point straight down and the clockwise will point straight up
            else:
                gamma.tangent_vector_counterclockwise=[0,-1]
                gamma.tangent_vector_clockwise=[0,1]
        else:
            # Counterclockwise vertex
            # If the counterclockwise vertex is on the real line, then the corresponding tangent vector points straight up
            if gamma.vertex_counterclockwise[1]==0:
                gamma.tangent_vector_counterclockwise=[0,1]
            # Otherwise we find the slope of the unit tangent vector at the counterclockwise vertex and compute its real and imaginary parts
            else:
                tangent_slope=(gamma.center-gamma.vertex_counterclockwise[0])/(gamma.radius^2-(gamma.vertex_counterclockwise[0]-gamma.center)^2)^(1/2)
                normalizer=(tangent_slope^2+1)^(1/2)
                # Note that these tangent vectors are oriented so that they point toward the interior of the side of the polygon
                if gamma.type=='Encloses I':
                    gamma.tangent_vector_counterclockwise=[1/normalizer,tangent_slope/normalizer]
                elif gamma.type=='Exposes I':
                    gamma.tangent_vector_counterclockwise=[-1/normalizer,-tangent_slope/normalizer]
            # Clockwise vertex
            # If the clockwise vertex is on the real line, then the corresponding tangent vector points straight up
            if gamma.vertex_clockwise[1]==0:
                gamma.tangent_vector_clockwise=[0,1]
            # Otherwise we find the slope of the unit tangent vector at the clockwise vertex and compute its real and imaginary parts
            else:
                tangent_slope=(gamma.center-gamma.vertex_clockwise[0])/(gamma.radius^2-(gamma.vertex_clockwise[0]-gamma.center)^2)^(1/2)
                normalizer=(tangent_slope^2+1)^(1/2)
                # Note that these tangent vectors are oriented so that they point toward the interior of the side of the polygon
                if gamma.type=='Encloses I':
                    gamma.tangent_vector_clockwise=[-1/normalizer,-tangent_slope/normalizer]
                elif gamma.type=='Exposes I':
                    gamma.tangent_vector_clockwise=[1/normalizer,tangent_slope/normalizer]                      

    # For each side of the polygon, we compute the internal angle of the counterclockwise-most vertex of this side
    # We will record the interior angles of each of the vertices of the polygon in a list
    interior_angles=[]
    for gamma0 in polygon:
        gamma0_tangent_vector_angle=ang(gamma0.tangent_vector_counterclockwise)
        # We treat the case when the vertex is at infinity separately so that we don't subract infinities; in this case we know the interior angle is zero
        if gamma0.vertex_counterclockwise[1]==infinity:
            interior_angles.append(0)
        else:
            # Find the other side of the polygon sharing the counterclockwise vertex of gamma0
            vertex_differences=[((gamma1.vertex_clockwise[0]-gamma0.vertex_counterclockwise[0])^2+(gamma1.vertex_clockwise[1]-gamma0.vertex_counterclockwise[1])^2) for gamma1 in polygon]
            j=vertex_differences.index(min(vertex_differences))
            gamma1=polygon[j]
            gamma1_tangent_vector_angle=ang(gamma1.tangent_vector_clockwise)
            # Find the difference in the angles of the tangent vectors corresponding to gamma0 and gamma1 and append it to interior_angles
            if gamma0_tangent_vector_angle>=gamma1_tangent_vector_angle:
                interior_angles.append(RLF(2*pi)*(gamma0_tangent_vector_angle-gamma1_tangent_vector_angle))
            else:
                interior_angles.append(RLF(2*pi)*(gamma0_tangent_vector_angle-gamma1_tangent_vector_angle+1))           
    num_vertices=len(polygon)
    # Application of Gauss-Bonnett Theorem
    volume=RLF(pi)*(num_vertices-2)-sum(interior_angles)
    return volume   

############################################################
############################################################









#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################










r''' 
ALGORITHM 2.1 from Sanderson 'Implementing Edwards's Algorithm for Computing the Veech Group of a Translation Surface' (with some improvements, e.g. Voronoi staples).
Takes as input any compact translation surface (X,omega) from the flatsurf package; outputs generators of the Veech group of (X,omega) if said group is a lattice; otherwise outputs all elements of Veech group whose Frobenius norms are bounded above by NormBound.
timeout is the maximum number of seconds we spend on exact comparisons of numbers before switching to approximate comparisons.
prec is the number of bits of precision if approximate comparisons have to be made.
'''

def veech_group(X,NormBound=32,timeout=1,prec=53):
    from flatsurf.geometry.polygon import polygons,PolygonPosition
    # Ensure we have a compact translation surface
    assert isinstance(X,TranslationSurface)
    assert X.is_finite() is True

    # Make sure the polygons are labeled beginning with 0
    poly_labels=[l for l in X.label_iterator()]
    poly_labels.sort()
    new_poly_labels=range(X.num_polygons())
    if poly_labels!=new_poly_labels:
        X=X.relabel({poly_labels[j]:new_poly_labels[j] for j in range(X.num_polygons())})[0]
    X0=deepcopy(X)  

    ##############################################################
    ##############################################################
    ## 1. Let $n:=1$, $M_0:=I_2$ and $TrivialStabilizer:=False$ ##
    ##############################################################
    ##############################################################

    # Apply matrices M0 to X until a surface whose Veech group trivially stabilizes i in the upper-half plane is found
    n=1
    M0=matrix([[1,0],[0,1]])
    TrivialStabilizer=False

    #########################################
    #########################################
    ## 2. While $TrivialStabilizer=False$: ##
    #########################################
    #########################################

    while TrivialStabilizer==False:

        #############################################################################
        ## 2.i If $n>1$, let $M_0:=\begin{pmatrix} 1 & 0 \\ 1/n & 0 \end{pmatrix}$ ##
        #############################################################################
        
        if n>1:
            M0=matrix([[1,0],[1/n,1]])

        ####################################################
        ## 2.ii Let $(X,\omega):=M_0\cdot (X_0,\omega_0)$ ##
        ####################################################

        X=M0*X0

        ##########################################################################################################
        ## 2.iii  Calculate Voronoi staples of $(X,\omega)$                                                     ##
        ## 2.iv Calculate the lengths of the longest Voronoi staples corresponding to each distinct cone angle  ##
        ##########################################################################################################

        # Finds vertex equivalence classes of X 
        vertex_equivalence_classes=[]
        while len(vertex_equivalence_classes)<X.num_singularities():
            for poly in range(X.num_polygons()):
                for vertex in range(X.polygon(poly).num_edges()):
                    vertex_equivalence_class_set=X.singularity(poly,vertex).vertex_set()
                    vertex_equivalence_class=[v for v in vertex_equivalence_class_set]
                    vertex_equivalence_class.sort()
                    already_accounted=False
                    for equiv_class in vertex_equivalence_classes:
                        if equiv_class==vertex_equivalence_class:
                            already_accounted=True
                            break
                    if already_accounted==False:
                        vertex_equivalence_classes.append(vertex_equivalence_class) 

        # This records information to determine which disjoint copy of the plane a marked segment should lie in
        # vertices_in_disjoint_planes[P][copy_of_plane] gives a list of vertex data (polygon,vertex) such that separatrices leaving these vertices should develop in the copy_of_plane^th copy of the plane corresponding to singularity P
        # There is ambiguity for some vertices, such as vertex 6 of the regular octagon...these ambiguous vertices where separatrices may develop in either the copy_of_plane^th or (copy_of_plane+1)^st plane lead the list corresponding to the (copy_of_plane+1)^st plane
        X.vertices_in_disjoint_planes=[]
        first_edge_from_singularity=[]
        for P in range(X.num_singularities()):
            # Cone angle of singularity P
            cone_angle_P=sum([numerical_approx(X.polygon(vertex_data[0]).angle(vertex_data[1])) for vertex_data in vertex_equivalence_classes[P]])
            copy_of_plane=0

            # Will contain lists where the (copy_of_plane)^th list contains vertex data (polygon,vertex) with vertices corresponding to the (copy_of_plane)^th copy of the plane that corresponds to singularity P
            vertices_in_disjoint_planes_P=[]

            # Find the first vertex_data in the equivalence class of P such that e1=(1,0) is in the interior of the polgon at vertex_data
            e1_in_interior=False
            i=0
            while e1_in_interior==False:
                vertex_data=vertex_equivalence_classes[P][i]
                right_most_edge=X.polygon(vertex_data[0]).edge(vertex_data[1])
                left_most_edge=-X.polygon(vertex_data[0]).edge((vertex_data[1]-1)%X.polygon(vertex_data[0]).num_edges())
                if ang(right_most_edge)==0 or (ang(right_most_edge)>ang(left_most_edge) and ang(left_most_edge)>0):
                    e1_in_interior=True
                i+=1            
            P_angle=ang(left_most_edge)

            # Will be the (copy_of_plane)^th list inside vertices_in_disjoint_planes_P
            vertices_in_disjoint_planes_P_copy_of_plane=[vertex_data]

            while copy_of_plane<cone_angle_P:
                vertex_data=X.opposite_edge(vertex_data[0],(vertex_data[1]-1)%X.polygon(vertex_data[0]).num_edges())
                P_angle+=X.polygon(vertex_data[0]).angle(vertex_data[1])
                if P_angle>copy_of_plane+1:
                    vertices_in_disjoint_planes_P.append(vertices_in_disjoint_planes_P_copy_of_plane)
                    vertices_in_disjoint_planes_P_copy_of_plane=[]
                    copy_of_plane+=1
                vertices_in_disjoint_planes_P_copy_of_plane.append(vertex_data)
            X.vertices_in_disjoint_planes.append(vertices_in_disjoint_planes_P)
        # Sort by cone angle
        X.vertices_in_disjoint_planes.sort(key=lambda x:len(x))

        # Cone angles of X in increasing order
        X.cone_angles=[round(angle.n()) for angle in X.angles()]
        X.cone_angles.sort()

        # Here we construct groups of permutations of singularities and cylclic groups of copies of the plane associated to a particular singularity; these will comprise the group Trans(O)
        # Records number of singularities of a particular order; will create symmetric groups on these numbers of elements and later take cartesian product of these groups for various permutations of singularities
        X.permutation_group_orders=[]
        angle_previous=None
        for i in range(len(X.cone_angles)):
            angle=X.cone_angles[i]
            if angle!=angle_previous:
                X.permutation_group_orders.append(X.cone_angles.count(angle))
            angle_previous=angle
        # Cartesian product of symmetric groups of orders from permutation_group_orders; will be used to permute singularities of the same order
        symm_groups=[]
        for i in X.permutation_group_orders:
            symm_groups.append(SymmetricGroup(i))
        X.permutations_of_singularities=cartesian_product(symm_groups)
        # Cartesian product of cyclic groups from orders of cone angles; will be used to check cyclic permutations of copies of the plane associated to a particular singularity
        cyc_groups=[]
        for j in X.cone_angles:
            cyc_groups.append(CyclicPermutationGroup(j))
        X.permutations_of_copies_of_plane=cartesian_product(cyc_groups)

        # Voronoi staples of X as determined by the Voronoi decomposition
        X.VoronoiStaples=Staples(X)

        # Records the max length of the Voronoi staples in each copy of the plane corresponding to each singularity P
        VoronoiStaples_max_lengths=[[max([sqrt((vs[0][0])^2+(vs[0][1])^2) for vs in X.VoronoiStaples[P][copy_of_plane]]) for copy_of_plane in range(len(X.VoronoiStaples[P]))] for P in range(X.num_singularities())]
        # For each distinct cone angle, we find the longest Voronoi staple emenating from a singularity of that particular cone angle
        previous_cone_angle=X.cone_angles[0]
        # Will record the length of longest Voronoi staple in each copy of the plane corresponding to all singularities of a particular cone angle
        max_length_in_each_plane=[]
        # List of lengths of longest Voronoi staple from a singularity of particular cone angle (in order of X.cone_angles), i.e. will record the max of max_length_in_each_plane for each distinct cone angle
        longest_VoronoiStaples=[]
        for P in range(X.num_singularities()):
            if X.cone_angles[P]==previous_cone_angle:
                max_length_in_each_plane=max_length_in_each_plane+VoronoiStaples_max_lengths[P]
            else:
                longest_VoronoiStaples.append(max(max_length_in_each_plane))
                max_length_in_each_plane=VoronoiStaples_max_lengths[P]
            previous_cone_angle=X.cone_angles[P]
            if P==X.num_singularities()-1:
                longest_VoronoiStaples.append(max(max_length_in_each_plane))

        # List of radii to find marked segments within; if there are n singularities of the same cone angle, then the radius r from longest_VoronoiStaples corresponding to this cone angle will appear n times
        X.radii=[]
        j=0
        previous_cone_angle=X.cone_angles[0]
        for P in range(X.num_singularities()):
            if X.cone_angles[P]!=previous_cone_angle:
                j+=1
            X.radii.append(longest_VoronoiStaples[j])
            previous_cone_angle=X.cone_angles[P]
        a=sqrt(2)

        ####################################################
        ## 2.v Calculate $\text{SL}^{\sqrt{2}}(X,\omega)$ ##        
        ####################################################

        [SL_a,matrices_checked_radii,matrices_to_check]=SLa(X,a,PreviousSLa=[matrix.identity(2)],MatricesChecked=[matrix.identity(2)],MatricesToCheck=[],PreviousNormBound=None)

        ###############################################################################################
        ## 2.vi If $\text{SL}^{\sqrt{2}}(X,\omega)\subseteq\{\pm I_2\}$, let TrivialStabilizer:=True ##
        ## 2.vii Increase $n$ by 1                                                                   ##     
        ###############################################################################################

        if len(SL_a)>2 or (len(SL_a)==2 and -matrix.identity(2) not in SL_a):
            n+=1
        else:
            TrivialStabilizer=True      

    ################################################################################################################################
    ################################################################################################################################
    ## 3. If $-I_2\in\text{SL}^{\sqrt{2}}(X,\omega)$, let $ContainsMinusIdentity:=True$; else, let $ContainsMinusIdentity:=False$ ##
    ################################################################################################################################
    ################################################################################################################################

    if -matrix.identity(2) in SL_a:
        ContainsMinusIdentity=True
    else:
        ContainsMinusIdentity=False 

    #######################################
    #######################################
    ## 4. Let $GuaranteedLattice:=False$ ##
    #######################################
    #######################################

    GuaranteedLattice=False

    ##########################################
    ##########################################
    ## 5. While $GuaranteedLattice:=False$: ##
    ##########################################
    ##########################################

    while GuaranteedLattice==False:

        #################################
        ## 5.i Double the norm bound a ##
        #################################
        previous_a=deepcopy(a)
        if 2*a<NormBound:
            if a==sqrt(2):
                a=2
            else:
                a=2*a
        else:
            a=NormBound

        ####################################################
        ## 5.ii & 5.iii Calculate $\text{SL}^a(X,\omega)$ ##
        ####################################################

        [SL_a,matrices_checked_radii,matrices_to_check]=SLa(X,a,PreviousSLa=SL_a,MatricesChecked=matrices_checked_radii,MatricesToCheck=matrices_to_check,PreviousNormBound=previous_a)

        ##############################################################################################
        ## 5.iv Construct $\Omega(\text{PSL^a(X,\omega)})$                                          ##
        ## 5.v If $\Omega(\text{PSL^a(X,\omega)})$ has no free sides, let $GuaranteedLattice:=True$ ##
        ##############################################################################################

        PSL_a=[matrix([[1,0],[0,1]])]
        if ContainsMinusIdentity==True:
            for M in SL_a:
                if M not in PSL_a and -M not in PSL_a:
                    PSL_a.append(M)
        else:
            PSL_a=SL_a

        # Below we determine whether or not the hyperbolic polygon determined by by PSL_a has finite hyperbolic area by considering whether or not it has free sides (we need not actually construct the polygon here)

        # The perpendicular bisectors corresponding to each M*I for M in PSL_a
        bisectors=[geodesic(matrix=M) for M in PSL_a if M!=matrix([[1,0],[0,1]])]
        # # Picture of the perpendicular bisectors found thus far
        # sum([g.plot() for g in bisectors])

        # When the cardinality of free_sides is finite (from the construction of difference(), this is when free_sides==[]), we'll know the hyperbolic polygon has finite area and hence the group is guaranteed to be a lattice
        free_sides=[[-infinity,infinity]]

        for gamma in bisectors:
            if gamma.type=='Vertical':
                if gamma.foot<0:
                    free_sides=difference(free_sides,[-infinity,gamma.foot],timeout=timeout,prec=prec)
                else:
                    free_sides=difference(free_sides,[gamma.foot,infinity],timeout=timeout,prec=prec)
            elif gamma.type=='Encloses I':
                free_sides=difference(free_sides,[-infinity,gamma.foot_left],timeout=timeout,prec=prec)
                free_sides=difference(free_sides,[gamma.foot_right,infinity],timeout=timeout,prec=prec)
            elif gamma.type=='Exposes I':
                free_sides=difference(free_sides,[gamma.foot_left,gamma.foot_right],timeout=timeout,prec=prec)
            if len(free_sides)==0:
                GuaranteedLattice=True
                break


    #####################################################################################
    #####################################################################################
    ## 6. If $GuaranteedLattice=False$, return $\text{SL}^a(X,\omega)$; else, continue ##
    #####################################################################################
    #####################################################################################                

        # If X is not found to be a lattice and we've reached the NormBound, ask the user whether to double NormBound and try again or end the computation
        if GuaranteedLattice==False and a>=NormBound:
            print([M0.inverse()*M*M0 for M in SL_a])
            user_answer=input('Determination of a lattice is inconclusive; an exhaustive list of Veech group elements with Frobenius norm bounded above by {} is given above.  Double the norm bound and continue? Y/N:'.format(a))
            if user_answer=='Y':
                NormBound=2*NormBound
            else:
                return [M0.inverse()*M*M0 for M in SL_a]    


    ##################################################################################
    ##################################################################################
    ## 7. Let $BoundedRegion:=\Omega(\text{PSL^a(X,\omega)}) \cap B(i,\log(\nu(a))$ ##
    ##################################################################################
    ##################################################################################

    # We actually calculate this in the next step; for now we assign arbitrary volumes so as to enter the next while loop
    BoundedRegion_volume_lower_bound=0
    Omega_SL_a_volume=100


    ################################################################################
    ################################################################################
    ## 8. While $vol(BoundedRegion)\le (1/2)*vol(\Omega(\text{PSL^a(X,\omega)})): ##
    ################################################################################
    ################################################################################

    first_run=True
    AllSidesRepresented=False
    while AllSidesRepresented==False or BoundedRegion_volume_lower_bound.n(prec)<=Omega_SL_a_volume.n(prec)/2:

        ######################################################################
        ## 8.i, 8.ii & 8.iii Double a and calculate $\text{SL}^a(X,\omega)$ ##
        ######################################################################

        # We only double a and compute new Veech group elements if this is not the first run; the first run skips this and calculates the volume of BoundedRegion below
        SL_a_previous=copy(SL_a)
        if first_run==False:
            a=2*a
            [SL_a,matrices_checked_radii,matrices_to_check]=SLa(X,a,PreviousSLa=SL_a,MatricesChecked=matrices_checked_radii,MatricesToCheck=matrices_to_check,PreviousNormBound=a/2)
        first_run=False     


        ####################################################################################################################################
        ## 8.iv Construct $\Omega(\text{PSL^a(X,\omega)})$ and let $BoundedRegion:=\Omega(\text{PSL^a(X,\omega)}) \cap B(i,\log(\nu(a))$  ##
        ####################################################################################################################################

        # To be used to construct the new hyperbolic polygon determined by SL_a, which we denote Omega_SL_a
        PSL_a=[matrix([[1,0],[0,1]])]
        if ContainsMinusIdentity==True:
            for M in SL_a:
                if M not in PSL_a and -M not in PSL_a:
                    PSL_a.append(M)
        else:
            PSL_a=SL_a          

        # We need only append geodesics corresponding to newly found matrices to the list bisectors
        for M in PSL_a:
            if M not in SL_a_previous:
                bisectors.append(geodesic(matrix=M))
        # Picture of the perpendicular bisectors found thus far
        # sum([g.plot() for g in bisectors])

        # The hyperbolic polygon determined by SL_a
        Omega_SL_a=Omega(bisectors,prec=prec,timeout=timeout)
        # Picture of Omega_SL_a
        # sum([g.plot() for g in Omega_SL_a])

        BoundingRadius=log(nu(a))

        # A computation with the integral definition of hyperbolic distance gives the following imaginary part of the Euclidean center (the center lies on the imaginary axis) and Euclidean radius of the hyperbolic ball centered at I with radius BoundingRadius
        BoundingBall_euc_center=cosh(BoundingRadius).n(prec)
        BoundingBall_euc_radius=sinh(BoundingRadius).n(prec)
        # Picture of Omega_SL_a and bounding ball
        # sum([g.plot() for g in Omega_SL_a]+[circle((0,BoundingBall_euc_center),BoundingBall_euc_radius,edgecolor='red')])
        # sum([g.plot() for g in Omega_SL_a]+[circle((0,BoundingBall_euc_center),BoundingBall_euc_radius,edgecolor='orange',linestyle='-')])

        # We will fill this list with the geodesics contributing to Omega_SL_a and some additional geodesics that depend on the bounding ball; see comments after the following for loop
        bounded_polygon_geodesics=[]

        # We determine whether or not all sides of the polygon are represented (if they all intersect the bounding ball)
        # If any side is not represented, we'll make this flag false and break the loop
        represented_flag=True
        # We keep track of which sides of the polygon intersect the ball, and where the intersections on the boundary of the ball occur; each entry will be a list [[x,y],gamma,vertex_orientation], where x,y are the real and imaginary parts of this intersection; gamma is the side in question; and vertex_orientation='counterclockwise' if the segment of gamma that the bounding ball 'cuts off' countains the counterclockwise vertex of gamma, and vertex_orientation='clockwise' otherwise
        ball_intersections=[]
        for gamma in Omega_SL_a:
            bounded_polygon_geodesics.append(gamma)
            # Gamma is vertical
            if gamma.type=='Vertical':
                # Check if the h-line that gamma lies on intersects the boundary of the bounding ball; note that the inequalities here are strict so that gamma may not be tangent to the ball but must actually intersect it
                if gamma.foot>-BoundingBall_euc_radius and gamma.foot<BoundingBall_euc_radius:
                    # Now we find the (imaginary parts of the) intersections and check if the finite vertices of gamma lie strictly between these intersections
                    lower_intersection=(-sqrt(BoundingBall_euc_radius^2-gamma.foot^2)+BoundingBall_euc_center)
                    upper_intersection=(sqrt(BoundingBall_euc_radius^2-gamma.foot^2)+BoundingBall_euc_center)
                    min_vertex=min(gamma.vertex_counterclockwise[1],gamma.vertex_clockwise[1])
                    max_vertex=max(gamma.vertex_counterclockwise[1],gamma.vertex_clockwise[1])
                    # If the intersection is not in the interior of the segment bounded by the vertices of gamma, then gamma is not represented
                    if min_vertex>=upper_intersection or max_vertex<=lower_intersection:
                        represented_flag=False
                        break
                    # Otherwise gamma is represented, and we see if the bounding ball `cuts off' either vertex                          
                    else:
                        if gamma.foot<0:
                            lower_vertex_orientation='counterclockwise'
                            upper_vertex_orientation='clockwise'
                        else:
                            lower_vertex_orientation='clockwise'
                            upper_vertex_orientation='counterclockwise'     
                    if min_vertex<lower_intersection:
                        ball_intersections.append([[gamma.foot,lower_intersection],gamma,lower_vertex_orientation])
                    if max_vertex>upper_intersection:
                        ball_intersections.append([[gamma.foot,upper_intersection],gamma,upper_vertex_orientation])
                else: 
                    represented_flag=False
                    break
            # Gamma is not vertical
            else: 
                # Setting equal the equations describing the the boundary of the ball and the semicircle describing gamma, we find that these curves intersect at x satisfying Ax^2+Bx+C=0 for A,B,C below
                A=4*(gamma.center^2+BoundingBall_euc_center^2)
                B=4*(gamma.center*((BoundingBall_euc_center^2-BoundingBall_euc_radius^2)-(gamma.center^2-gamma.radius^2))-2*BoundingBall_euc_center^2*gamma.center)
                C=((BoundingBall_euc_center^2-BoundingBall_euc_radius^2)-(gamma.center^2-gamma.radius^2))^2+4*BoundingBall_euc_center^2*(gamma.center^2-gamma.radius^2)
                # If the discriminant of this polynomial is positive, then the h-line corresponding to gamma intersects the ball (it must intersect the boundary of the ball twice, hence the inequality is strict); otherwise the side is not represented by the ball
                if B^2-4*A*C>0:
                    # The quadratic formula gives the real parts of the intersections
                    left_intersection=(-B-sqrt(B^2-4*A*C))/(2*A)
                    right_intersection=(-B+sqrt(B^2-4*A*C))/(2*A)
                    # We use these to find the imaginary parts of the intersections
                    left_intersection_imgy_pt=sqrt(gamma.radius^2-(left_intersection-gamma.center)^2)
                    right_intersection_imgy_pt=sqrt(gamma.radius^2-(right_intersection-gamma.center)^2)
                    # Real parts of the left and right vertices of gamma
                    left_vertex=min(gamma.vertex_counterclockwise[0],gamma.vertex_clockwise[0])
                    right_vertex=max(gamma.vertex_counterclockwise[0],gamma.vertex_clockwise[0])
                    # If the intersection is not in the interior of the segment bounded by the vertices of gamma, then gamma is not represented
                    if left_vertex>=right_intersection or right_vertex<=left_intersection:
                        represented_flag=False
                        break
                    # Otherwise gamma is represented, and we see if the bounding ball `cuts off' either vertex                          
                    else:
                        if gamma.type=='Encloses I':
                            left_vertex_orientation='counterclockwise'
                            right_vertex_orientation='clockwise'
                        elif gamma.type=='Exposes I':
                            left_vertex_orientation='clockwise'
                            right_vertex_orientation='counterclockwise'
                        if left_vertex<left_intersection:
                            ball_intersections.append([[left_intersection,left_intersection_imgy_pt],gamma,left_vertex_orientation])
                        if right_vertex>right_intersection:
                            ball_intersections.append([[right_intersection,right_intersection_imgy_pt],gamma,right_vertex_orientation])
                else:
                    represented_flag=False
                    break
        # If our flag is still true, then all sides of the polygon are represented by the bounding ball
        if represented_flag==True:
            AllSidesRepresented=True                        

        # If AllSidesRepresented is True, then for each side of the polygon, find where (if at all) the boundary of the bounding ball intersects the interior of the side...
        # If it intersects two sides that share a vertex, find the h-line between these two intersections...
        # This h-line lies inside the closure of the ball since a Euclidean ball in the upper-half plane is h-convex...
        # Find the area of the hyperbolic polygon formed as the region enclosed by the previous sides of the polygon plus these new h-lines
        # This gives a lower bound on the hyperbolic area of Omega_SL_a intersected with the ball

        # Note that we only need to compute this if BoundedRegion_volume_lower_bound<=Omega_SL_a_volume/2; once BoundedRegion_volume_lower_bound>Omega_SL_a_volume/2 for some r, then this is true for all larger r' since the bounding ball is growing and the Dirichlet polygon Omega_SL_a is shrinking as r increases

        if AllSidesRepresented==True and BoundedRegion_volume_lower_bound.n(prec)<=Omega_SL_a_volume.n(prec)/2:
            # Begin by finding h-lines passing through points of ball_intersections corresponding to sides that share a vertex
            skip=[]
            for j in range(len(ball_intersections)):
                if j not in skip:
                    intersection0=ball_intersections[j]
                    skip.append(j)
                    gamma0=intersection0[1]
                    # Find the vertex of gamma0 that the bounding ball 'cuts off'
                    if intersection0[-1]=='counterclockwise':
                        vertex=gamma0.vertex_counterclockwise
                    else:
                        vertex=gamma0.vertex_clockwise
                    # Find which other (distinct) intersection also corresponds to the bounding ball 'cutting off' vertex
                    for k in range(len(ball_intersections)):
                        if k not in skip:
                            intersection1=ball_intersections[k]
                            gamma1=intersection1[1]
                            ### LOOK INTO THIS; THE COMPARISONS WITH 0 DON'T ALWAYS WORK, E.G. FOR X=mcmullen_L(1,1,1,1) ###
                            if (intersection0[-1]=='counterclockwise' and (((gamma1.vertex_clockwise[0]-vertex[0]).n(prec)==0.n(prec) and (gamma1.vertex_clockwise[1]-vertex[1]).n(prec)==0.n(prec)) or (gamma1.vertex_clockwise[1]==infinity and vertex[1]==infinity))) or (intersection0[-1]=='clockwise' and (((gamma1.vertex_counterclockwise[0]-vertex[0]).n(prec)==0.n(prec) and (gamma1.vertex_counterclockwise[1]-vertex[1]).n(prec)==0.n(prec)) or (gamma1.vertex_counterclockwise[1]==infinity and vertex[1]==infinity))):
                                skip.append(k)
                                # Below, we construct the h-line between these two intersections and append it to the sides comprising the original polygon
                                intersection0_x=intersection0[0][0]
                                intersection0_y=intersection0[0][1]
                                intersection1_x=intersection1[0][0]
                                intersection1_y=intersection1[0][1]
                                # If the real parts of the two intersections are the same, then the h-line is vertical; we don't bother with finding the actual vertices here since the function Omega() will automatically find them
                                if intersection0_x==intersection1_x:
                                    h_line=geodesic(foot=intersection0_x)
                                # Otherwise we find the left and right feet of the h-line passing through both intersections; again we don't bother with finding the actual vertices here since the function Omega() will automatically find them
                                else:
                                    # Center and radius of the h-line between intersection0 and intersection1
                                    center=(intersection0_y^2+intersection0_x^2-(intersection1_y^2+intersection1_x^2))/(2*(intersection0_x-intersection1_x))
                                    radius=sqrt(intersection0_y^2+(intersection0_x-center)^2)
                                    # Left and right feet of the h-line
                                    left=(center-radius)
                                    right=(center+radius)
                                    h_line=geodesic(foot_left=left,foot_right=right)
                                # By default we set the minimum distance from I to the h-line equal to infinity as a specific value isn't needed 
                                h_line.distance_from_I=infinity
                                bounded_polygon_geodesics.append(h_line)
            Omega_SL_a_bounding_ball=Omega(bounded_polygon_geodesics,prec=prec,timeout=timeout)
            # Picture of Omega_SL_a_bounding_ball and bounding ball
            # sum([g.plot() for g in Omega_SL_a_bounding_ball]+[circle((0,BoundingBall_euc_center),BoundingBall_euc_radius,edgecolor='red')])

            BoundedRegion_volume_lower_bound=PolygonVolume(Omega_SL_a_bounding_ball,prec=prec)
            Omega_SL_a_volume=PolygonVolume(Omega_SL_a,prec=prec)


    ######################################################################################################################################
    ######################################################################################################################################
    ## 9. Let $SideRepresentatives$ be the elements of $\text{PSL^a(X,\omega)}$ that form the sides of $\Omega(\text{PSL^a(X,\omega)})$ ##
    ######################################################################################################################################
    ######################################################################################################################################

    # List of all matrices corresponding to sides of Omega_SL_a
    all_matrices=[gamma.matrix for gamma in Omega_SL_a]
    # We don't need to include a matrix if its inverse is already include, so we get rid of these non-essential matrices
    SideRepresentatives=[]
    for M in all_matrices:
        if M.inverse() not in SideRepresentatives and -M.inverse() not in SideRepresentatives:
            SideRepresentatives.append(M)   

    ###########################################################################################################################################################
    ###########################################################################################################################################################
    ## 10. If $ContainsMinusIdentity=True$, let $AssociatedGenerators:=SideRepresentatives\cup\{-I\}$; else, let $AssociatedGenerators:=SideRepresentatives$ ##
    ###########################################################################################################################################################
    ###########################################################################################################################################################

    AssociatedGenerators=SideRepresentatives
    if ContainsMinusIdentity==True:
        AssociatedGenerators.append(-matrix.identity(2))

    #######################################################################
    #######################################################################
    ## 11. Let $Generators:=M_0^{-1}\cdot AssociatedGenerators\cdot M_0$ ##
    #######################################################################
    #######################################################################

    Generators=[M0.inverse()*M*M0 for M in AssociatedGenerators]        

    #############################
    #############################
    ## 12. Return $Generators$ ##
    #############################
    #############################

    print('The Veech group is a lattice generated by the following matrices:')
    return Generators


























