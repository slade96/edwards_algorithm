############################################################
############################################################
##  Classes and functions used within veech_group() below ##
############################################################
############################################################

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
			self.distance_from_I=chi_2(frobenius_norm)
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


def chi_1_inv(x):
	return sqrt(x^2+x^(-2))

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

# Takes as its input the Frobenius norm of a matrix M and outputs the minimum distance between I and the perpendicular bisector between I and M*I in the upper-half plane
def chi_2(x):
	return -ln(sqrt((1/2)*(x^2-sqrt(x^4-4))))

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


# Calculates the volume of a hyperbolic polygon, where the polygon is a list of geodesics cmprising the sides of the polygon
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
ALGORITHM 7.2
Takes as input any compact translation surface (X,omega) with rational side lengths from the flatsurf package with polygon labels 0,1,...,n.
iteration_limit is the number of times we double the value r in step 8 of Algorithm 7.1
timeout is the maximum number of seconds we spend on exact comparisons of numbers before switching to approximate comparisons
prec is the number of bits of precision if approximate comparisons have to be made
'''

def veech_group(X,iteration_limit=10,timeout=1,prec=53):

	from flatsurf.geometry.polygon import polygons,PolygonPosition
	assert isinstance(X,TranslationSurface)
	assert X.is_finite() is True

	# Make sure the polygons are labeled beginning with 0
	num_polygons=X.num_polygons()
	poly_labels=[l for l in X.label_iterator()]
	poly_labels.sort()
	new_poly_labels=range(num_polygons)
	if poly_labels!=new_poly_labels:
		X=X.relabel({poly_labels[j]:new_poly_labels[j] for j in range(num_polygons)})[0]

	# We iterate the first five steps of the algorithm until we find a new surface whose Veech group trivially stabilizes I in the upper-half plane (see step 5)
	n=1
	trivial_stabilizer=False
	while trivial_stabilizer==False:
		if n==1:
			M0=matrix([[1,0],[0,1]])
		else:
			M0=matrix([[1,0],[1/n,1]])
		X=M0*X



		###################################################
		###################################################
		## 7.2.0: Compute the Voronoi decomposition of X ##
		###################################################
		###################################################

		# Finds vertex equivalence classes of X
		num_singularities=X.num_singularities()
		vertex_equivalence_classes=[]
		vertex_equivalence_classes_sets=[]
		while len(vertex_equivalence_classes)<num_singularities:
			for polygon in range(num_polygons):
				for vertex in range(X.polygon(polygon).num_edges()):
					vertex_is_equivalent_to=X.singularity(polygon,vertex).vertex_set()
					already_accounted=False
					for i in vertex_equivalence_classes_sets:
						if i==vertex_is_equivalent_to:
							already_accounted=True
							break
					if already_accounted==False:
						vertex_equivalence_classes_sets.append(vertex_is_equivalent_to)
						vertex_equivalence_classes.append([v for v in vertex_is_equivalent_to])

		# Gives the vertices of each polygon of X in AA (algebraic reals) so that the Voronoi decomposition can be computed with respect to these points
		vertices=[[[AA(c) for c in vert] for vert in X.polygon(poly).vertices()] for poly in range(num_polygons)]
		# Gives a list whose entries are the Voronoi decompositions with respect to the vertices of each polygon
		VD=[VoronoiDiagram(vertices[poly]) for poly in range(num_polygons)]
		# A list of lists; each entry of the outer list corresponds to a polygon of X; each entry of the inner lists is the (potentially) unbounded polyhedron in the Voronoi decomposition associated to a particular vertex of said polygon
		VD_polyhedra_unbdd=[[VD[poly].regions()[VD[poly].points()[v]] for v in range(X.polygon(poly).num_edges())] for poly in range(num_polygons)]
		# Same structure as the list above, but now the polyhedra are intersected with the proper polygon of X, so they are bounded
		VD_polyhedra=[[VD_polyhedra_unbdd[poly][vert].intersection(Polyhedron(vertices=vertices[poly])) for vert in range(X.polygon(poly).num_edges())] for poly in range(num_polygons)]

		###################################################
		###################################################










		######################################################################
		######################################################################
		## 7.2.1: Calculate $\rho=\rho(X,\omega)$ using the Voronoi 2-cells ##
		######################################################################
		######################################################################

		# Creates a list of the max distances in each polygon of the Voronoi decomposition from the singularity of X (main_vertex below) to the boundary of the polygon
		distances_from_singularities=[]
		for poly in range(num_polygons):
			for vert in range(X.polygon(poly).num_edges()):
				Voro_poly=VD_polyhedra[poly][vert]
				Voro_poly_verts=Voro_poly.vertices_list()
				Voro_poly_num_edges=len(Voro_poly_verts)

				# Gives the singularity of X (called main_vertex) that the polygon Voro_poly was constructed about
				vertices_tuple=[tuple(v) for v in vertices[poly]]
				Voro_poly_verts_tuple=[tuple(v) for v in Voro_poly_verts]
				main_vertex=[list(set(vertices_tuple).intersection(Voro_poly_verts_tuple))[0][0],list(set(vertices_tuple).intersection(Voro_poly_verts_tuple))[0][1]]

				distances_from_main_vertex=[sqrt((Voro_poly_verts[v][0]-main_vertex[0])^2+(Voro_poly_verts[v][1]-main_vertex[1])^2) for v in range(Voro_poly_num_edges)]
				distances_from_singularities.append(max(distances_from_main_vertex))
		rho=max(distances_from_singularities)

		##################################################################
		##################################################################









		####################################################################
		####################################################################
		## 7.2.2: Let $r=2\rho$ and let $b=\chi_1^{-1}(2\rho/r)=\sqrt{2}$ ##
		####################################################################
		####################################################################

		r=2*rho
		b=sqrt(2)

		####################################################################
		####################################################################










		##############################################################
		##############################################################
		## 7.2.3: Calculate $MP_P^r(X,\omega)$ for all $P\in\Sigma$ ##
		##############################################################
		##############################################################

		# This records information to determine which disjoint copy of the plane an element of $MP_P^r(X,\omega)$ should lie in
		# vertices_in_disjoint_planes[P][n] gives a list of vertex data (polygon,vertex) such that separatrices leaving these vertices should develop in the n^th copy of the plane corresponding to singularity P
		# There is ambiguity for some vertices, such as vertex 6 of the regular octagon...these ambiguous vertices where separatrices may develop in either the n^th or (n+1)^st plane lead the list corresponding to the (n+1)^st plane
		vertices_in_disjoint_planes=[]
		first_edge_from_singularity=[]
		cone_angles=[]
		for P in range(num_singularities):
			vertex_angles=[X.polygon(vertex_data[0]).angle(vertex_data[1]).numerical_approx() for vertex_data in vertex_equivalence_classes[P]]
			cone_angles.append(round(sum(vertex_angles)))
			copy_of_plane=0
			
			# Will contain lists where the (copy_of_plane)^th list contains vertex data (polygon,vertex) with vertices corresponding to the (copy_of_plane)^th copy of the plane that corresponds to singularity P
			vertices_in_disjoint_planes_P=[]

			# Canonicalizes the set of vertices corresponding to singularity P by vertex data (polygon,vertex): sorts by lowest polygon index then lowest vertex index
			singularity_P=sorted(vertex_equivalence_classes[P])

			# Ensures that the vector e1=(1,0) is in the interior of the polgon at vertex_data
			e1_in_interior=False
			i=0
			while e1_in_interior==False:
				vertex_data=singularity_P[i]
				right_most_edge=X.polygon(vertex_data[0]).edge(vertex_data[1])
				left_most_edge=-X.polygon(vertex_data[0]).edge((vertex_data[1]-1)%X.polygon(vertex_data[0]).num_edges())
				if ang(right_most_edge)==0 or (ang(right_most_edge)>ang(left_most_edge) and ang(left_most_edge)>0):
					e1_in_interior=True
				i+=1			

			P_angle=ang(left_most_edge)

			# Will be the (copy_of_plane)^th list inside vertices_in_disjoint_planes_P
			vertices_in_disjoint_planes_P_copy_of_plane=[vertex_data]

			while copy_of_plane<cone_angles[P]:
				vertex_data=X.opposite_edge(vertex_data[0],(vertex_data[1]-1)%X.polygon(vertex_data[0]).num_edges())
				P_angle+=X.polygon(vertex_data[0]).angle(vertex_data[1])
				if P_angle>copy_of_plane+1:
					vertices_in_disjoint_planes_P.append(vertices_in_disjoint_planes_P_copy_of_plane)
					vertices_in_disjoint_planes_P_copy_of_plane=[]
					copy_of_plane+=1
				vertices_in_disjoint_planes_P_copy_of_plane.append(vertex_data)
			vertices_in_disjoint_planes.append(vertices_in_disjoint_planes_P)



		### To do: create an optional parameter to add to previously computed marked_periods(radius0) ###
		def marked_periods(radius):
			MPr=[]
			for P in range(num_singularities):
				MPr_P=[[] for i in range(cone_angles[P])]
				for copy_of_plane in range(len(vertices_in_disjoint_planes[P])):
					for vertex_data in vertices_in_disjoint_planes[P][copy_of_plane]:			
						polygon_start=vertex_data[0]
						vertex_start=vertex_data[1]
						saddle_connections=X.saddle_connections(radius^2,polygon_start,vertex_start)
						saddle_connections_directions=[sc.direction() for sc in saddle_connections]
						# First see if the current vertex is ambiguous (in the sense that saddle connections from it may be assigned to the current copy of the plane or the previous copy)
						if vertex_data==vertices_in_disjoint_planes[P][copy_of_plane][0]:
							right_most_edge=X.polygon(polygon_start).edge(vertex_start)
							left_most_edge=-X.polygon(polygon_start).edge((vertex_start-1)%X.polygon(polygon_start).num_edges())
							right_most_edge_angle=ang(right_most_edge)
							left_most_edge_angle=ang(left_most_edge)
							for i in range(len(saddle_connections)):
								sc=saddle_connections[i]
								holonomy=sc.holonomy()
								holonomy_angle=ang(holonomy)
								# If sc leaves its singularity clockwise of the vector (1,0), then its holonomy vector is assigned to the previous copy of the plane
								if right_most_edge_angle>left_most_edge_angle and holonomy_angle>=right_most_edge_angle:
									cumulative_angle=(copy_of_plane-1)%cone_angles[P]+holonomy_angle
									MPr_P[copy_of_plane-1].append([holonomy,sc,cumulative_angle])
								# Otherwise its holonomy vector is assigned to the current copy of the plane
								else:
									cumulative_angle=copy_of_plane+holonomy_angle
									MPr_P[copy_of_plane].append([holonomy,sc,cumulative_angle])
						# These are the unambiguous vertices
						else:	
							for i in range(len(saddle_connections)):
								sc=saddle_connections[i]
								holonomy=sc.holonomy()
								holonomy_angle=ang(holonomy)
								cumulative_angle=copy_of_plane+holonomy_angle
								MPr_P[copy_of_plane].append([holonomy,sc,cumulative_angle])
				MPr.append(MPr_P)
				# Sort holonomies in each copy of the plane by cumulative counterclockwise angle from (1,0) on 0th copy of plane
				for copy_of_plane in range(len(MPr[P])):
					MPr[P][copy_of_plane].sort(key=lambda x:x[2])
			return MPr

		MP_2rho=marked_periods(2*rho)
		# Sort MP_2rho by cone angles of singularities so we may easily check this list against (permutations of) F_M_MPr_bounded_by_2rho below
		MP_2rho_sorted=sorted(MP_2rho,key=len)
		cone_angles_sorted=sorted(cone_angles)

		# This is the same as MP_2rho_sorted, but only with the holonomy information (i.e., we drop the saddle connection data and angle data); this will be used in the function A_r() below
		MP_2rho_sorted_holonomies=[[[MP_2rho_sorted[P][copy_of_plane][i][0] for i in range(len(MP_2rho_sorted[P][copy_of_plane]))] for copy_of_plane in range(len(MP_2rho_sorted[P]))] for P in range(num_singularities)]
		# This is just a list of the saddle connections from MP_2rho_sorted and is not broken up into sub-lists corresponding to singularities or copies of the plane; this again will be used in the function A_r() below to check the Z2-action
		MP_2rho_sorted_saddle_connections=[MP_2rho_sorted[P][copy_of_plane][i][1] for P in range(num_singularities) for copy_of_plane in range(len(MP_2rho_sorted[P])) for i in range(len(MP_2rho_sorted[P][copy_of_plane]))]
		num_holonomies_in_MP_2rho=len(MP_2rho_sorted_saddle_connections)

		##############################################################
		##############################################################










		#################################################################################################################################
		#################################################################################################################################
		## 7.2.4: Calculate $A_r=\{M\in\Gamma(X,\omega)\ |\ ||M||\le b\}=\text{SO}(2,\mathbb{R})\cap\Gamma(X,\omega)$ using Theorem 18 ##
		#################################################################################################################################
		#################################################################################################################################

		def A_r(radius,norm_bound,previous_Ar=[matrix([[1,0],[0,1]])],matrices_checked=[matrix([[1,0],[0,1]])]):	
			# matrices_to_check will be a list comprised of matrices in $\text{SL}_2\mathbb{R}$ whose Frobenius norms are bounded above by norm_bound and whose inverses take a basis of vectors in marked_periods(2*rho) to another basis of vectors in marked_periods(radius)
			matrices_to_check=deepcopy(matrices_checked)
			MPr=marked_periods(radius)
			MPr_sorted=sorted(MPr,key=len)

			# Basis from projection of MP_2rho onto complex plane
			v0=MP_2rho_sorted[0][0][0][0]
			v1=MP_2rho_sorted[0][0][1][0]
			T=matrix([[v0[0],v1[0]],[v0[1],v1[1]]])
			# Note v0,v1 correspond to a singularity of minimal cone angle; below we will find another pair of marked periods that can be sent to v0,v1; these must correspond to a singularity of minimal cone angle as well
			num_sings_of_min_angle=cone_angles_sorted.count(cone_angles_sorted[0])
			for P in range(num_sings_of_min_angle):
				all_projections_to_plane=[[MPr_sorted[P][copy_of_plane][mp][0] for mp in range(len(MPr_sorted[P][copy_of_plane]))] for copy_of_plane in range(len(MPr_sorted[P]))]
				for copy_of_plane in range(len(MPr_sorted[P])):
					# We consider only projected marked periods from (copy_of_plane)^th and adjacent copies of the plane onto the comlex plane and choose our new basis from these holonomies
					specific_projections_to_plane=[]
					for i in range(3):
						specific_projections_to_plane.extend(all_projections_to_plane[(copy_of_plane-1+i)%len(MPr_sorted[P])])
					unique_projections_to_plane=[]
					for holonomy in specific_projections_to_plane:
						if holonomy not in unique_projections_to_plane:
							unique_projections_to_plane.append(holonomy)
					for i in range(len(MPr_sorted[P][copy_of_plane])):
						w0=MPr_sorted[P][copy_of_plane][i][0]
						for w1 in unique_projections_to_plane:
							S=matrix([[w0[0],w1[0]],[w0[1],w1[1]]])
							# M_inverse sends v0,v1 to w0,w1, respectively 
							M_inverse=S*T.inverse()							
							if det(M_inverse)==1:
								M=M_inverse.inverse()
								frobenius_norm_squared=M[0][0]^2+M[0][1]^2+M[1][0]^2+M[1][1]^2
								if frobenius_norm_squared<=norm_bound^2 and M not in matrices_to_check:
									matrices_to_check.append(M)

			# Here we construct groups of permutations of singularities and cylclic groups of copies of the plane associated to a particular singularity; these will be used below to check MP_2rho_sorted against various permutations of the bounded image of MPr_sorted under a matrix M in matrices_to_check

			# Records number of singularities of a particular order; will create symmetric groups on these numbers of elements and later take cartesian product of these groups for various permutations of singularities
			permutation_group_orders=[]
			angle_previous=None
			for i in range(len(cone_angles_sorted)):
				angle=cone_angles_sorted[i]
				if angle!=angle_previous:
					permutation_group_orders.append(cone_angles_sorted.count(angle))
				angle_previous=angle
			# Cartesian product of symmetric groups of orders from permutation_group_orders; will be used to permute singularities of the same order
			symm_groups=[]
			for i in permutation_group_orders:
				symm_groups.append(SymmetricGroup(i))
			permutations_of_singularities=cartesian_product(symm_groups)
			# Cartesian product of cyclic groups from orders of cone angles; will be used to check cyclic permutations of copies of the plane associated to a particular singularity
			cyc_groups=[]
			for j in cone_angles_sorted:
				cyc_groups.append(CyclicPermutationGroup(j))
			permutations_of_copies_of_plane=cartesian_product(cyc_groups)

			# A_r will keep the matrices with Frobenius norm bounded by b that are in the Veech group
			A_r=previous_Ar
			# We find some products of matrices already found to be in the Veech group
			for j in range(len(previous_Ar)):
				M=previous_Ar[j]
				for k in range(j,len(previous_Ar)):
					N=previous_Ar[k]
					L=M*N
					L_frobenius_norm_squared=L[0][0]^2+L[0][1]^2+L[1][0]^2+L[1][1]^2
					if L_frobenius_norm_squared<=norm_bound^2 and L not in A_r:
						A_r.append(L)
						if L.inverse() not in A_r:
							A_r.append(L.inverse())
					K=M*N
					K_frobenius_norm_squared=K[0][0]^2+K[0][1]^2+K[1][0]^2+K[1][1]^2
					if K_frobenius_norm_squared<=norm_bound^2 and K not in A_r:
						A_r.append(K)
						if K.inverse() not in A_r:
							A_r.append(K.inverse())


			# We will apply a matrix M from matrices_to_check to each point in MPr and keep only those bounded by $2*\rho$; we call this F_M_MPr_bounded_by_2rho
			for M in matrices_to_check:
				if M not in A_r and M not in matrices_checked:
					# We count how many elements are in F_M_MPr_bounded_by_2rho and compare with the number in MP_2rho as a quick test
					num_holonomies_in_F_M_MPr_bounded_by_2rho=0
					# The inverse of M is in A_r if and only if M is in A_r, so we need only check one of these matrices; at the end of this function, we append the inverse of M to A_r if M is found to be in A_r
					matrices_checked.append(M)
					if M.inverse() not in matrices_checked:
						matrices_checked.append(M.inverse())
					F_M_MPr_bounded_by_2rho=[]
					for P in range(num_singularities):
						number_of_planes=len(MPr[P])
						F_M_MPr_P_bounded_by_2rho=[[] for c in range(number_of_planes)]
						angle_of_image_of_e1=ang(M*vector((1,0)))
						for copy_of_plane in range(number_of_planes):
							for i in range(len(MPr[P][copy_of_plane])):
								mp=MPr[P][copy_of_plane][i][0]
								point=M*vector([mp[0],mp[1]])
								if point[0]^2+point[1]^2<=(2*rho)^2:
									num_holonomies_in_F_M_MPr_bounded_by_2rho+=1
									angle_of_point=ang(point)
									# If the counter-clockwise angle of `point' from angle_of_image_of_slice is greater than or equal to counterclockwise_angle_of_image_slice_from_original_slice, then we assign this marked period to this copy_of_plane; otherwise,it gets assigned to the next copy_of_plane (modulo the number of copies of the plane)
									# We store this data as a list: [holonomy, preimage's corresponding saddle connection on X (as recorded in flatsurf), cumulative angle of holonomy from (1,0) in a the 0th copy of the plane, cumulative angle of preimage of this marked period from (1,0) on 0th copy of plane]
									if angle_of_point>=angle_of_image_of_e1:
										F_M_MPr_P_bounded_by_2rho[copy_of_plane].append([point,MPr[P][copy_of_plane][i][1],copy_of_plane+angle_of_point,MPr[P][copy_of_plane][i][2]])
									else:
										F_M_MPr_P_bounded_by_2rho[(copy_of_plane+1)%number_of_planes].append([point,MPr[P][copy_of_plane][i][1],(copy_of_plane+1)%number_of_planes+angle_of_point,MPr[P][copy_of_plane][i][2]])

						# Sort holonomies in each copy of the plane by counterclockwise angle from original slice in chosen copy of plane
						for copy_of_plane in range(len(MPr[P])):
							F_M_MPr_P_bounded_by_2rho[copy_of_plane].sort(key=lambda x:x[2])
						F_M_MPr_bounded_by_2rho.append(F_M_MPr_P_bounded_by_2rho)
					# Sort the previous list by cone angles of singularities so we may easily check this list (and permutations of it) against MP_2rho_sorted
					F_M_MPr_bounded_by_2rho_sorted=sorted(F_M_MPr_bounded_by_2rho,key=len)

					# If the number of elements in F_M_MPr_bounded_by_2rho does not equal the number of elements in MP_2rho, we know M is not in A_r and need not check anything further
					if num_holonomies_in_F_M_MPr_bounded_by_2rho==num_holonomies_in_MP_2rho:
						# Now we test various permutations of singularities of the same cone angle within F_M_MPr_bounded_by_2rho_sorted as well as cyclic permutations of copies of the plane corresponding to each singularity
						# We first check if the holonomies of a permutation of F_M_MPr_bounded_by_2rho_sorted match those of MP_2rho_sorted; if so, we then test to see if the Z2-action is the same
						for permute_singularities in permutations_of_singularities.list():
							F_M_first_copy=deepcopy(F_M_MPr_bounded_by_2rho_sorted)
							place=0
							for n in range(len(permutation_group_orders)):
								num_sings_of_same_angle=permutation_group_orders[n]
								F_M_first_copy[place:place+num_sings_of_same_angle]=permute_singularities[n](F_M_first_copy[place:place+num_sings_of_same_angle])
								place+=num_sings_of_same_angle
							# For each permutation of singularities, we cyclically permute the copies of the plane associated to a particular singularity
							for permute_planes in permutations_of_copies_of_plane.list():
								F_M_copy=deepcopy(F_M_first_copy)
								for P in range(num_singularities):
									permute_planes_P=permute_planes[P]
									F_M_copy[P]=permute_planes_P(F_M_copy[P])
								# Now we check to see if the holonomies in each copy of the plane corresponding to each singularity of of MP_2rho_sorted match those of F_M_copy
								# List of holonomies associated to a particular singularity and copy of the plane
								F_M_copy_holonomies=[[[F_M_copy[P][copy_of_plane][i][0] for i in range(len(F_M_copy[P][copy_of_plane]))] for copy_of_plane in range(len(F_M_copy[P]))] for P in range(num_singularities)]
								# List of all saddle connection data from X; this is not associated to a particular singularity and copy of the plane, but corresponds to the same ordering of the information in F_M_copy_holonomies
								F_M_copy_saddle_connections=[F_M_copy[P][copy_of_plane][i][1] for P in range(num_singularities) for copy_of_plane in range(len(F_M_copy[P])) for i in range(len(F_M_copy[P][copy_of_plane]))]
								# If the holonomies are the same, then we check the Z2-action
								if F_M_copy_holonomies==MP_2rho_sorted_holonomies:
									Z2_action_agrees_flag=True
									for i in range(num_holonomies_in_MP_2rho):
										# We find the saddle connection on X corresponding to a particular holonomy in MP_2rho
										mp_saddle_connection=MP_2rho_sorted_saddle_connections[i]
										# For the same holonomy in F_M_copy (same even w/r/t singularity and copy of the plane since F_M_copy_holonomies=MP_2rho_sorted_holonomies and the orderings of F_M_copy_saddle_connections and MP_2rho_sorted_saddle_connections are the same), we find the saddle connection on X corresponding to the preimage under M of this holonomy vector
										Mmp_preimage_saddle_connection=F_M_copy_saddle_connections[i]
										# These are the identical but oppositely directed saddle connections of those found above
										mp_saddle_connection_inverse=mp_saddle_connection.invert()
										Mmp_preimage_saddle_connection_inverse=Mmp_preimage_saddle_connection.invert()
										# For the Z2-action to match, these inverse saddle connections must correspond to the same holonomies in F_M_copy_holonomies and MP_2rho_sorted_holonomies; otherwise the Z2-action is different and F_M_copy does not equal MP_2rho
										if MP_2rho_sorted_saddle_connections.index(mp_saddle_connection_inverse)!=F_M_copy_saddle_connections.index(Mmp_preimage_saddle_connection_inverse):
											Z2_action_agrees_flag=False
											break
									if Z2_action_agrees_flag==True:
										A_r.append(M)
										# M.inverse() is in A_r if and only if M is in A_r
										if M.inverse() not in A_r:
											A_r.append(M.inverse())
			return([A_r,matrices_to_check])

		[Ar,matrices_checked_r]=A_r(r,b)

		#################################################################################################################################
		#################################################################################################################################










		#######################################################################################################################################################################################################################################################
		#######################################################################################################################################################################################################################################################
		## 7.2.5: If $M\in A_r$ with $M\notin\{\pm\text{Id}_{\text{S0}(2,\mathbb{R})}\}$, then find an $M_0\in\text{SL}(2,\mathbb{R})$ such that $\Gamma(M_0\cdot (X,\omega))\cap\text{SO}(2,\mathbb{R})\subseteq\{\pm\text{Id}_{\text{SO}(2,\mathbb{R})}\}$ ##
		## 7.2.6: Let $(X_0,\omega_0)=M_0\cdot (X,\omega)$																																																	 ##
		#######################################################################################################################################################################################################################################################
		#######################################################################################################################################################################################################################################################

		if len(Ar)>2 or (len(Ar)==2 and matrix([[-1,0],[0,-1]]) not in Ar):
			n+=1
		else:
			trivial_stabilizer=True

		#######################################################################################################################################################################################################################################################
		#######################################################################################################################################################################################################################################################










	#############################################################################################################
	#############################################################################################################
	## 7.2.7: Execute Algorithm 7.1 with input: $(X_0,\omega_0)$ and rename the output: $AssociatedGenerators$ ##
	#############################################################################################################
	#############################################################################################################

	# We have exited the above while loop, i.e. have found a matrix $M0$ such that the Veech group of $M0\cdot X$ has trivial stabilizer 
	# Here we begin to run Algorithm 7.1 of Edwards, beginning from step 7.1.5 as steps 7.1.1--7.1.4 were done above in Algorithm 7.2

	#############################################################################################################
	#############################################################################################################










	##########################################################################
	##########################################################################
	## 7.1.5: If $-\text{Id}\in A_r$, then let $ContainsMinusIdentity=TRUE$ ##
	## 7.1.6: else let $ContainsMinusIdentity=FALSE$                        ##
	##########################################################################
	##########################################################################

	MinusIdentity=matrix([[-1,0],[0,-1]])

	if MinusIdentity in Ar:
		ContainsMinusIdentity=True
	else:
		ContainsMinusIdentity=False 

	##########################################################################
	##########################################################################










	###########################################
	###########################################
	## 7.1.7: Let $ContainmentVolume=\infty$ ##
	###########################################
	###########################################

	ContainmentVolume=infinity

	###########################################
	###########################################










	########################################################################################################
	########################################################################################################
	## 7.1.8: Do while $ContainmentVolume=\infty$:                                                        ##
	##  (a) Double the value of $r$ and let $b=\chi_1^{-1}(2\rho/r)$                                      ##
	##  (b) Calculate $MP_P^r(X,\omega)$ for all $P\in\Sigma$                                             ##
	##  (c) Use Theorem 18 to complete the set $A_{r/2}$ to $A_r=\{M\in\Gamma(X,\omega)\ |\ ||M||\le b\}$ ##
	##  (d) Construct $\Omega(\overline{A}_r)$                                                            ##
	##  (e) Let $ContainmentVolume=\nu_{\mathbb{H}}(\Omega(\overline{A}_r))$                              ##
	########################################################################################################
	########################################################################################################

	iteration=0
	while ContainmentVolume==infinity and iteration<iteration_limit:
		r=2*r
		b=chi_1_inv(2*rho/r)
		[Ar,matrices_checked_r]=A_r(r,b,previous_Ar=Ar,matrices_checked=matrices_checked_r)
		Ar_mod_MinusIdentity=[matrix([[1,0],[0,1]])]
		if ContainsMinusIdentity==True:
			for M in Ar:
				if M not in Ar_mod_MinusIdentity and -M not in Ar_mod_MinusIdentity:
					Ar_mod_MinusIdentity.append(M)
		else:
			Ar_mod_MinusIdentity=Ar

		# Below we determine whether or not the hyperbolic polygon determined by the group generated by Ar_mod_MinusIdentity has finite hyperbolic area by considering whether or not it has free sides
		# We defer the actual construction of this polygon and the calculation of its area until the loop in step 13 where it is necessary

		# The perpendicular bisectors corresponding to each M*I for M in Ar_mod_MinusIdentity
		bisectors=[geodesic(matrix=M) for M in Ar_mod_MinusIdentity if M!=matrix([[1,0],[0,1]])]
		# # Picture of the perpendicular bisectors found thus far
		# sum([g.plot() for g in bisectors])

		# When the cardinality of free_sides is finite (from the construction of difference(), this is when free_sides==[]), we'll know ContainmentVolume is finite
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
				ContainmentVolume='Is Finite'
				break

		iteration+=1

	###########################################################
	# HERE WE END THE FUNCTION IF ContainmentVolume==infinity #
	###########################################################
	if ContainmentVolume==infinity:
		print('Iteration limit reached.  An exhaustive list of elements of the Veech group with Frobenius norm bounded above by {} is given below:'.format(b))
		return [M0.inverse()*M*M0 for M in Ar]

	########################################################################################################
	########################################################################################################










	###########################################
	###########################################
	## 7.1.9: Let $BoundingRadius=\chi_2(b)$ ##
	###########################################
	###########################################

	BoundingRadius=chi_2(b)

	###########################################
	###########################################









	###############################################################################
	###############################################################################
	## 7.1.10: Let $BoundedPiece=\Omega(\overline{A_r})\cap B(i,BoundingRadius)$ ##
	###############################################################################
	###############################################################################

	# Since steps 11 and 12 declare $AllSidesRepresented=False$ and $ParabolicCycles=False$, we don't actually need to find BoundedPiece yet; the loop in step 13 will be entered regardless
	# To enter the while loop in step 13 below, we artificially let 0 be a lower bound on the volume of BoundedPiece and 100 the volume of the Dirichlet region corresponding to Ar, and we defer these computations until said loop
	BoundedPiece_Volume_lower_bound=0
	Omega_Ar_Volume=100

	###############################################################################
	###############################################################################









	#############################################
	#############################################
	## 7.1.11: Let $AllSidesRepresented=False$ ##
	#############################################
	#############################################

	AllSidesRepresented=False

	#############################################
	#############################################









	#########################################
	#########################################
	## 7.1.12: Let $ParabolicCycles=False$ ##
	#########################################
	#########################################

	# The test to make ParabolicCycles=True in step 7.1.13 can be shown to be extraneous

	#########################################
	#########################################










	###############################################################################################################################################################################################################################
	###############################################################################################################################################################################################################################
	## 7.1.13: Do while $AllSidesRepresented=False$ or $ParabolicCycles=False$ or $\nu_{\mathbb{H}}(BoundedPiece)\le\nu_{\mathbb{H}}(\Omega(\overline{A}_r))$:																	 ##
	##  (a) Double the value of $r$ and let $b=\chi_1^{-1}(2\rho/r)$																																							 ##
	##  (b) Recalculate $MP_P^r(X,\omega)$ for each $P\in\Sigma$, compute $A_r$ and $\Omega(\overline{A_r})$																													 ##
	##  (c) Calculate $BoundingRadius=\chi_2(b)$																																												 ##
	##  (d) Let $BoundedPiece=\Omega(\overline{A_r})\cap B(i,BoundingRadius)$																																					 ##
	##  (e) If the only sides of $\Omega(\overline{A_r})$ not contained in $BoundedPiece$ have one endpoint on the line at infinity and the other an interior point of $B(i,BoundedRadius)$, then let $AllSidesRepresented=True$ ##
	##    (i)   Calculate the ideal vertex cycles associated to the side pairing transformations on the sides containing an ideal vertex endpoint																				 ##
	##    (ii)  If all ideal vertex cycles are parabolic, then let $ParabolicCycles=True$																																		 ##
	##    (iii) else let $ParabolicCycles=False$																																												 ##
	##  (f) else let $AllSidesRepresented=False$ and let $ParabolicCycles=False$																																			     ##
	###############################################################################################################################################################################################################################
	###############################################################################################################################################################################################################################

	# The condition that ParabolicCycles be true is extraneous; it can be shown that it is sufficient to have AllSidesRepresented and the volume of the bounded piece greater than half the volume of the polygon generated by Omega_Ar
	# For the first run of the while loop, we restrain from doubling r as it may not be necessary; fewer computations of Ar speeds up the algoritm considerably
	first_run=True
	while AllSidesRepresented==False or BoundedPiece_Volume_lower_bound.n(prec)<=(Omega_Ar_Volume.n(prec)/2):
		Ar_previous=copy(Ar)
		if first_run==False:
			r=2*r
			b=chi_1_inv(2*rho/r)
			[Ar,matrices_checked_r]=A_r(r,b,previous_Ar=Ar,matrices_checked=matrices_checked_r)
		first_run=False

		# To be used to construct the new hyperbolic polygon determined by Ar, which we denote Omega_Ar
		Ar_mod_MinusIdentity=[matrix([[1,0],[0,1]])]
		if ContainsMinusIdentity==True:
			for M in Ar:
				if M not in Ar_mod_MinusIdentity and -M not in Ar_mod_MinusIdentity:
					Ar_mod_MinusIdentity.append(M)
		else:
			Ar_mod_MinusIdentity=Ar

		# We need only append geodesics corresponding to newly found matrices to the list bisectors
		for M in Ar_mod_MinusIdentity:
			if M not in Ar_previous:
				bisectors.append(geodesic(matrix=M))
		# Picture of the perpendicular bisectors found thus far
		# sum([g.plot() for g in bisectors])

		# The hyperbolic polygon determined by Ar
		Omega_Ar=Omega(bisectors,prec=prec,timeout=timeout)
		# Picture of Omega_Ar
		# sum([g.plot() for g in Omega_Ar])

		BoundingRadius=chi_2(b)

		# A computation with the integral definition of hyperbolic distance gives the following imaginary part of the Euclidean center (the center lies on the imaginary axis) and Euclidean radius of the hyperbolic ball centered at I with radius BoundingRadius
		BoundingBall_euc_center=cosh(BoundingRadius).n(prec)
		BoundingBall_euc_radius=sinh(BoundingRadius).n(prec)
		# Picture of Omega_Ar and bounding ball
		# sum([g.plot() for g in Omega_Ar]+[circle((0,BoundingBall_euc_center),BoundingBall_euc_radius,edgecolor='red')])

		# We will fill this list with the geodesics contributing to Omega_Ar and some additional geodesics that depend on the bounding ball; see comments after the following for loop
		bounded_polygon_geodesics=[]

		# We determine whether or not all sides of the polygon are represented, as described in step 13(e) above
		# If any side is not represented, we'll make this flag false and break the loop
		represented_flag=True
		# We keep track of which sides of the polygon intersect the ball, and where the intersections on the boundary of the ball occur; each entry will be a list [[x,y],gamma,vertex_orientation], where x,y are the real and imaginary parts of this intersection; gamma is the side in question; and vertex_orientation='counterclockwise' if the segment of gamma that the bounding ball 'cuts off' countains the counterclockwise vertex of gamma, and vertex_orientation='clockwise' otherwise
		ball_intersections=[]
		for gamma in Omega_Ar:
			bounded_polygon_geodesics.append(gamma)
			# Gamma is vertical
			if gamma.type=='Vertical':
				# Check if the h-line that gamma lies on intersects the boundary of the bounding ball; note that the inequalities here are strict so that gamma may not be tangent to the ball but must actually intersects it
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
					left_interesction_imgy_pt=sqrt(gamma.radius^2-(left_intersection-gamma.center)^2)
					right_interesction_imgy_pt=sqrt(gamma.radius^2-(right_intersection-gamma.center)^2)
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
							ball_intersections.append([[left_intersection,left_interesction_imgy_pt],gamma,left_vertex_orientation])
						if right_vertex>right_intersection:
							ball_intersections.append([[right_intersection,right_interesction_imgy_pt],gamma,right_vertex_orientation])
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
		# This gives a lower bound on the hyperbolic area of Omega_Ar intersected with the ball

		# Note that we only need to compute this if BoundedPiece_Volume_lower_bound<=Omega_Ar_Volume/2; once BoundedPiece_Volume_lower_bound>Omega_Ar_Volume/2 for some r, then this is true for all larger r' since the bounding ball is growing and the Dirichlet polygon Omega_Ar is shrinking as r increases

		if AllSidesRepresented==True and BoundedPiece_Volume_lower_bound.n(prec)<=Omega_Ar_Volume.n(prec)/2:
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
							if (intersection0[-1]=='counterclockwise' and ((gamma1.vertex_clockwise[0].n(prec)==vertex[0].n(prec) and gamma1.vertex_clockwise[1].n(prec)==vertex[1].n(prec)) or (gamma1.vertex_clockwise[1]==infinity and vertex[1]==infinity))) or (intersection0[-1]=='clockwise' and ((gamma1.vertex_counterclockwise[0].n(prec)==vertex[0].n(prec) and gamma1.vertex_counterclockwise[1].n(prec)==vertex[1].n(prec)) or (gamma1.vertex_counterclockwise[1]==infinity and vertex[1]==infinity))):
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
			Omega_Ar_bounding_ball=Omega(bounded_polygon_geodesics,prec=prec,timeout=timeout)
			# Picture of Omega_Ar_bounding_ball and bounding ball
			# sum([g.plot() for g in Omega_Ar_bounding_ball]+[circle((0,BoundingBall_euc_center),BoundingBall_euc_radius,edgecolor='red')])

			BoundedPiece_Volume_lower_bound=PolygonVolume(Omega_Ar_bounding_ball,prec=prec)
			Omega_Ar_Volume=PolygonVolume(Omega_Ar,prec=prec)

	###############################################################################################################################################################################################################################
	###############################################################################################################################################################################################################################










	#########################################################################################################################################################
	#########################################################################################################################################################
	## 7.1.14: Let $SidePairingRepresentatives$ be the set of elements in $A_r$ associated to the side pairing transformations of $\Omega(\overline{A}_r)$ ##
	#########################################################################################################################################################
	#########################################################################################################################################################

	# List of all matrices corresponding to sides of Omega_Ar
	all_matrices=[gamma.matrix for gamma in Omega_Ar]
	# We don't need to include a matrix if its inverse is already include, so we get rid of these non-essential matrices
	SidePairingRepresentatives=[]
	for M in all_matrices:
		if M.inverse() not in SidePairingRepresentatives and -M.inverse() not in SidePairingRepresentatives:
			SidePairingRepresentatives.append(M)

	#########################################################################################################################################################
	#########################################################################################################################################################









	###################################################################################################################################################################
	###################################################################################################################################################################
	## 7.1.15: If $ContainsMinusIdentity=TRUE$, then let $Generators=SidePairingRepresentatives\cup\{-\text{Id}\}$; else let $Generators=SidePairingRepresentatives$ ##
	###################################################################################################################################################################
	###################################################################################################################################################################

	# Step 7.2.7 renames the output of Alorithm 7.1 `AssociatedGenerators;' so we do so here
	AssociatedGenerators=SidePairingRepresentatives
	if ContainsMinusIdentity==True:
		AssociatedGenerators.append(-matrix.identity(2))

	# This ends Algorithm 7.1, and we return to Algorithm 7.2 below

	###################################################################################################################################################################
	###################################################################################################################################################################








	#########################################################################
	#########################################################################
	## 7.2.8: Let $Generators=M_0^{-1}\cdot AssociatedGenerators\cdot M_0$ ##
	#########################################################################
	#########################################################################


	Generators=[M0.inverse()*M*M0 for M in AssociatedGenerators]

	#########################################################################
	#########################################################################









	################################
	################################
	## 7.2.9: Output $Generators$ ##
	################################
	################################

	print('The Veech group is a lattice generated by the following matrices:')
	return Generators

	################################
	################################