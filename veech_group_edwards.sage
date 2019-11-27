r''' 
ALGORITHM 7.1
Takes as input translation surface (X,omega) such that $\Gamma(X,omega)\cap \text{SO}(2,\mathbb{R})\subseteq \{\pm \text{Id}\}$

X must be a translation surface from the flatsurf package
Example: 
sage: from flatsurf import *
sage: X=translation_surfaces.veech_2n_gon(5)
'''

# # Surface built from 2 regular pentagons and 2 squares
# p0=polygons.regular_ngon(5)
# p1=matrix([[-1,0],[0,-1]])*p0
# p2=polygons.regular_ngon(4)
# p3=polygons.regular_ngon(4)
# S=Surface_list(p0.base_ring())
# S.add_polygon(p0)
# S.add_polygon(p1)
# S.add_polygon(p2)
# S.add_polygon(p3)
# S.change_edge_gluing(0,1,1,1)
# S.change_edge_gluing(0,2,1,2)
# S.change_edge_gluing(0,3,1,3)
# S.change_edge_gluing(0,4,1,4)
# S.change_edge_gluing(0,0,2,2)
# S.change_edge_gluing(1,0,3,0)
# S.change_edge_gluing(2,0,3,2)
# S.change_edge_gluing(2,1,2,3)
# S.change_edge_gluing(3,1,3,3)
# X=TranslationSurface(S)


# Will eventually uncomment this function when all of the subsections of the algorithm work
# def veech_group(X):
from flatsurf.geometry.polygon import polygons,PolygonPosition
assert isinstance(X,TranslationSurface)










###############################################
###############################################
## 0: Compute the Voronoi decomposition of X ##
###############################################
###############################################

# Finds vertex equivalence classes of X
X._num_polygons=X.num_polygons()
X._num_singularities=X.num_singularities()
vertex_equivalence_classes=[]
vertex_equivalence_classes_sets=[]
while len(vertex_equivalence_classes)<X._num_singularities:
	for polygon in range(X._num_polygons):
		for vertex in range(X.polygon(polygon).num_edges()):
			vertex_is_equivalent_to=X.singularity(polygon,vertex).vertex_set()
			already_accounted=0
			for i in vertex_equivalence_classes_sets:
				if i==vertex_is_equivalent_to:
					already_accounted=1
					break
			if already_accounted==0:
				vertex_equivalence_classes_sets.append(vertex_is_equivalent_to)
				vertex_equivalence_classes.append([v for v in vertex_is_equivalent_to])

# Gives the vertices of each polygon of X in AA (algebraic reals) so that the Voronoi decomposition can be computed with respect to these points
vertices=[[[AA(c) for c in vert] for vert in X.polygon(poly).vertices()] for poly in range(X._num_polygons)]
# Gives a list whose entries are the Voronoi decompositions with respect to the vertices of each polygon
VD=[VoronoiDiagram(vertices[poly]) for poly in range(X._num_polygons)]
# A list of lists; each entry of the outer list corresponds to a polygon of X; each entry of the inner lists is the (potentially) unbounded polyhedron in the Voronoi decomposition associated to a particular vertex of said polygon
VD_polyhedra_unbdd=[[VD[poly].regions()[VD[poly].points()[v]] for v in range(X.polygon(poly).num_edges())] for poly in range(X._num_polygons)]
# Same structure as the list above, but now the polyhedra are intersected with the proper polygon of X, so they are bounded
VD_polyhedra=[[VD_polyhedra_unbdd[poly][vert].intersection(Polyhedron(vertices=vertices[poly])) for vert in range(X.polygon(poly).num_edges())] for poly in range(X._num_polygons)]

###############################################
###############################################










##################################################################
##################################################################
## 1: Calculate $\rho=\rho(X,\omega)$ using the Voronoi 2-cells ##
##################################################################
##################################################################

# Creates a list of the max distances in each polygon of the Voronoi decomposition from the singularity of X (main_vertex below) to the boundary of the polygon
distances_from_singularities=[]
for poly in range(X._num_polygons):
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










################################################################
################################################################
## 2: Let $r=2\rho$ and let $b=\chi_1^{-1}(2\rho/r)=\sqrt{2}$ ##
################################################################
################################################################

def chi_1_inv(x):
	return sqrt(x^2+x^(-2))

r=2*rho
b=sqrt(2)

################################################################
################################################################










##########################################################
##########################################################
## 3: Calculate $MP_P^r(X,\omega)$ for all $P\in\Sigma$ ##
##########################################################
##########################################################

# This records information to determine which disjoint copy of the plane an element of $MP_P^r(X,\omega)$ should lie in
# vertices_in_disjoint_planes[P][n] gives a list of vertex data (polygon,vertex) such that separatrices leaving these vertices should develop in the n^th copy of the plane corresponding to singularity P
# There is ambiguity for some vertices, such as vertex 5 of the regular hexagon...these ambiguous vertices where separatrices may develop in either the n^th or (n+1)^st plane lead the list corresponding to the (n+1)^st plane
vertices_in_disjoint_planes=[]
first_edge_from_singularity=[]
cone_angles=[]
for P in range(X._num_singularities):
	vertex_angles=[X.polygon(vertex_data[0]).angle(vertex_data[1]) for vertex_data in vertex_equivalence_classes[P]]
	cone_angles.append(sum(vertex_angles))
	copy_of_plane=0
	
	# Will contain lists where the (copy_of_plane)^th list contains vertex data (polygon,vertex) with vertices corresponding to the (copy_of_plane)^th copy of the plane that corresponds to singularity P
	vertices_in_disjoint_planes_P=[]

	# Canonicalizes the set of vertices corresponding to singularity P by vertex data (polygon,vertex): sorts by lowest polygon index then lowest vertex index
	singularity_P=sorted(vertex_equivalence_classes[P])

	vertex_data=singularity_P[0]
	first_edge_from_singularity.append(X.polygon(vertex_data[0]).edge(vertex_data[1]))
	P_angle=X.polygon(vertex_data[0]).angle(vertex_data[1])

	# Will be the (copy_of_plane)^th list inside vertices_in_disjoint_planes
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

# To be used below so that holonomies may be sorted counter-clockwise
def ang(holonomy):
	x=holonomy[0]
	y=holonomy[1]
	if x>0:
		if y>=0:
			theta=arctan(y/x)
		else:
			theta=arctan(y/x)+2*pi
	elif x<0:
			theta=arctan(y/x)+pi
	else:
		if y>0:
			theta=pi/2
		else: theta=3*pi/2
	return theta

### To do: create an optional parameter to add to previously computed marked_periods(radius0) ###
def marked_periods(radius):
	MPr=[]
	for P in range(X._num_singularities):
		MPr_P=[[] for i in range(cone_angles[P])]
		for vertex_data in vertex_equivalence_classes[P]:
			copy_of_plane=[n for n in range(len(vertices_in_disjoint_planes[P])) if vertex_data in vertices_in_disjoint_planes[P][n]][0]
			polygon_start=vertex_data[0]
			vertex_start=vertex_data[1]
			saddle_connections=X.saddle_connections(radius^2,polygon_start,vertex_start)
			saddle_connections_directions=[sc.direction() for sc in saddle_connections]
			# First see if the current vertex is ambiguous (in the sense that saddle connections from it may be assigned to the current copy of the plane or the previous copy)
			if vertex_data==vertices_in_disjoint_planes[P][copy_of_plane][0]:
				for i in range(len(saddle_connections)):
					sc=saddle_connections[i]
					v0=first_edge_from_singularity[P]
					v1=saddle_connections_directions[i]
					M=matrix([[v0[0],v1[0]],[v0[1],v1[1]]])
					holonomy=sc.holonomy()
					# If sc leaves its singularity clockwise of where the current copy of the plane begins, then its holonomy vector is assigned to the previous copy of the plane
					if det(M)<0:
						MPr_P[copy_of_plane-1].append([holonomy,sc,(copy_of_plane-1)*2*pi+ang(holonomy)])
					# Otherwise its holonomy vector is assigned to the current copy of the plane
					else:
						MPr_P[copy_of_plane].append([holonomy,sc,copy_of_plane*2*pi+ang(holonomy)])
			# These are the unambiguous vertices
			else:	
				for i in range(len(saddle_connections)):
					sc=saddle_connections[i]
					holonomy=sc.holonomy()
					MPr_P[copy_of_plane].append([holonomy,sc,copy_of_plane*2*pi+ang(holonomy)])
		MPr.append(MPr_P)
		# Sort holonomies in each copy of the plane by cumulative counterclockwise angle from (1,0) on 0th copy of plane
		for copy_of_plane in range(len(MPr[P])):
			MPr[P][copy_of_plane].sort(key=lambda x:x[2])
	return MPr

MP_2rho=marked_periods(2*rho)
# Sort MP_2rho by cone angles of singularities so we may easily check this list against (permutations of) F_M_MPr_bounded_by_2rho below
MP_2rho_sorted=sorted(MP_2rho,key=len)
cone_angles_sorted=sorted(cone_angles)

##########################################################
##########################################################










#############################################################################################################################
#############################################################################################################################
## 4: Calculate $A_r=\{M\in\Gamma(X,\omega)\ |\ ||M||\le b\}=\text{SO}(2,\mathbb{R})\cap\Gamma(X,\omega)$ using Theorem 18 ##
#############################################################################################################################
#############################################################################################################################

### To do: create an optional parameter to add to previously computed A_r(radius0,norm_bound0) ###
def A_r(radius,norm_bound):	
	# matrices_to_check will be a list comprised of matrices in $\text{SL}_2\mathbb{R}$ whose Frobenius norms are bounded above by norm_bound and whose inverses take a basis of vectors in marked_periods(2*rho) to another basis of vectors in marked_periods(radius)
	matrices_to_check=[]
	MPr=marked_periods(radius)
	MPr_sorted=sorted(MPr,key=len)

	# Basis from projection of MP_2rho onto complex plane
	v0=MP_2rho_sorted[0][0][0][0]
	v1=MP_2rho_sorted[0][0][1][0]
	T=matrix([[v0[0],v1[0]],[v0[1],v1[1]]])
	# Note v0,v1 correspond to a singularity of minimal angle; below we will find another pair of marked periods that can be sent to v0,v1; these must correspond to a singularity of minimal angle as well
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

	# This will keep the matrices with Frobenius norm bounded by b that are in the Veech group
	A_r=[]

	# We will apply a matrix M from matrices_to_check to each point in MPr and keep only those bounded by $2*\rho$; we call this F_M_MPr_bounded_by_2rho
	for M in matrices_to_check:
		F_M_MPr_bounded_by_2rho=[]
		for P in range(X._num_singularities):
			number_of_planes=len(MPr[P])
			F_M_MPr_P_bounded_by_2rho=[[] for c in range(number_of_planes)]
			angle_of_image_of_first_vector=ang(M*MPr[P][0][0][0])
			for copy_of_plane in range(number_of_planes):
				for i in range(len(MPr[P][copy_of_plane])):
					mp=MPr[P][copy_of_plane][i][0]
					point=M*vector([mp[0],mp[1]])
					if point[0]^2+point[1]^2<=(2*rho)^2:
						angle=ang(point)
						# If the angle of point under is between angle_of_image_of_first_vector and 2*pi, we assign this marked period to copy_of_plane; otherwise, it gets assigned to the next copy_of_the_plane (mod the number of copies of the plane)
						# We store this data as a list: [holonomy, preimage's corresponding saddle connection on X (as recorded in flatsurf), angle of holonomy from (1,0) in a single copy of the plane, cumulative angle of preimage of this marked period from (1,0) on 0th copy of plane]
						if angle>=angle_of_image_of_first_vector:
							F_M_MPr_P_bounded_by_2rho[copy_of_plane].append([point,MPr[P][copy_of_plane][i][1],angle,MPr[P][copy_of_plane][i][2]])
						else:
							F_M_MPr_P_bounded_by_2rho[(copy_of_plane+1)%number_of_planes].append([point,MPr[P][copy_of_plane][i][1],angle,MPr[P][copy_of_plane][i][2]])
			# Sort holonomies in each copy of the plane by counterclockwise angle from (1,0) in chosen copy of plane
			for copy_of_plane in range(len(MPr[P])):
				F_M_MPr_P_bounded_by_2rho[copy_of_plane].sort(key=lambda x:x[2])
			F_M_MPr_bounded_by_2rho.append(F_M_MPr_P_bounded_by_2rho)
			# Sort the previous list by cone angles of singularities so we may easily check this list (and permutations of it) against MP_2rho_sorted
			F_M_MPr_bounded_by_2rho_sorted=sorted(F_M_MPr_bounded_by_2rho,key=len)

		# Now we test various permutations of singularities of the same cone angle within F_M_MPr_bounded_by_2rho_sorted as well as cyclic permutations of copies of the plane corresponding to each singularity
		# We first check if the holonomies of a permutation of F_M_MPr_bounded_by_2rho_sorted match those of MP_2rho_sorted; if so, we then test to see if the $\mathbb{Z}_2$-action is the same
		for permute_singularities in permutations_of_singularities.list():
					F_M_first_copy=deepcopy(F_M_MPr_bounded_by_2rho_sorted)
					place=0
					for n in range(len(permutation_group_orders)):
						num_sings_of_same_angle=permutation_group_orders[n]
						F_M_first_copy[place:place+num_sings_of_same_angle]=permute_singularities[n](F_M_first_copy[place:place+num_sings_of_same_angle])
						place+=num_sings_of_same_angle
					for permute_planes in permutations_of_copies_of_plane.list():
						F_M_copy=deepcopy(F_M_first_copy)
						for P in range(X._num_singularities):
							permute_planes_P=permute_planes[P]
							F_M_copy[P]=permute_planes_P(F_M_copy[P])
						M_is_a_keeper=True
						for P0 in range(X._num_singularities):
							for copy_of_plane0 in range(len(F_M_copy[P0])):
								for i in range(len(F_M_copy[P0][copy_of_plane0])):
									# Check that the holonomies associated to a particular singularity and copy of the plane match; if so, check the $\mathbb{Z}_2$-action
									if F_M_copy[P0][copy_of_plane0][i][0]==MP_2rho_sorted[P0][copy_of_plane0][i][0]:
										F_M_saddle_connection_on_X=F_M_copy[P0][copy_of_plane0][i][1]
										F_M_saddle_connection_preimage_angle=F_M_copy[P0][copy_of_plane0][i][3]
										MPr_saddle_connection_on_X=MP_2rho_sorted[P0][copy_of_plane0][i][1]
										MPr_saddle_connection_angle=MP_2rho_sorted[P0][copy_of_plane0][i][2]
										F_M_saddle_connection_on_X_opposite=F_M_saddle_connection_on_X.invert()
										MPr_saddle_connection_on_X_opposite=MPr_saddle_connection_on_X.invert()
										F_M_found_opposite=False
										MPr_found_opposite=False
										for P1 in range(X._num_singularities):
											for copy_of_plane1 in range(len(MPr[P1])):
												for j in range(len(MPr[P1][copy_of_plane1])):
													if MPr[P1][copy_of_plane1][j][1]==F_M_saddle_connection_on_X_opposite:
														### THINK ABOUT THIS IF X HAS MULTIPLE SINGULARITIES OF DIFFERENT CONE ANGLES ###
														F_M_angle_difference=(int((MPr[P1][copy_of_plane1][j][2]-F_M_saddle_connection_preimage_angle)/pi)%(2*len(MPr[P1])))*pi
														F_M_Z2_index=[P1,F_M_angle_difference]
														F_M_found_opposite=True
													if MPr[P1][copy_of_plane1][j][1]==MPr_saddle_connection_on_X_opposite:
														MPr_angle_difference=(int((MPr[P1][copy_of_plane1][j][2]-MPr_saddle_connection_angle)/pi)%(2*len(MPr[P1])))*pi
														MPr_Z2_index=[P1,MPr_angle_difference]
														MPr_found_opposite=True
													if F_M_found_opposite==True and MPr_found_opposite==True:
														if F_M_Z2_index!=MPr_Z2_index:
															M_is_a_keeper=False
														break
												if M_is_a_keeper==False:
													break
											if M_is_a_keeper==False:
												break
									else:
										M_is_a_keeper=False
										break
								if M_is_a_keeper==False:
									break
							if M_is_a_keeper==False:
								break
						if M_is_a_keeper==True:
							A_r.append(M)
							break
					if M_is_a_keeper==True:
						break	
	return(A_r)

Ar=A_r(r,b)

#############################################################################################################################
#############################################################################################################################










######################################################################
######################################################################
## 5: If $-\text{Id}\in A_r$, then let $ContainsMinusIdentity=TRUE$ ##
## 6: else let $ContainsMinusIdentity=FALSE$                        ##
######################################################################
######################################################################

MinusIdentity=matrix([[-1,0],[0,-1]])

if MinusIdentity in Ar:
	ContainsMinusIdentity=True
else:
	ContainsMinusIdentity=False 

######################################################################
######################################################################










#######################################
#######################################
## 7: Let $ContainmentVolume=\infty$ ##
#######################################
#######################################

ContainmentVolume=infinity

#######################################
#######################################










########################################################################################################
########################################################################################################
## 8: Do while $ContainmentVolume=\infty$:                                                            ##
##  (a) Double the value of $r$ and let $b=\chi_1^{-1}(2\rho/r)$                                      ##
##  (b) Calculate $MP_P^r(X,\omega)$ for all $P\in\Sigma$                                             ##
##  (c) Use Theorem 18 to complete the set $A_{r/2}$ to $A_r=\{M\in\Gamma(X,\omega)\ |\ ||M||\le b\}$ ##
##  (d) Construct $\Omega(\overline{A}_r)$                                                            ##
##  (e) Let $ContainmentVolume=\nu_{\mathbb{H}}(\Omega(\overline{A}_r))$                              ##
########################################################################################################
########################################################################################################

while ContainmentVolume=infinity:
	r=2*r
	b=chi_1_inv(2*rho/r)
	Ar=A_r(r,b)
	Ar_mod_MinusIdentity=[matrix([[1,0],[0,1]])]
	if ContainsMinusIdentity==True:
		for M in Ar:
			if M not in Ar_mod_MinusIdentity and -M not in Ar_mod_MinusIdentity:
				Ar_mod_MinusIdentity.append(M)
	else:
		Ar_mod_MinusIdentity=Ar

########################################################################################################
########################################################################################################