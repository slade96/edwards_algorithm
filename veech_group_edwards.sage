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

iteration_limit=3
# Will eventually uncomment this function when all of the subsections of the algorithm work
# def veech_group(X,iteration_limit=3):
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
# There is ambiguity for some vertices, such as vertex 6 of the regular octagon...these ambiguous vertices where separatrices may develop in either the n^th or (n+1)^st plane lead the list corresponding to the (n+1)^st plane
vertices_in_disjoint_planes=[]
first_edge_from_singularity=[]
cone_angles=[]
for P in range(X._num_singularities):
	vertex_angles=[X.polygon(vertex_data[0]).angle(vertex_data[1]).numerical_approx() for vertex_data in vertex_equivalence_classes[P]]
	cone_angles.append(int(sum(vertex_angles)))
	copy_of_plane=0
	
	# Will contain lists where the (copy_of_plane)^th list contains vertex data (polygon,vertex) with vertices corresponding to the (copy_of_plane)^th copy of the plane that corresponds to singularity P
	vertices_in_disjoint_planes_P=[]

	# Canonicalizes the set of vertices corresponding to singularity P by vertex data (polygon,vertex): sorts by lowest polygon index then lowest vertex index
	singularity_P=sorted(vertex_equivalence_classes[P])

	vertex_data=singularity_P[0]
	first_edge_from_singularity.append(X.polygon(vertex_data[0]).edge(vertex_data[1]))
	P_angle=X.polygon(vertex_data[0]).angle(vertex_data[1])

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
		angle_of_slice_of_planes=ang(first_edge_from_singularity[P])
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
					holonomy_angle=ang(holonomy)
					# If sc leaves its singularity clockwise of where the current copy of the plane is sliced, then its holonomy vector is assigned to the previous copy of the plane
					if det(M)<0:
						# MPr_P[copy_of_plane-1].append([holonomy,sc,(copy_of_plane-1)*2*pi+ang(holonomy)])
						if holonomy_angle>angle_of_slice_of_planes:
							cumulative_angle_from_initial_slice=2*pi*(copy_of_plane-1)+(holonomy_angle-angle_of_slice_of_planes)
							MPr_P[copy_of_plane-1].append([holonomy,sc,cumulative_angle_from_initial_slice])
						else:
							cumulative_angle_from_initial_slice=2*pi*(copy_of_plane)+(holonomy_angle-angle_of_slice_of_planes)
							MPr_P[copy_of_plane-1].append([holonomy,sc,cumulative_angle_from_initial_slice])
					# Otherwise its holonomy vector is assigned to the current copy of the plane
					else:
						if holonomy_angle>=angle_of_slice_of_planes:
							cumulative_angle_from_initial_slice=2*pi*(copy_of_plane)+(holonomy_angle-angle_of_slice_of_planes)
							MPr_P[copy_of_plane].append([holonomy,sc,cumulative_angle_from_initial_slice])
						else:
							cumulative_angle_from_initial_slice=2*pi*(copy_of_plane+1)+(holonomy_angle-angle_of_slice_of_planes)
							MPr_P[copy_of_plane].append([holonomy,sc,cumulative_angle_from_initial_slice])
			# These are the unambiguous vertices
			else:	
				for i in range(len(saddle_connections)):
					sc=saddle_connections[i]
					holonomy=sc.holonomy()
					holonomy_angle=ang(holonomy)
					# MPr_P[copy_of_plane].append([holonomy,sc,copy_of_plane*2*pi+ang(holonomy)])
					if holonomy_angle>=angle_of_slice_of_planes:
						cumulative_angle_from_initial_slice=2*pi*(copy_of_plane)+(holonomy_angle-angle_of_slice_of_planes)
						MPr_P[copy_of_plane].append([holonomy,sc,cumulative_angle_from_initial_slice])
					else:
						cumulative_angle_from_initial_slice=2*pi*(copy_of_plane+1)+(holonomy_angle-angle_of_slice_of_planes)
						MPr_P[copy_of_plane].append([holonomy,sc,cumulative_angle_from_initial_slice])
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

### NOT PICKING UP ALL THE INVERSE MATRICES OF THE SURFACE IN THE COMMENTS AT THE TOP

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
			angle_of_image_of_slice=ang(M*MPr[P][0][0][0])
			angle_of_original_slice=ang(first_edge_from_singularity[P])
			if angle_of_image_of_slice>=angle_of_original_slice:
				counterclockwise_angle_of_image_slice_from_original_slice=angle_of_image_of_slice-angle_of_original_slice
			else:
				counterclockwise_angle_of_image_slice_from_original_slice=2*pi+(angle_of_image_of_slice-angle_of_original_slice)
			for copy_of_plane in range(number_of_planes):
				for i in range(len(MPr[P][copy_of_plane])):
					mp=MPr[P][copy_of_plane][i][0]
					point=M*vector([mp[0],mp[1]])
					if point[0]^2+point[1]^2<=(2*rho)^2:
						angle_of_point=ang(point)
						if angle_of_point>=angle_of_original_slice:
							counterclockwise_angle_of_point_from_original_slice=angle_of_point-angle_of_original_slice
						else:
							counterclockwise_angle_of_point_from_original_slice=2*pi+(angle_of_point-angle_of_original_slice)
						# If the counter-clockwise angle of `point' from angle_of_image_of_slice is greater than or equal to counterclockwise_angle_of_image_slice_from_original_slice, then we assign this marked period to this copy_of_plane; otherwise,it gets assigned to the next copy_of_plane (modulo the number of copies of the plane)
						# We store this data as a list: [holonomy, preimage's corresponding saddle connection on X (as recorded in flatsurf), angle of holonomy from angle_of_original slice in a single copy of the plane, cumulative angle of preimage of this marked period from slice on 0th copy of plane]
						if counterclockwise_angle_of_point_from_original_slice>=counterclockwise_angle_of_image_slice_from_original_slice:
							F_M_MPr_P_bounded_by_2rho[copy_of_plane].append([point,MPr[P][copy_of_plane][i][1],counterclockwise_angle_of_point_from_original_slice,MPr[P][copy_of_plane][i][2]])
						else:
							F_M_MPr_P_bounded_by_2rho[(copy_of_plane+1)%number_of_planes].append([point,MPr[P][copy_of_plane][i][1],counterclockwise_angle_of_point_from_original_slice,MPr[P][copy_of_plane][i][2]])
			# Sort holonomies in each copy of the plane by counterclockwise angle from original slice in chosen copy of plane
			for copy_of_plane in range(len(MPr[P])):
				F_M_MPr_P_bounded_by_2rho[copy_of_plane].sort(key=lambda x:x[2])
			F_M_MPr_bounded_by_2rho.append(F_M_MPr_P_bounded_by_2rho)
		# Sort the previous list by cone angles of singularities so we may easily check this list (and permutations of it) against MP_2rho_sorted
		F_M_MPr_bounded_by_2rho_sorted=sorted(F_M_MPr_bounded_by_2rho,key=len)

		# Now we test various permutations of singularities of the same cone angle within F_M_MPr_bounded_by_2rho_sorted as well as cyclic permutations of copies of the plane corresponding to each singularity
		# We first check if the holonomies of a permutation of F_M_MPr_bounded_by_2rho_sorted match those of MP_2rho_sorted; if so, we then test to see if the $\mathbb{Z}_2$-action is the same
	# 	for permute_singularities in permutations_of_singularities.list():
	# 				F_M_first_copy=deepcopy(F_M_MPr_bounded_by_2rho_sorted)
	# 				place=0
	# 				for n in range(len(permutation_group_orders)):
	# 					num_sings_of_same_angle=permutation_group_orders[n]
	# 					F_M_first_copy[place:place+num_sings_of_same_angle]=permute_singularities[n](F_M_first_copy[place:place+num_sings_of_same_angle])
	# 					place+=num_sings_of_same_angle
	# 				for permute_planes in permutations_of_copies_of_plane.list():
	# 					F_M_copy=deepcopy(F_M_first_copy)
	# 					for P in range(X._num_singularities):
	# 						permute_planes_P=permute_planes[P]
	# 						F_M_copy[P]=permute_planes_P(F_M_copy[P])
	# 					M_is_a_keeper=True
	# 					for P0 in range(X._num_singularities):
	# 						for copy_of_plane0 in range(len(F_M_copy[P0])):
	# 							# Check that the holonomies associated to a particular singularity and copy of the plane match; if so, check the $\mathbb{Z}_2$-action
	# 							F_M_holonomies=[F_M_copy[P0][copy_of_plane0][i][0] for i in range(len(F_M_copy[P0][copy_of_plane0]))]
	# 							MP_2rho_holonomies=[MP_2rho_sorted[P0][copy_of_plane0][i][0] for i in range(len(MP_2rho_sorted[P0][copy_of_plane0]))]
	# 							if F_M_holonomies==MP_2rho_holonomies:
	# 								for i in range(len(F_M_copy[P0][copy_of_plane0])):
	# 									# if F_M_copy[P0][copy_of_plane0][i][0]==MP_2rho_sorted[P0][copy_of_plane0][i][0]:
	# 									F_M_saddle_connection_on_X=F_M_copy[P0][copy_of_plane0][i][1]
	# 									F_M_saddle_connection_preimage_angle=F_M_copy[P0][copy_of_plane0][i][3]
	# 									MPr_saddle_connection_on_X=MP_2rho_sorted[P0][copy_of_plane0][i][1]
	# 									MPr_saddle_connection_angle=MP_2rho_sorted[P0][copy_of_plane0][i][2]
	# 									F_M_saddle_connection_on_X_opposite=F_M_saddle_connection_on_X.invert()
	# 									MPr_saddle_connection_on_X_opposite=MPr_saddle_connection_on_X.invert()
	# 									F_M_found_opposite=False
	# 									MPr_found_opposite=False
	# 									for P1 in range(X._num_singularities):
	# 										for copy_of_plane1 in range(len(MPr[P1])):
	# 											for j in range(len(MPr[P1][copy_of_plane1])):
	# 												if MPr[P1][copy_of_plane1][j][1]==F_M_saddle_connection_on_X_opposite:
	# 													### THINK ABOUT THIS IF X HAS MULTIPLE SINGULARITIES OF DIFFERENT CONE ANGLES ###
	# 													F_M_angle_difference=(int((MPr[P1][copy_of_plane1][j][2]-F_M_saddle_connection_preimage_angle)/pi)%(2*len(MPr[P1])))*pi
	# 													F_M_Z2_index=[P1,F_M_angle_difference]
	# 													F_M_found_opposite=True
	# 												if MPr[P1][copy_of_plane1][j][1]==MPr_saddle_connection_on_X_opposite:
	# 													MPr_angle_difference=(int((MPr[P1][copy_of_plane1][j][2]-MPr_saddle_connection_angle)/pi)%(2*len(MPr[P1])))*pi
	# 													MPr_Z2_index=[P1,MPr_angle_difference]
	# 													MPr_found_opposite=True
	# 												if F_M_found_opposite==True and MPr_found_opposite==True:
	# 													if F_M_Z2_index!=MPr_Z2_index:
	# 														M_is_a_keeper=False
	# 													break
	# 											if M_is_a_keeper==False:
	# 												break
	# 										if M_is_a_keeper==False:
	# 											break
	# 									# else:
	# 									# 	M_is_a_keeper=False
	# 									# 	break
	# 								if M_is_a_keeper==False:
	# 									break
	# 						if M_is_a_keeper==False:
	# 							break
	# 					if M_is_a_keeper==True:
	# 						A_r.append(M)
	# 						break
	# 				if M_is_a_keeper==True:
	# 					break	
	# return(A_r)

		M_is_a_keeper=[[True for i in range(len(MP_2rho_sorted[P]))] for P in range(len(MP_2rho_sorted))]
		for permute_singularities in permutations_of_singularities.list():
					F_M_first_copy=deepcopy(F_M_MPr_bounded_by_2rho_sorted)
					place=0
					for n in range(len(permutation_group_orders)):
						num_sings_of_same_angle=permutation_group_orders[n]
						F_M_first_copy[place:place+num_sings_of_same_angle]=permute_singularities[n](F_M_first_copy[place:place+num_sings_of_same_angle])
						place+=num_sings_of_same_angle
					for permute_planes in permutations_of_copies_of_plane.list():
						# This will let us end a particular test if it is guaranteed not to work 
						break_flag=False
						M_is_a_keeper_flag=[[False for i in range(len(MP_2rho_sorted[P]))] for P in range(len(MP_2rho_sorted))]
						F_M_copy=deepcopy(F_M_first_copy)
						for P in range(X._num_singularities):
							permute_planes_P=permute_planes[P]
							F_M_copy[P]=permute_planes_P(F_M_copy[P])
						# M_is_a_keeper=True
						for P0 in range(X._num_singularities):
							for copy_of_plane0 in range(len(F_M_copy[P0])):
								# Check that the holonomies associated to a particular singularity and copy of the plane match; if so, check the $\mathbb{Z}_2$-action
								F_M_holonomies=[F_M_copy[P0][copy_of_plane0][i][0] for i in range(len(F_M_copy[P0][copy_of_plane0]))]
								MP_2rho_holonomies=[MP_2rho_sorted[P0][copy_of_plane0][i][0] for i in range(len(MP_2rho_sorted[P0][copy_of_plane0]))]
								if F_M_holonomies==MP_2rho_holonomies:
									for i in range(len(F_M_copy[P0][copy_of_plane0])):
										# if F_M_copy[P0][copy_of_plane0][i][0]==MP_2rho_sorted[P0][copy_of_plane0][i][0]:
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
														if F_M_Z2_index==MPr_Z2_index:
															M_is_a_keeper_flag[P0][copy_of_plane0]=True
														break
								if M_is_a_keeper_flag[P0][copy_of_plane0]==False:
									break_flag=True
									break
						if M_is_a_keeper_flag==M_is_a_keeper:
							A_r.append(M)
							break




					# 							if M_is_a_keeper==False:
					# 								break
					# 						if M_is_a_keeper==False:
					# 							break
					# 					# else:
					# 					# 	M_is_a_keeper=False
					# 					# 	break
					# 				if M_is_a_keeper==False:
					# 					break
					# 		if M_is_a_keeper==False:
					# 			break
					# 	if M_is_a_keeper==True:
					# 		A_r.append(M)
					# 		break
					# if M_is_a_keeper==True:
					# 	break	
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

# Takes as its input the Frobenius norm of a matrix M and outputs the minimum distance between I and the perpendicular bisector between I and M*I
def chi_2(x):
	return -ln(sqrt((1/2)*(x^2-sqrt(x^4-4))))

iteration=0
while ContainmentVolume==infinity and iteration<iteration_limit:
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

	x_min=-5
	x_max=5
	y_min=0
	y_max=5

	perp_bisector_plots=[]
	x=var('x')
	y=var('y')
	# When the cardinality of free_sides is finite, we'll know ContainmentVolume is finite
	free_sides=RealSet(-infinity,infinity)

	# These will record perpendicular bisectors between I and M*I for each M in Ar_mod_MinusIdentity
	# Each geodesic is recorded as [feet,polygon vertices,matrix,(1/2)d(I,M*I)] where... 
	# ...feet is a list [l,r] of left and right feet of the geodesic (where r=infinity if it is vertical);
	# ...polygon vertices is None unless it is verified that the geodesic is a side of the polygon formed by intersecting all halfplanes containing I, in which case we record vertices as [v0,v1] where v0 is the left-most vertex (or has smallest imaginary part if the geodesic is vertical);
	# ...matrix is the matrix M in M*I;
	# ...(1/2)d(I,M*I) is the minimum distance between I and the perpendicular bisector between I and M*I (this is giveen by chi_2(||M||); see pp. 98-100 of Edwards)
	vertical_geodesics=[]
	# Geodesics such that I is in the interior of the Euclidean semicircle defining the geodesic
	geodesics_enclosing_I=[]
	# Geodesics such that I is in the exterior of the Euclidean semicircle defining the geodesic
	geodesics_exposing_I=[]
	for M in Ar_mod_MinusIdentity:
		frobenius_norm=sqrt(M[0][0]^2+M[0][1]^2+M[1][0]^2+M[1][1]^2)
		if M!=matrix([[1,0],[0,1]]):
			M_act_on_I=(QQbar(M[0][0])*QQbar(I)+QQbar(M[0][1]))/(QQbar(M[1][0])*QQbar(I)+QQbar(M[1][1]))
			x_M=M_act_on_I.real()
			y_M=M_act_on_I.imag()
			perp_bisector_is_vertical=False
			# If M_act_on_I has imaginary part 1, then the perpendicular bisector of the geodesic connecting I and M_act_on_I is a vertical line
			if y_M==1:
				perp_bisector_is_vertical=True
				perp_bisector_l_foot=x_M/2
				perp_bisector_r_foot=infinity
			# If M_act_on_I is on the imaginary axis, we can find the perpendicular bisector of the geodesic connecting I and M_act_on_I immediately
			elif x_M==0:
				perp_bisector_l_foot=-y_M^(1/2)
				perp_bisector_r_foot=y_M^(1/2)
				perp_bisector_center=(0,0)
				perp_bisector_radius=(perp_bisector_r_foot-perp_bisector_l_foot)/2
			# Otherwise we find the perpendicular bisector of the geodesic connecting I and M_act_on_I by mapping this geodesic under a Mobius transformation $phi$ to the imaginary axis, finding the perpendicular bisector of this image, and then finding $phi^{-1}$ of this perpendicular bisector
			else:
				# Center and radius, respecively, of semicircular geodesic connecting I and M_act_on_I
				c=(x_M^2+y_M^2-1)/(2*x_M)
				R=(1+c^2)^(1/2)
				# Left foot of geodesic connecting I and M_act_on_I
				l=c-R
				# Imaginary part of hyperbolic midpoint between $I$ and $phi(M_act_on_I)$, where $phi$ is the Mobius transformation taking the semicircular geodesic connecting I and M_act_on_I with left and right feet l and r, resp., to imaginary axis where l -> 0, i -> i, and r -> infinity
				# The map $phi$ is given by $phi(z)=(z-l)/(l*z+1)$, and $phi^{-1}(z)=(z+l)/(-l*z+1)$
				# imag_mdpt is the squareroot of the imaginary part of $phi(M_act_on_I)$
				imag_mdpt=((y_M*(1+l^2))/(l^2*(x_M^2+y_M^2)+2*l*x_M+1))^(1/2)
				# The left and right feet of the perpendicular bisector of the geodesic between $I$ and $phi(M_act_on_I)$ are -imag_mdpt and imag_mdpt
				# The feet of the perpendicular bisector of the geodesic between I and M_act_on_I are given by $phi^{-1}(-imag_mdpt)$ and $phi^{-1}(imag_mdpt)$
				foot0=(-imag_mdpt+l)/(l*imag_mdpt+1)
				foot1=(imag_mdpt+l)/(-l*imag_mdpt+1)
				perp_bisector_center=((foot1+foot0)/2,0)
				perp_bisector_radius=((foot1-foot0).abs())/2
				if foot0<foot1:
					perp_bisector_l_foot=foot0
					perp_bisector_r_foot=foot1
				else:
					perp_bisector_l_foot=foot1
					perp_bisector_r_foot=foot0
			if perp_bisector_is_vertical==True:
				vertical_geodesics.append([[perp_bisector_l_foot,perp_bisector_r_foot],[[perp_bisector_l_foot,0],[perp_bisector_r_foot,0]],M,chi_2(frobenius_norm)])
				perp_bisector_plots.append(implicit_plot(x-perp_bisector_l_foot,(x,x_min,x_max),(y,y_min,y_max)))
				# If the vertical perpendicular bisector is left of i, then we remove $(-infinity,perp_bisector_l_foot)$ from free_sides
				if perp_bisector_l_foot<0:
					free_sides=free_sides.difference(-infinity,perp_bisector_l_foot)
				# If the vertical perpendicular bisector is right of i, then we remove $(perp_bisector_l_foot,infinity)$ from free_sides	
				else:
					free_sides=free_sides.difference(perp_bisector_l_foot,infinity)
			else:
				perp_bisector_plots.append(arc(perp_bisector_center,perp_bisector_radius,sector=(0,pi)))
				# If the geodesic [perp_bisector_l_foot,perp_bisector_r_foot] `encloses' i, then we remove $(-infinity,perp_bisector_l_foot)\cup (perp_bisector_r_foot,infinity)$ from free_sides
				if perp_bisector_l_foot*perp_bisector_r_foot<-1:
					geodesics_enclosing_I.append([[perp_bisector_l_foot,perp_bisector_r_foot],[[perp_bisector_l_foot,0],[perp_bisector_r_foot,0]],M,chi_2(frobenius_norm)])
					free_sides=free_sides.difference(-infinity,perp_bisector_l_foot)
					free_sides=free_sides.difference(perp_bisector_r_foot,infinity)
				# If the geodesic [perp_bisector_l_foot,perp_bisector_r_foot] does not `enclose' i, then we remove the interval $(perp_bisector_l_foot,perp_bisector_r_foot)$ from free_sides
				else:
					geodesics_exposing_I.append([[perp_bisector_l_foot,perp_bisector_r_foot],[[perp_bisector_l_foot,0],[perp_bisector_r_foot,0]],M,chi_2(frobenius_norm)])
					free_sides=free_sides.difference(perp_bisector_l_foot,perp_bisector_r_foot)
			
	free_sides_cardinality=free_sides.cardinality()
	if free_sides_cardinality in ZZ:
		ContainmentVolume='Is Finite'

	Dir_dom_plot=sum(perp_bisector_plots)
	iteration+=1


# This function constructs the hyperbolic polygon formed as the intersection of half-planes containing I defined by h-lines
def Omega(vertical_geods,geods_enclosing_I,geods_exposing_I):
	# This will be needed below to find all geodesics exposing I that form sides of polygon
	FreeSides=RealSet(-infinity,infinity)

	# Sort the geodesics by distance from I 
	vertical_geods_copy=sorted(vertical_geods,key=lambda x:x[3])
	geods_enclosing_I_copy=sorted(geods_enclosing_I,key=lambda x:x[3])
	geods_exposing_I_copy=sorted(geods_exposing_I,key=lambda x:x[3])

	polygon=[]

	# Find the geodesics enclosing I (if any) that are closest to I, and include them as sides of polygon
	if len(geods_enclosing_I_copy)>0:
		enclosing_min_distance=geods_enclosing_I_copy[0][3]
		for gamma in geods_enclosing_I:
			if gamma[3]==enclosing_min_distance:
				geods_enclosing_I_copy.remove(gamma)
				# Left and right feet of gamma
				left=gamma[0][0]
				right=gamma[0][1]
				# Center and radius of semicircle describing gamma
				center=(right+left)/2
				radius=(right-left)/2
				for polygon_side in polygon:
					# Check if polygon_side is vertical
					if polygon_side[0][1]==infinity:
						# Check if gamma intersects polygon_side
						if polygon_side[0][0]>left and polygon_side[0][0]<right:
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							vertex_x=polygon_side[0][0]
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_y>polygon_side[1][0][0] and vertex_y<polygon_side[1][1][0]:
								# The 'upper-most' vertex of polygon_side becomes this intersection
								polygon_side[1][1]=[vertex_x,vertex_y]
								# If the intersection is left of the imaginary axis, then the 'left-most' vertex of gamma becomes this intersection
								if vertex_x<0:
									gamma[1][0]=[vertex_x,vertex_y]
								# If the intersection is right of the imaginary axis, then the 'right-most' vertex of gamma becomes this intersection
								else:
									gamma[1][1]=[vertex_x,vertex_y]
					# Consider when polygon_side is not vertical
					else: 
						l_polygon_side=polygon_side[0][0]
						r_polygon_side=polygon_side[0][1]
						# Check if polygon_side encloses I
						if l_polygon_side*r_polygon_side<-1:
							# Check if gamma and polygon_side intersect
							if (l_polygon_side>left and l_polygon_side<right and r_polygon_side>right):
								# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
								center_polygon_side=(r_polygon_side+l_polygon_side)/2
								radius_polygon_side=(r_polygon_side-l_polygon_side)/2
								vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
								vertex_y=sqrt(radius^2-(vertex_x-center)^2)
								# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
								if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:								
									# The 'right-most' vertex of polygon_side becomes this intersection
									polygon_side[1][1]=[vertex_x,vertex_y]
									# The 'left-most' vertex of gamma becomes this intersection
									gamma[1][0]=[vertex_x,vertex_y]
							elif (l_polygon_side<left and r_polygon_side>left and r_polygon_side<right):
								# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
								center_polygon_side=(r_polygon_side+l_polygon_side)/2
								radius_polygon_side=(r_polygon_side-l_polygon_side)/2
								vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
								vertex_y=sqrt(radius^2-(vertex_x-center)^2)
								# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
								if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:								
									# The 'left-most' vertex of polygon_side becomes this intersection
									polygon_side[1][0]=[vertex_x,vertex_y]
									# The 'right-most' vertex of gamma becomes this intersection
									gamma[1][1]=[vertex_x,vertex_y]
						# Otherwise polygon_side exposes I
						else:
							# Check if gamma and polygon_side intersect
							if (l_polygon_side>left and l_polygon_side<right and r_polygon_side>right):
								# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
								center_polygon_side=(r_polygon_side+l_polygon_side)/2
								radius_polygon_side=(r_polygon_side-l_polygon_side)/2
								vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
								vertex_y=sqrt(radius^2-(vertex_x-center)^2)
								# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
								if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:								
									# The 'right-most' vertex of polygon_side becomes this intersection
									polygon_side[1][1]=[vertex_x,vertex_y]
									# The 'right-most' vertex of gamma becomes this intersection
									gamma[1][1]=[vertex_x,vertex_y]
							elif (l_polygon_side<left and r_polygon_side>left and r_polygon_side<right):
								# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
								center_polygon_side=(r_polygon_side+l_polygon_side)/2
								radius_polygon_side=(r_polygon_side-l_polygon_side)/2
								vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
								vertex_y=sqrt(radius^2-(vertex_x-center)^2)
								# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
								if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:								
									# The 'left-most' vertex of polygon_side becomes this intersection
									polygon_side[1][0]=[vertex_x,vertex_y]
									# The 'left-most' vertex of gamma becomes this intersection
									gamma[1][0]=[vertex_x,vertex_y]
				FreeSides=FreeSides.difference(-infinity,left)
				FreeSides=FreeSides.difference(right,infinity)															
				polygon.append(gamma)

	# Find the geodesics exposing I that are closest to I, and include them as sides of polygon
	exposing_min_distance=geods_exposing_I_copy[0][3]
	for gamma in geods_exposing_I:
		if gamma[3]==exposing_min_distance:
			geods_exposing_I_copy.remove(gamma)
			# Left and right feet of gamma
			left=gamma[0][0]
			right=gamma[0][1]
			# Center and radius of semicircle describing gamma
			center=(right+left)/2
			radius=(right-left)/2
			for polygon_side in polygon:
				# Check if polygon_side is vertical
				if polygon_side[0][1]==infinity:
					# Check if gamma intersects polygon_side
					if polygon_side[0][0]>left and polygon_side[0][0]<right:
						# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
						vertex_x=polygon_side[0][0]
						vertex_y=sqrt(radius^2-(vertex_x-center)^2)
						# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
						if vertex_y>polygon_side[1][0][0] and vertex_y<polygon_side[1][1][0]:
							# The 'lower-most' vertex of polygon_side becomes this intersection
							polygon_side[1][0]=[vertex_x,vertex_y]
							# If the intersection is left of the imaginary axis, then the 'left-most' vertex of gamma becomes this intersection
							if vertex_x<0:
								gamma[1][0]=[vertex_x,vertex_y]
							# If the intersection is right of the imaginary axis, then the 'right-most' vertex of gamma becomes this intersection
							else:
								gamma[1][1]=[vertex_x,vertex_y]
				# Consider when polygon_side is not vertical
				else: 
					l_polygon_side=polygon_side[0][0]
					r_polygon_side=polygon_side[0][1]
					# Check if polygon_side encloses I
					if l_polygon_side*r_polygon_side<-1:
						# Check if gamma and polygon_side intersect
						if (l_polygon_side>left and l_polygon_side<right and r_polygon_side>right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:												
								# The 'left-most' vertex of polygon_side becomes this intersection
								polygon_side[1][0]=[vertex_x,vertex_y]
								# The 'left-most' vertex of gamma becomes this intersection
								gamma[1][0]=[vertex_x,vertex_y]
						elif (l_polygon_side<left and r_polygon_side>left and r_polygon_side<right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:											
								# The 'right-most' vertex of polygon_side becomes this intersection
								polygon_side[1][1]=[vertex_x,vertex_y]
								# The 'right-most' vertex of gamma becomes this intersection
								gamma[1][1]=[vertex_x,vertex_y]
					# Otherwise polygon_side exposes I
					else:
						# Check if gamma and polygon_side intersect
						if (l_polygon_side>left and l_polygon_side<right and r_polygon_side>right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:														
								# The 'left-most' vertex of polygon_side becomes this intersection
								polygon_side[1][0]=[vertex_x,vertex_y]
								# The 'right-most' vertex of gamma becomes this intersection
								gamma[1][1]=[vertex_x,vertex_y]
						elif (l_polygon_side<left and r_polygon_side>left and r_polygon_side<right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:																
								# The 'right-most' vertex of polygon_side becomes this intersection
								polygon_side[1][1]=[vertex_x,vertex_y]
								# The 'left-most' vertex of gamma becomes this intersection
								gamma[1][0]=[vertex_x,vertex_y]						
			FreeSides=FreeSides.difference(left,right)
			polygon.append(gamma)

	# Find the vertical sides contributing to polygon, if any
	if len(vertical_geods_copy)>=2:
		closest_vertical_geods=vertical_geods_copy[:2]
		for gamma in closest_vertical_geods:
			append_to_polygon=False
			# The finite foot determining gamma
			left=gamma[0][0]
			# Check if geods_enclosing_I is not empty; if this is the case, then gamma must intersect some segment contributing to the side polygon of some polygon_side enclosing I in order for gamma to be included in polygon (if not, then gamma lies too far left or right of polygon to contribute to its sides)
			if len(geods_enclosing_I)>0:
				for polygon_side in polygon:
					# Check if gamma intersects the segment of polygon_side contributing to polygon
					if left>polygon_side[1][0][0] and left<polygon_side[1][1][0]:
						l_polygon_side=polygon_side[0][0]
						r_polygon_side=polygon_side[0][1]
						center_polygon_side=(r_polygon_side+l_polygon_side)/2
						radius_polygon_side=(r_polygon_side-l_polygon_side)/2
						vertex_x=left
						vertex_y=sqrt(radius_polygon_side^2-(left-center_polygon_side)^2)
						# Check if polygon_side encloses I
						if l_polygon_side*r_polygon_side<-1:
							append_to_polygon=True
							# If gamma is left of the imaginary axis, then the 'left-most' vertex of polygon_side and the 'upper-most' vertex of gamma become the intersection of gamma and polygon_side
							if left<0:
								FreeSides=FreeSides.difference(-infinity,left)
								polygon_side[1][0]=[vertex_x,vertex_y]
								gamma[1][1]=[vertex_x,vertex_y]
							# If gamma is right of the imaginary axis, then the 'right-most' vertex of polygon_side and the 'upper-most' vertex of gamma become the intersection of gamma and polygon_side
							else:
								FreeSides=FreeSides.difference(left,infinity)
								polygon_side[1][1]=[vertex_x,vertex_y]
								gamma[1][1]=[vertex_x,vertex_y]
						# Otherwise polygon_side exposes I
						else:
							# If gamma is left of the imaginary axis, then the 'left-most' vertex of polygon_side and the 'lower-most' vertex of gamma become the intersection of gamma and polygon_side
							if left<0:
								FreeSides=FreeSides.difference(-infinity,left)								
								polygon_side[1][0]=[vertex_x,vertex_y]
								gamma[1][0]=[vertex_x,vertex_y]
							# If gamma is right of the imaginary axis, then the 'right-most' vertex of polygon_side and the 'lower-most' vertex of gamma become the intersection of gamma and polygon_side
							else:
								FreeSides=FreeSides.difference(left,infinity)
								polygon_side[1][1]=[vertex_x,vertex_y]
								gamma[1][0]=[vertex_x,vertex_y]
			# If geods_enclosing_I is empty, then gamma is included in polygon by default (otherwise the polygon would have infinite area as it would always have free sides)
			else:
				append_to_polygon=True
				for polygon_side in polygon:
					# Check if gamma intersects the segment of polygon_side contributing to polygon
					if left>polygon_side[1][0][0] and l<polygon_side[1][1][0]:
						l_polygon_side=polygon_side[0][0]
						r_polygon_side=polygon_side[0][1]
						center_polygon_side=(r_polygon_side+l_polygon_side)/2
						radius_polygon_side=(r_polygon_side-l_polygon_side)/2
						vertex_x=left
						vertex_y=sqrt(radius_polygon_side^2-(left-center_polygon_side)^2)
						# Here each polygon_side exposes I by assumption
						# If gamma is left of the imaginary axis, then the 'left-most' vertex of polygon_side and the 'lower-most' vertex of gamma become the intersection of gamma and polygon_side
						if left<0:
							FreeSides=FreeSides.difference(-infinity,left)							
							polygon_side[1][0]=[vertex_x,vertex_y]
							gamma[1][0]=[vertex_x,vertex_y]
						# If gamma is right of the imaginary axis, then the 'right-most' vertex of polygon_side and the 'lower-most' vertex of gamma become the intersection of gamma and polygon_side
						else:
							FreeSides=FreeSides.difference(left,infinity)							
							polygon_side[1][1]=[vertex_x,vertex_y]
							gamma[1][0]=[vertex_x,vertex_y]											
			if append_to_polygon==True:
				polygon.append(gamma)

	# Now we see which other geodesics might contribute to the sides of polygon; these sides are not the closest sides of polygon to I, as those were all found above
	while FreeSides.cardinality()==infinity:
		geods_enclosing_I_copy2=deepcopy(geods_enclosing_I_copy)
		for gamma in geods_enclosing_I_copy2:
			append_to_polygon=False
			# Left and right feet of gamma
			left=gamma[0][0]
			right=gamma[0][1]
			# Center and radius of semicircle describing gamma
			center=(right+left)/2
			radius=(right-left)/2
			for polygon_side in polygon:
				# Check if polygon_side is vertical
				if polygon_side[0][1]==infinity:
					# Check if gamma intersects polygon_side
					if polygon_side[0][0]>left and polygon_side[0][0]<right:
						# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
						vertex_x=polygon_side[0][0]
						vertex_y=sqrt(radius^2-(vertex_x-center)^2)
						# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
						if vertex_y>polygon_side[1][0][0] and vertex_y<polygon_side[1][1][0]:
							geods_enclosing_I_copy.remove(gamma)
							append_to_polygon=True
							# The 'upper-most' vertex of polygon_side becomes this intersection
							polygon_side[1][1]=[vertex_x,vertex_y]
							# If the intersection is left of the imaginary axis, then the 'left-most' vertex of gamma becomes this intersection
							if vertex_x<0:
								gamma[1][0]=[vertex_x,vertex_y]
							# If the intersection is right of the imaginary axis, then the 'right-most' vertex of gamma becomes this intersection
							else:
								gamma[1][1]=[vertex_x,vertex_y]
				# Consider when polygon_side is not vertical
				else: 
					l_polygon_side=polygon_side[0][0]
					r_polygon_side=polygon_side[0][1]
					# Check if polygon_side encloses I
					if l_polygon_side*r_polygon_side<-1:
						# Check if gamma and polygon_side intersect
						if (l_polygon_side>left and l_polygon_side<right and r_polygon_side>right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:
								geods_enclosing_I_copy.remove(gamma)
								append_to_polygon=True
								# The 'right-most' vertex of polygon_side becomes this intersection
								polygon_side[1][1]=[vertex_x,vertex_y]
								# The 'left-most' vertex of gamma becomes this intersection
								gamma[1][0]=[vertex_x,vertex_y]
						elif (l_polygon_side<left and r_polygon_side>left and r_polygon_side<right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:
								geods_enclosing_I_copy.remove(gamma)
								append_to_polygon=True
								# The 'left-most' vertex of polygon_side becomes this intersection
								polygon_side[1][0]=[vertex_x,vertex_y]
								# The 'right-most' vertex of gamma becomes this intersection
								gamma[1][1]=[vertex_x,vertex_y]
					# Otherwise polygon_side exposes I
					else:
						# Check if gamma and polygon_side intersect
						if (l_polygon_side>left and l_polygon_side<right and r_polygon_side>right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:
								geods_enclosing_I_copy.remove(gamma)
								append_to_polygon=True
								# The 'right-most' vertex of polygon_side becomes this intersection
								polygon_side[1][1]=[vertex_x,vertex_y]
								# The 'right-most' vertex of gamma becomes this intersection
								gamma[1][1]=[vertex_x,vertex_y]
						elif (l_polygon_side<left and r_polygon_side>left and r_polygon_side<right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:
								geods_enclosing_I_copy.remove(gamma)
								append_to_polygon=True
								# The 'left-most' vertex of polygon_side becomes this intersection
								polygon_side[1][0]=[vertex_x,vertex_y]
								# The 'left-most' vertex of gamma becomes this intersection
								gamma[1][0]=[vertex_x,vertex_y]		
			if append_to_polygon==True:
				polygon.append(gamma)

		# Here we find geodesics gamma exposing I that do not intersect an already existing side of polygon, but such that the interval between the feet of gamma is contained in a free side of the polygon, and this interval is maximal in the sense that for all other gamma' satisfying this condition, the interval beween the feet of gamma' does not contain the interval between the feet of gamma
		for i in range(FreeSides.n_components()):
			geods_exposing_I_copy2=deepcopy(geods_exposing_I_copy)
			# We find the i^th connected component of FreeSides and make sure it is not a singleton
			J=FreeSides.get_interval(i)
			J_dummy=RealSet(J)
			if J_dummy.cardinality()==infinity:
				lower_bound=J.lower()
				# Find all geodesics gamma whose left foot is at the lower bound of the interval J of FreeSides
				gammas_with_left_foot_at_lower_bound=[]
				for gamma in geods_exposing_I_copy2:
				# for gamma in geods_exposing_I_copy:
					if gamma[0][0]==lower_bound:
						gammas_with_left_foot_at_lower_bound.append(gamma)
				if len(gammas_with_left_foot_at_lower_bound)>0:
					# Sort these by right foot
					gammas_with_left_foot_at_lower_bound.sort(key=lambda x:x[0][1])
					# The geodesic with largest right foot is maximal in the sense described above and will contribute to the side of polygon
					gamma=gammas_with_left_foot_at_lower_bound[-1]
					FreeSides=FreeSides.difference(gamma[0][0],gamma[0][1])
					geods_exposing_I_copy.remove(gamma)
					polygon.append(gamma)

		geods_exposing_I_copy2=deepcopy(geods_exposing_I_copy)
		for gamma in geods_exposing_I_copy2:
			append_to_polygon=False
			# Left and right feet of gamma
			left=gamma[0][0]
			right=gamma[0][1]
			# Center and radius of semicircle describing gamma
			center=(right+left)/2
			radius=(right-left)/2
			for polygon_side in polygon:
				# Check if polygon_side is vertical
				if polygon_side[0][1]==infinity:
					# Check if gamma intersects polygon_side
					if polygon_side[0][0]>left and polygon_side[0][0]<right:
						# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
						vertex_x=polygon_side[0][0]
						vertex_y=sqrt(radius^2-(vertex_x-center)^2)
						# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
						if vertex_y>polygon_side[1][0][0] and vertex_y<polygon_side[1][1][0]:
							FreeSides=FreeSides.difference(left,right)
							geods_exposing_I_copy.remove(gamma)
							append_to_polygon=True
							# The 'lower-most' vertex of polygon_side becomes this intersection
							polygon_side[1][0]=[vertex_x,vertex_y]
							# If the intersection is left of the imaginary axis, then the 'left-most' vertex of gamma becomes this intersection
							if vertex_x<0:
								gamma[1][0]=[vertex_x,vertex_y]
							# If the intersection is right of the imaginary axis, then the 'right-most' vertex of gamma becomes this intersection
							else:
								gamma[1][1]=[vertex_x,vertex_y]
				# Consider when polygon_side is not vertical
				else: 
					l_polygon_side=polygon_side[0][0]
					r_polygon_side=polygon_side[0][1]
					# Check if polygon_side encloses I
					if l_polygon_side*r_polygon_side<-1:
						# Check if gamma and polygon_side intersect
						if (l_polygon_side>left and l_polygon_side<right and r_polygon_side>right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:
								FreeSides=FreeSides.difference(left,right)
								geods_exposing_I_copy.remove(gamma)								
								append_to_polygon=True
								# The 'left-most' vertex of polygon_side becomes this intersection
								polygon_side[1][0]=[vertex_x,vertex_y]
								# The 'left-most' vertex of gamma becomes this intersection
								gamma[1][0]=[vertex_x,vertex_y]
						elif (l_polygon_side<left and r_polygon_side>left and r_polygon_side<right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:
								FreeSides=FreeSides.difference(left,right)
								geods_exposing_I_copy.remove(gamma)								
								append_to_polygon=True
								# The 'right-most' vertex of polygon_side becomes this intersection
								polygon_side[1][1]=[vertex_x,vertex_y]
								# The 'right-most' vertex of gamma becomes this intersection
								gamma[1][1]=[vertex_x,vertex_y]
					# Otherwise polygon_side exposes I
					else:
						# Check if gamma and polygon_side intersect
						if (l_polygon_side>left and l_polygon_side<right and r_polygon_side>right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:
								FreeSides=FreeSides.difference(left,right)
								geods_exposing_I_copy.remove(gamma)								
								append_to_polygon=True								
								# The 'left-most' vertex of polygon_side becomes this intersection
								polygon_side[1][0]=[vertex_x,vertex_y]
								# The 'right-most' vertex of gamma becomes this intersection
								gamma[1][1]=[vertex_x,vertex_y]
						elif (l_polygon_side<left and r_polygon_side>left and r_polygon_side<right):
							# Find the intersection vertex_x+I*vertex_y of gamma and polygon_side
							center_polygon_side=(r_polygon_side+l_polygon_side)/2
							radius_polygon_side=(r_polygon_side-l_polygon_side)/2
							vertex_x=(radius_polygon_side^2-center_polygon_side^2-(radius^2-center^2))/(2*(center-center_polygon_side))
							vertex_y=sqrt(radius^2-(vertex_x-center)^2)
							# Check that the intersection is in the interior of the segment of polygon_side contributing to the side of the polygon
							if vertex_x>polygon_side[1][0][0] and vertex_x<polygon_side[1][1][0]:
								FreeSides=FreeSides.difference(left,right)
								geods_exposing_I_copy.remove(gamma)								
								append_to_polygon=True							
								# The 'right-most' vertex of polygon_side becomes this intersection
								polygon_side[1][1]=[vertex_x,vertex_y]
								# The 'left-most' vertex of gamma becomes this intersection
								gamma[1][0]=[vertex_x,vertex_y]	
			if append_to_polygon==True:
				polygon.append(gamma)

	return polygon

a=Omega(vertical_geodesics,geodesics_enclosing_I,geodesics_exposing_I)
b=[a[i][0] for i in range(len(a))]
c=[[((b[i][0]+b[i][1])/2,0),(b[i][1]-b[i][0])/2] for i in range(len(a))]
arcs=[arc(c[i][0],c[i][1],sector=(0,pi)) for i in range(len(a))]

########################################################################################################
########################################################################################################