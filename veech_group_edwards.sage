r''' 
ALGORITHM 7.1
Takes as input translation surface (X,omega) such that $\Gamma(X,omega)\cap \text{SO}(2,\mathbb{R})\subseteq \{\pm \text{Id}\}$

X must be a translation surface from the flatsurf package
Example: 
sage: from flatsurf import *
sage: X=translation_surfaces.veech_2n_gon(5)
'''

# Will eventually uncomment this function when all of the subsections of the algorithm work
# def veech_group(X):
from flatsurf.geometry.polygon import polygons,PolygonPosition
assert isinstance(X,TranslationSurface)


###############################################
## 0: Compute the Voronoi decomposition of X ##
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


##################################################################
## 1: Calculate $\rho=\rho(X,\omega)$ using the Voronoi 2-cells ##
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


################################################################
## 2: Let $r=2\rho$ and let $b=\chi_1^{-1}(2\rho/r)=\sqrt{2}$ ##
################################################################
def chi_1_inv(x):
	return sqrt(x^2+x^(-2))

r=2*rho
b=sqrt(2)


##########################################################
## 3: Calculate $MP_P^r(X,\omega)$ for all $P\in\Sigma$ ##
##########################################################
# This records information to determine which disjoint copy of the plane an element of $MP_P^r$ should lie in
# vertices_in_disjoint_planes[P][n] gives a list of vertex data (polygon,vertex) such that separatrices leaving these vertices should develop in the n^th copy of the plane corresponding to singularity P
# There is ambiguity for some vertices, such as vertex 5 of the regular hexagon...these ambiguous vertices where separatrices may develop in either the n^th or (n+1)^st plane lead the list corresponding to the (n+1)^st plane
vertices_in_disjoint_planes=[]
first_edge_from_singularity=[]
for P in range(X._num_singularities):
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

	while copy_of_plane<X.angles()[P]:
		vertex_data=X.opposite_edge(vertex_data[0],(vertex_data[1]-1)%X.polygon(vertex_data[0]).num_edges())
		P_angle+=X.polygon(vertex_data[0]).angle(vertex_data[1])
		if P_angle>copy_of_plane+1:
			vertices_in_disjoint_planes_P.append(vertices_in_disjoint_planes_P_copy_of_plane)
			vertices_in_disjoint_planes_P_copy_of_plane=[]
			copy_of_plane+=1
		vertices_in_disjoint_planes_P_copy_of_plane.append(vertex_data)
	vertices_in_disjoint_planes.append(vertices_in_disjoint_planes_P)


# marked_periods(radius)[P][copy_of_plane] is the set of marked periods from singularity P of X bounded by radius that correspond to the (copy_of_plane)^th copy of the plane
### To do: create optional parameter to add to previously computed marked_periods(r0) ###
def marked_periods(radius):
	MPr=[]
	for P in range(X._num_singularities):
		MPr_P=[[] for i in range(X.angles()[P])]
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
					# If sc leaves its singularity clockwise of where the current copy of the plane begins, then its holonomy vector is assigned to the previous copy of the plane
					if det(M)<0:
						MPr_P[copy_of_plane-1].append([sc.holonomy(),sc])
					# Otherwise its holonomy vector is assigned to the current copy of the plane
					else:
						MPr_P[copy_of_plane].append([sc.holonomy(),sc])

			# These are the unambiguous vertices
			else:	
				for i in range(len(saddle_connections)):
					sc=saddle_connections[i]
					MPr_P[copy_of_plane].append([sc.holonomy(),sc])
		MPr.append(MPr_P)
	return MPr




#############################################################################################################################
## 4: Calculate $A_r=\{M\in\Gamma(X,\omega)\ |\ ||M||\le b\}=\text{SO}(2,\mathbb{R})\cap\Gamma(X,\omega)$ using Theorem 18 ##
#############################################################################################################################
### To do: create an optional parameter to add to previously computed A_r(radius0,norm_bound0) ###


def A_r(radius,norm_bound):	
	# matrices_to_check will be a list comprised of matrices in $\text{SL}_2\mathbb{Z}$ whose Frobenius norms are bounded above by norm_bound and which take a basis of vectors in marked_periods(radius) to another basis of vectors in marked_periods(radius)
	# MPr_bounded_by_2rho will be the list of marked periods bounded by $2*\rho$
	matrices_to_check=[]
	MPr=marked_periods(radius)
	MP_2rho=marked_periods(2*rho)
	MPr_bounded_by_2rho=[]
	for P in range(X._num_singularities):
		MPr_P_bounded_by_2rho=[]
		for copy_of_plane in range(len(MPr[P])):
			# Each set of holonomies in a copy of the plane is `dictionary' sorted, so vectors are listed from bottom left-most to top right-most
			MP_2rho[P][copy_of_plane].sort()
			v0=MP_2rho[P][copy_of_plane][0][0]
			v1=MP_2rho[P][copy_of_plane][1][0]
			T=matrix([[v0[0],v1[0]],[v0[1],v1[1]]])
			for k in range(len(MPr[P][copy_of_plane])):
				w0=MPr[P][copy_of_plane][k][0]
				for l in range(len(MPr[P][copy_of_plane])):
					w1=MPr[P][copy_of_plane][l][0]
					M=matrix([[w0[0],w1[0]],[w0[1],w1[1]]])
					if M.det()!=0:
						N=M*T.inverse()
						frobenius_norm_squared=N[0][0]^2+N[0][1]^2+N[1][0]^2+N[1][1]^2
						if N not in matrices_to_check and N.det()==1 and frobenius_norm_squared<=norm_bound^2:
							matrices_to_check.append(N)
		# Sort MP_2rho by cone angles of singularities so we may easily check this list against (permutations of) F_M_MPr_bounded_by_2rho below
		MPr_bounded_by_2rho_sorted=sorted(MP_2rho,key=len)

	# We will apply a matrix M to each point in MPr and keep only those bounded by $2*\rho$; we call this F_M_MPr_bounded_by_2rho
	A_r=[]
	for M in matrices_to_check:
		F_M_MPr_bounded_by_2rho=[]
		for P in range(X._num_singularities):
			F_M_MPr_P_bounded_by_2rho=[]
			for copy_of_plane in range(len(MPr[P])):
				F_M_MPr_P_copy_of_plane_bounded_by_2rho=[]
				for i in range(len(MPr[P][copy_of_plane])):
					mp=MPr[P][copy_of_plane][i][0]
					point=M*vector([mp[0],mp[1]])
					if point[0]^2+point[1]^2<=(2*rho)^2:
						F_M_MPr_P_copy_of_plane_bounded_by_2rho.append([point,MPr[P][copy_of_plane][i][1]])
				# Each set of holonomies in a copy of the plane is `dictionary' sorted, so vectors are listed from bottom left-most to top right-most
				F_M_MPr_P_bounded_by_2rho.append(sorted(F_M_MPr_P_copy_of_plane_bounded_by_2rho))
			F_M_MPr_bounded_by_2rho.append(F_M_MPr_P_bounded_by_2rho)
			# Sort the previous list by cone angles of singularities so we may easily check this list (and permutations of it) against MPr_bounded_by_2rho_sorted
			F_M_MPr_bounded_by_2rho_sorted=sorted(F_M_MPr_bounded_by_2rho)

		# Here we test various permutations of singularities of like cone angles within F_M_MPr_bounded_by_2rho_sorted as well as cyclic permutations of copies of the plane corresponding to each singularity
		# We first check if the holonomies of a permutation of F_M_MPr_bounded_by_2rho_sorted match those of MPr_bounded_by_2rho_sorted; ...
		# ...if so, we then test to see if the $\mathbb{Z}_2$-action is the same

		### NEED TO THINK MORE ABOUT PERMUTATIONS OF COPIES OF THE PLANE; AS IS, THE CODE DOESN'T KEEP M=[[1,0],[2,1]] WHEN X IS THE L-SURFACE ###

		# List of sizes of cone angles after sorted
		cone_angles=[len(MPr_bounded_by_2rho_sorted[P]) for P in range(X._num_singularities)]
		# Records number of singularities of a particular order; will create symmetric groups on these numbers of elements and later take cartesian product of these groups for various permutations of singularities
		permutation_group_orders=[]
		angle_previous=None
		for i in range(len(cone_angles)):
			angle=cone_angles[i]
			if angle!=angle_previous:
				permutation_group_orders.append(cone_angles.count(angle))
			angle_previous=angle
		# Cartesian product of symmetric groups of orders from permutation_group_orders; will be used to permute singularities of the same order
		permutations_of_singularities=SymmetricGroup(permutation_group_orders[0])
		for i in permutation_group_orders[1:]:
			permutations_of_singularities=permutations_of_singularities.cartesian_product(SymmetricGroup(i))
		# Product of cyclic groups from orders of cone angles; will be used to check permutations of copies of the plane associated to a particular singularity
		permutations_of_copies_of_plane=CyclicPermutationGroup(cone_angles[0])
		for j in cone_angles[1:]:
			permutations_of_copies_of_plane=permutations_of_copies_of_plane.cartesian_product(CyclicPermutationGroup(j))
		# Permute singularities and copies of planes of F_M_MPr_bounded_by_2rho_sorted and see if the sets of holomonies (associated with a singularity of specific order and ordering of copies of planes) matches the set of holonomies of MPr_bounded_by_2rho_sorted
		for permute_singularities in permutations_of_singularities.list():
			F_M_first_copy=deepcopy(F_M_MPr_bounded_by_2rho_sorted)
			place=0
			for n in range(len(permutation_group_orders)):
				num_sings_of_same_angle=permutation_group_orders[n]
				F_M_first_copy[place:place+num_sings_of_same_angle]=permute_singularities(F_M_first_copy[place:place+num_sings_of_same_angle])
				place+=n
			for permute_planes in permutations_of_copies_of_plane.list():
				F_M_copy=deepcopy(F_M_first_copy)
				if X._num_singularities==1:
					F_M_copy[0]=permute_planes(F_M_copy[0])
				else:
					for P in range(X._num_singularities):
						permute_planes_P=permute_planes[P]
						F_M_copy[P]=permute_planes_P(F_M_copy[P])
				M_is_a_keeper=True
				for P0 in range(X._num_singularities):
					for copy_of_plane0 in range(len(F_M_copy[P0])):
						for i in range(len(F_M_copy[P0][copy_of_plane0])):
							# Check that the holonomies associated to a particular singularity and copy of the plane match; if so, check the $\mathbb{Z}_2$-action
							if F_M_copy[P0][copy_of_plane0][i][0]==MPr_bounded_by_2rho_sorted[P0][copy_of_plane0][i][0]:
								F_M_saddle_connection_on_X=F_M_copy[P0][copy_of_plane0][i][1]
								MPr_saddle_connection_on_X=MPr_bounded_by_2rho_sorted[P0][copy_of_plane0][i][1]
								F_M_saddle_connection_on_X_opposite=F_M_saddle_connection_on_X.invert()
								MPr_saddle_connection_on_X_opposite=MPr_saddle_connection_on_X.invert()
								F_M_found_opposite=False
								MPr_found_opposite=False
								for P1 in range(X._num_singularities):
									for copy_of_plane1 in range(len(F_M_copy[P1])):
										for j in range(len(F_M_copy[P1][copy_of_plane1])):
											if F_M_copy[P1][copy_of_plane1][j][1]==F_M_saddle_connection_on_X_opposite:
												F_M_Z2_index=[P1,copy_of_plane1,j]
												F_M_found_opposite=True
											if MPr_bounded_by_2rho_sorted[P1][copy_of_plane1][j][1]==MPr_saddle_connection_on_X_opposite:
												MPr_Z2_index=[P1,copy_of_plane1,j]
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
	return(A_r)
