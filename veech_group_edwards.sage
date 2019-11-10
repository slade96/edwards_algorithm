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
### ALSO NEED TO KEEP TRACK OF $\mathbb{Z}_2$-ACTION ###
### To do: create optional parameter to add to previously computed marked_periods(r0) ###
### Can maybe speed this up by computing X.saddle_connections(radius^2,polygon,vertex) and using direction and length associated to each entry of outputted list; this might not use exact arithmetic though ###
# def marked_periods(radius):
# 	MPr=[]
# 	for P in range(X._num_singularities):
# 		MPr_P=[[] for i in range(X.angles()[P])]
# 		for vertex_data in vertex_equivalence_classes[P]:
# 			copy_of_plane=[n for n in range(len(vertices_in_disjoint_planes[P])) if vertex_data in vertices_in_disjoint_planes[P][n]][0]
# 			polygon=vertex_data[0]
# 			vertex=vertex_data[1]
# 			saddle_connections=X.saddle_connections(radius^2,polygon,vertex)
# 			saddle_connections_directions=[sc.direction() for sc in saddle_connections]
# 			poly=X.polygon(polygon)
# 			point=poly.vertex(vertex)
# 			V=poly.parent().vector_space()
# 			saddle_connection_endpoints=[]
# 			for direction in saddle_connections_directions:
# 				is_vertex=0
# 				pnt=point
# 				polygon2=polygon
# 				poly2=poly
# 				holonomy=(0,0)
# 				while is_vertex==0:
# 					flow_to_bdry=poly2.flow_to_exit(V(pnt),V(direction))
# 					holonomy=V(holonomy)+V(flow_to_bdry[0])-V(pnt)
# 					if flow_to_bdry[1].is_vertex()==True:
# 						is_vertex=1
# 						saddle_connection_endpoints.append(holonomy)
# 					else:
# 						continue_flow_from=X.opposite_edge(polygon2,flow_to_bdry[1].get_edge())
# 						pnt=X.edge_transformation(polygon2,flow_to_bdry[1].get_edge())(flow_to_bdry[0])
# 						polygon2=continue_flow_from[0]
# 						poly2=X.polygon(polygon2)
# 			if vertex_data==vertices_in_disjoint_planes[P][copy_of_plane][0]:
# 				for i in range(len(saddle_connections_directions)):
# 					v0=first_edge_from_singularity[P]
# 					v1=saddle_connections_directions[i]
# 					M=matrix([[v0[0],v1[0]],[v0[1],v1[1]]])
# 					if det(M)<0:
# 						MPr_P[copy_of_plane-1].append(saddle_connection_endpoints[i])
# 					else:
# 						MPr_P[copy_of_plane].append(saddle_connection_endpoints[i])
# 			else:	
# 				for endpoint in saddle_connection_endpoints:
# 					MPr_P[copy_of_plane].append(endpoint)
# 		MPr.append(MPr_P)		
# 	return MPr 


# marked_periods4(radius)[P][copy_of_plane][i] outputs a list [holonomy,start_data,end_data] with a holonomy vector in the (copy_of_plane)^th copy of the plane associated to singularity P along with the start data (polygon,vertex) and end data (polygon,vertex) of the corresponding saddle connection in X (to be used for determining the $\mathbb{Z}_2$-action)
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
						MPr_P[copy_of_plane-1].append([sc.holonomy(),sc.start_data(),sc.end_data()])
					# Otherwise its holonomy vector is assigned to the current copy of the plane
					else:
						MPr_P[copy_of_plane].append([sc.holonomy(),sc.start_data(),sc.end_data()])

			# These are the unambiguous vertices
			else:	
				for i in range(len(saddle_connections)):
					sc=saddle_connections[i]
					MPr_P[copy_of_plane].append([sc.holonomy(),sc.start_data(),sc.end_data()])
		MPr.append(MPr_P)
	return MPr

# Z2_action(marked_periods)[P0][copy_of_plane0][i] gives [holonomy,(P1,copy_of_plane1)] where holonomy is the i^th holonomy vector in the (copy_of_plane0)^th copy of the plane associated to singularity P0, and (P1,copy_of_plane1) is the singularity and copy of plane where the holonomy vector of the image under the $\mathbb{Z}_2$-action of the original holonomy ends up; note that this image is -holonomy
def Z2_action(marked_periods):
	marked_periods_copy=deepcopy(marked_periods)
	for P0 in range(X._num_singularities):
		for copy_of_plane0 in range(len(marked_periods[P0])):
			for i in range(len(marked_periods[P0][copy_of_plane0])):
				sc0=marked_periods[P0][copy_of_plane0][i]
				end=False
				for P1 in range(X._num_singularities):
					for copy_of_plane1 in range(len(marked_periods[P1])):
						for k in range(len(marked_periods[P1][copy_of_plane1])):
							sc1=marked_periods[P1][copy_of_plane1][k]
							if sc0[0]==-sc1[0] and sc0[1]==sc1[2] and sc0[2]==sc1[1]:
								marked_periods_copy[P0][copy_of_plane0][i]=[sc0[0],(P1,copy_of_plane1)]
								end=True
								break
						if end==True:
							break
					if end==True:
						break
			marked_periods_copy[P0][copy_of_plane0].sort()
	return(marked_periods_copy)


def marked_periods2(radius):
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
# Perhaps the $\mathbb{Z}_2$-action data is necessary when X is the L-surface: when radius=2*r and norm_bound=chi_1_inv(2*rho/radius), the matrix [[1,1],[0,1]] is outputted from this current A_r() but is not in the Veech group of X
def A_r(radius,norm_bound):
	# matrices_to_check will be a list comprised of matrices in $\text{SL}_2\mathbb{Z}$ whose Frobenius norm is bounded above by norm_bound and which take a basis of vectors in MPr_bounded_by_2rho to another basis of vectors in MPr_bounded_by_2rho
	matrices_to_check=[]
	MPr=marked_periods(radius)
	MPr_bounded_by_2rho=[]
	preemptive_test=[]
	for P in range(X._num_singularities):
		MPr_P_bounded_by_2rho=[]
		for copy_of_plane in range(len(MPr[P])):
			v0=MPr[P][copy_of_plane][0]
			v1=MPr[P][copy_of_plane][1]
			MPr_P_copy_of_plane_bounded_by_2rho=[]
			preemptive_test_copy_of_plane=[]
			for mp0 in MPr[P][copy_of_plane]:
				if mp0[0]^2+mp0[1]^2<=(2*rho)^2:
					MPr_P_copy_of_plane_bounded_by_2rho.append(tuple(mp0))
					preemptive_test_copy_of_plane.append(tuple(mp0))
				for mp1 in MPr[P][copy_of_plane]:
					M=matrix([[mp0[0],mp1[0]],[mp0[1],mp1[1]]])
					if M.det()!=0:
						w0=M.inverse()*vector(v0)
						w1=M.inverse()*vector(v1)
						N=matrix([[w0[0],w1[0]],[w0[1],w1[1]]])
						frobenius_norm_squared=w0[0]^2+w1[0]^2+w0[1]^2+w1[1]^2
						if N not in matrices_to_check and N.det()==1 and frobenius_norm_squared<=norm_bound^2:
							matrices_to_check.append(N)
			MPr_P_bounded_by_2rho.append(set(MPr_P_copy_of_plane_bounded_by_2rho))
			preemptive_test.append(set(preemptive_test_copy_of_plane))
		MPr_bounded_by_2rho.append(MPr_P_bounded_by_2rho)

	# The matrices_such_that_points_match list will record matrices from matrices_to_check such that points in MPr_bounded_by_2rho match points in F_M_MPr_bounded_by_2rho (up to permutation of singularities and cyclic permutations of copys of planes for a particular singularity); this does not yet check the preservation of the $\mathbb{Z}_2$-action
	matrices_such_that_points_match=[]			
	for M in matrices_to_check:
		MPr_bounded_by_2rho_alterable=copy(MPr_bounded_by_2rho)
		F_M_MPr_bounded_by_2rho=[]
		end0=0
		end1=0
		end2=0

		for P in range(X._num_singularities):
			F_M_MPr_P_bounded_by_2rho=[]
			for copy_of_plane in range(len(MPr[P])):
				F_M_MPr_P_copy_of_plane_bounded_by_2rho=[]
				F_M_preemptive_test=[]
				for mp in MPr[P][copy_of_plane]:
					point=M*vector([mp[0],mp[1]])
					if point[0]^2+point[1]^2<=(2*rho)^2:
						F_M_MPr_P_copy_of_plane_bounded_by_2rho.append(tuple(point))
						F_M_preemptive_test.append(tuple(point))
				if set(F_M_preemptive_test) not in preemptive_test:
					end0=1
					break
				F_M_MPr_P_bounded_by_2rho.append(set(F_M_MPr_P_copy_of_plane_bounded_by_2rho))				
			if end0==1:
				end1=1
				break
			F_M_MPr_bounded_by_2rho.append(F_M_MPr_P_bounded_by_2rho)

		if end1!=1:
			for P in range(X._num_singularities):
				length_flag=len(MPr_bounded_by_2rho_alterable)
				for Q in range(len(MPr_bounded_by_2rho_alterable)):
					if F_M_MPr_bounded_by_2rho[P]==MPr_bounded_by_2rho_alterable[Q]:
						MPr_bounded_by_2rho_alterable.remove(MPr_bounded_by_2rho_alterable[Q])
						break
					else:
						if len(F_M_MPr_bounded_by_2rho[P])==len(MPr_bounded_by_2rho_alterable[Q]):
							for copy_of_plane in range(len(MPr_bounded_by_2rho_alterable[Q])-1):
								hold=MPr_bounded_by_2rho_alterable[Q][0]
								MPr_bounded_by_2rho_alterable[Q].remove(hold)
								MPr_bounded_by_2rho_alterable[Q].append(hold)
								if F_M_MPr_bounded_by_2rho[P]==MPr_bounded_by_2rho_alterable[Q0]:
									MPr_bounded_by_2rho_alterable.remove(MPr_bounded_by_2rho_alterable[Q0])
									break
					if len(MPr_bounded_by_2rho_alterable)==length_flag:
						end2=1
				if end2==1:
					break
			if len(MPr_bounded_by_2rho_alterable)==0:
				matrices_such_that_points_match.append(M)

	# Of the matrices in matrices_such_that_points_match, A_r will keep only those matrices that preserve the $\mathbb{Z}_2$-action on MPr_bounded_by_2rho
	# A_r=[]
	# MPr_bounded_by_2rho_Z2_action=[[[[MPr_bounded_by_2rho[P][copy_of_plane][point],a,b] for point in range(len(MPr_bounded_by_2rho[P][copy_of_plane]))] for copy_of_plane in range(len(MPr_bounded_by_2rho[P]))] for P in range(len(MPr_bounded_by_2rho))]
	return(matrices_such_that_points_match)

def A_r2(radius,norm_bound):
	# matrices_to_check will be a list comprised of matrices in $\text{SL}_2\mathbb{Z}$ whose Frobenius norm is bounded above by norm_bound and which take a basis of vectors in MPr_bounded_by_2rho to another basis of vectors in MPr_bounded_by_2rho
	matrices_to_check=[]
	MPr=marked_periods(radius)
	MPr_holonomies=[[[MPr[P][copy_of_plane][holonomy] for holonomy in range(len(MPr[P][copy_of_plane]))] for copy_of_plane in range(len(MPr[P]))] for P in range(X._num_singularities)]
	MPr_holonomies_bounded_by_2rho=[]
	preemptive_test=[]
	for P in range(X._num_singularities):
		MPr_P_holonomies_bounded_by_2rho=[]
		for copy_of_plane in range(len(MPr[P])):
			v0=MPr_holonomies[P][copy_of_plane][0]
			v1=MPr_holonomies[P][copy_of_plane][1]
			MPr_P_copy_of_plane_holonomies_bounded_by_2rho=[]
			preemptive_test_copy_of_plane=[]
			for mp0 in MPr_holonomies[P][copy_of_plane]:
				if mp0[0]^2+mp0[1]^2<=(2*rho)^2:
					MPr_P_copy_of_plane_holonomies_bounded_by_2rho.append(tuple(mp0))
					preemptive_test_copy_of_plane.append(tuple(mp0))
				for mp1 in MPr_holonomies[P][copy_of_plane]:
					M=matrix([[mp0[0],mp1[0]],[mp0[1],mp1[1]]])
					if M.det()!=0:
						w0=M.inverse()*vector(v0)
						w1=M.inverse()*vector(v1)
						N=matrix([[w0[0],w1[0]],[w0[1],w1[1]]])
						frobenius_norm_squared=w0[0]^2+w1[0]^2+w0[1]^2+w1[1]^2
						if N not in matrices_to_check and N.det()==1 and frobenius_norm_squared<=norm_bound^2:
							matrices_to_check.append(N)
			MPr_P_holonomies_bounded_by_2rho.append(set(MPr_P_copy_of_plane_holonomies_bounded_by_2rho))
			preemptive_test.append(set(preemptive_test_copy_of_plane))
		MPr_holonomies_bounded_by_2rho.append(MPr_P_bounded_by_2rho)

	# The matrices_such_that_holonomies_match list will record matrices from matrices_to_check such that points in MPr_bounded_by_2rho match points in F_M_MPr_bounded_by_2rho (up to permutation of singularities and cyclic permutations of copys of planes for a particular singularity); this does not yet check the preservation of the $\mathbb{Z}_2$-action
	matrices_such_that_holonomies_match=[]			
	for M in matrices_to_check:
		MPr_bounded_by_2rho_alterable=copy(MPr_holonomies_bounded_by_2rho)
		F_M_MPr_bounded_by_2rho=[]
		end0=0
		end1=0
		end2=0

		for P in range(X._num_singularities):
			F_M_MPr_P_bounded_by_2rho=[]
			for copy_of_plane in range(len(MPr[P])):
				F_M_MPr_P_copy_of_plane_bounded_by_2rho=[]
				F_M_preemptive_test=[]
				for mp in MPr[P][copy_of_plane]:
					point=M*vector([mp[0],mp[1]])
					if point[0]^2+point[1]^2<=(2*rho)^2:
						F_M_MPr_P_copy_of_plane_bounded_by_2rho.append(tuple(point))
						F_M_preemptive_test.append(tuple(point))
				if set(F_M_preemptive_test) not in preemptive_test:
					end0=1
					break
				F_M_MPr_P_bounded_by_2rho.append(set(F_M_MPr_P_copy_of_plane_bounded_by_2rho))				
			if end0==1:
				end1=1
				break
			F_M_MPr_bounded_by_2rho.append(F_M_MPr_P_bounded_by_2rho)

		if end1!=1:
			for P in range(X._num_singularities):
				length_flag=len(MPr_bounded_by_2rho_alterable)
				for Q in range(len(MPr_bounded_by_2rho_alterable)):
					if F_M_MPr_bounded_by_2rho[P]==MPr_bounded_by_2rho_alterable[Q]:
						MPr_bounded_by_2rho_alterable.remove(MPr_bounded_by_2rho_alterable[Q])
						break
					else:
						if len(F_M_MPr_bounded_by_2rho[P])==len(MPr_bounded_by_2rho_alterable[Q]):
							for copy_of_plane in range(len(MPr_bounded_by_2rho_alterable[Q])-1):
								hold=MPr_bounded_by_2rho_alterable[Q][0]
								MPr_bounded_by_2rho_alterable[Q].remove(hold)
								MPr_bounded_by_2rho_alterable[Q].append(hold)
								if F_M_MPr_bounded_by_2rho[P]==MPr_bounded_by_2rho_alterable[Q0]:
									MPr_bounded_by_2rho_alterable.remove(MPr_bounded_by_2rho_alterable[Q0])
									break
					if len(MPr_bounded_by_2rho_alterable)==length_flag:
						end2=1
				if end2==1:
					break
			if len(MPr_bounded_by_2rho_alterable)==0:
				matrices_such_that_holonomies_match.append(M)

	# Of the matrices in matrices_such_that_points_match, A_r will keep only those matrices that preserve the $\mathbb{Z}_2$-action on MPr_bounded_by_2rho
	A_r=[]
	


	return(matrices_such_that_holonomies_match)


def A_r3(radius,norm_bound):
	# matrices_to_check will be a list comprised of matrices in $\text{SL}_2\mathbb{Z}$ whose Frobenius norms are bounded above by norm_bound and which take a basis of vectors in marked_periods(radius) to another basis of vectors in marked_periods(radius)
	#### CHECK THIS TO MAKE SURE matrices_to_check IS EXHAUSTIVE; IT DOESN'T SEEM TO FIND [[1,0],[2,1]] FOR THE L-SHAPED SURFACE
	matrices_to_check=[]
	MPr=marked_periods(radius)
	MPr_bounded_by_2rho=[]
	# preemptive_test=[]
	for P in range(X._num_singularities):
		MPr_P_bounded_by_2rho=[]
		for copy_of_plane in range(len(MPr[P])):
			v0=MPr[P][copy_of_plane][0][0]
			v1=MPr[P][copy_of_plane][1][0]
			MPr_P_copy_of_plane_bounded_by_2rho=[]
			# preemptive_test_copy_of_plane=[]
			for i in range(len(MPr[P][copy_of_plane])):
				mp0=MPr[P][copy_of_plane][i][0]
				if mp0[0]^2+mp0[1]^2<=(2*rho)^2:
					# MPr_P_copy_of_plane_bounded_by_2rho.append(tuple(mp0))
					MPr_P_copy_of_plane_bounded_by_2rho.append(MPr[P][copy_of_plane][i])
					# preemptive_test_copy_of_plane.append(tuple(mp0))
				for k in range(len(MPr[P][copy_of_plane])):
					mp1=MPr[P][copy_of_plane][k][0]
					M=matrix([[mp0[0],mp1[0]],[mp0[1],mp1[1]]])
					if M.det()!=0:
						w0=M.inverse()*vector(v0)
						w1=M.inverse()*vector(v1)
						N=matrix([[w0[0],w1[0]],[w0[1],w1[1]]])
						frobenius_norm_squared=w0[0]^2+w1[0]^2+w0[1]^2+w1[1]^2
						if N not in matrices_to_check and N.det()==1 and frobenius_norm_squared<=norm_bound^2:
							matrices_to_check.append(N)
			MPr_P_bounded_by_2rho.append((MPr_P_copy_of_plane_bounded_by_2rho))
			# preemptive_test.append(set(preemptive_test_copy_of_plane))
		MPr_bounded_by_2rho.append(MPr_P_bounded_by_2rho)

	# The matrices_such_that_points_match list will record matrices from matrices_to_check such that points in MPr_bounded_by_2rho match points in F_M_MPr_bounded_by_2rho (up to permutation of singularities and cyclic permutations of copys of planes for a particular singularity); this does not yet check the preservation of the $\mathbb{Z}_2$-action
	matrices_such_that_points_match=[]			
	for M in matrices_to_check:
		# MPr_bounded_by_2rho_alterable=deepcopy(MPr_bounded_by_2rho)
		F_M_MPr_bounded_by_2rho=[]
		end0=0
		end1=0
		end2=0

		for P in range(X._num_singularities):
			F_M_MPr_P_bounded_by_2rho=[]
			for copy_of_plane in range(len(MPr[P])):
				F_M_MPr_P_copy_of_plane_bounded_by_2rho=[]
				# F_M_preemptive_test=[]
				for i in range(len(MPr[P][copy_of_plane])):
					mp=MPr[P][copy_of_plane][i][0]
					point=M*vector([mp[0],mp[1]])
					if point[0]^2+point[1]^2<=(2*rho)^2:
						# F_M_MPr_P_copy_of_plane_bounded_by_2rho.append(tuple(point))
						F_M_MPr_P_copy_of_plane_bounded_by_2rho.append([point,MPr[P][copy_of_plane][i][1],MPr[P][copy_of_plane][i][2]])
						# F_M_preemptive_test.append(tuple(point))
				# if set(F_M_preemptive_test) not in preemptive_test:
					# end0=1
					# break
				# F_M_MPr_P_bounded_by_2rho.append(set(F_M_MPr_P_copy_of_plane_bounded_by_2rho))
				F_M_MPr_P_bounded_by_2rho.append(F_M_MPr_P_copy_of_plane_bounded_by_2rho)				
			# if end0==1:
				# end1=1
				# break
			F_M_MPr_bounded_by_2rho.append(F_M_MPr_P_bounded_by_2rho)
### Apply Z2-action to both MPr_bounded_by_2rho and F_M_MPr_bounded_by_2rho; see if same; may have to permute (P,copy_of_plane) data to check if the same ###

		A_r=[]
		MPr_bounded_by_2rho_Z2_action=Z2_action(MPr_bounded_by_2rho)
		F_M_MPr_bounded_by_2rho_Z2_action=Z2_action(F_M_MPr_bounded_by_2rho)
		# if MPr_bounded_by_2rho_Z2_action==F_M_MPr_bounded_by_2rho_Z2_action:
		# 	A_r.append(M)
		# else:

		# Sorts lists by size of cone angles
		MPr_bounded_by_2rho_Z2_action_sorted=sorted(MPr_bounded_by_2rho_Z2_action,key=len)
		F_M_MPr_bounded_by_2rho_Z2_action_sorted=sorted(F_M_MPr_bounded_by_2rho_Z2_action,key=len)
		# List of sizes of cone angles after sorted
		cone_angles=[len(MPr_bounded_by_2rho_Z2_action_sorted[P]) for P in range(X._num_singularities)]
		# Records number of singularities of a particular order; will create symmetric groups on these number of elements and later take cartesian product of these groups for various permutations singularities
		permutation_group_orders=[]
		angle_previous=None
		for i in range(len(cone_angles)):
			angle=cone_angles[i]
			if angle!=angle_previous:
				permutation_group_orders.append(cone_angles.count(angle))
			angle_previous=angle
		# Product of symmetric groups of orders from permutation_group_orders
		permutations_of_singularities=SymmetricGroup(permutation_group_orders[0])
		for i in permutation_group_orders[1:]:
			permutations_of_singularities=permutations_of_singularities.cartesian_product(SymmetricGroup(i))			
		for permutation_sing in permutations_of_singularities.list():
			place=0
			for n in range(len(permutation_group_orders)):
				num_sings_of_same_angle=permutation_group_orders[n]
				F_M_MPr_bounded_by_2rho_Z2_action_sorted[place:place+num_sings_of_same_angle]=permutation_sing(F_M_MPr_bounded_by_2rho_Z2_action_sorted[place:place+num_sings_of_same_angle])
				place+=n 
			### We've permuted singularities; first check if holonomies match and then if the sorted sets match up to permutation of singularities in Z2-action data and cyclic permutation of copy of planes in Z2-action data





		# MPr_cone_angles_list=[len(MPr_bounded_by_2rho_Z2_action[P]) for P in range(X._num_singularities)]
		# F_M_MPr_cone_angles_list=[len(F_M_MPr_bounded_by_2rho_Z2_action[P]) for P in range(X._num_singularities)]



		# if end1!=1:
		# 	for P0 in range(X._num_singularities):
		# 		length_flag=len(MPr_bounded_by_2rho_Z2_action)
		# 		for P1 in range(len(MPr_bounded_by_2rho_Z2_action)):
		# 			if F_M_MPr_bounded_by_2rho_Z2_action[P0]==MPr_bounded_by_2rho_Z2_action[P1]:
		# 				MPr_bounded_by_2rho_Z2_action.remove(MPr_bounded_by_2rho_Z2_action[P1])
		# 				break
		# 			elif len(F_M_MPr_bounded_by_2rho_Z2_action[P0])==len(MPr_bounded_by_2rho_Z2_action[P1]):
		# 				for copy_of_plane in range(len(MPr_bounded_by_2rho_Z2_action[P1])-1):
		# 					hold=MPr_bounded_by_2rho_Z2_action[P1][0]
		# 					MPr_bounded_by_2rho_Z2_action[P1].remove(hold)
		# 					MPr_bounded_by_2rho_Z2_action[P1].append(hold)
		# 					num_F_M_holonomies=len(F_M_MPr_bounded_by_2rho_Z2_action[P0][copy_of_plane])
		# 					if len(MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane])==num_F_M_holonomies:
		# 						for i in range(num_F_M_holonomies):
		# 							if [F_M_MPr_bounded_by_2rho_Z2_action[P0][copy_of_plane][sc][0] for sc in range(len(F_M_MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane]))]==[MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane][(sc+i)%num_F_M_holonomies][0] for sc in range(len(MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane]))]:
		# 								if [F_M_MPr_bounded_by_2rho_Z2_action[P0][copy_of_plane][sc][1:] for sc in range(len(F_M_MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane]))]==[MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane][sc+i][1:] for sc in range(len(MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane]))]:
		# 									MPr_bounded_by_2rho_Z2_action.remove(MPr_bounded_by_2rho_Z2_action[Q0])
		# 									break
		# 								# else:
		# 									# for P in range(X._num_singularities-1)
		# 			if len(MPr_bounded_by_2rho_Z2_action)==length_flag:
		# 				end2=1
		# 		if end2==1:
		# 			break
		# 	if len(MPr_bounded_by_2rho_Z2_action)==0:
		# 		matrices_such_that_points_match.append(M)

	### As a first check, just keep track of various matrices_such_that_points_match; for each keep track of corresponding F_M_MPr_bounded_by_2rho_Z2_action, and afterwards check if Z2-action mathces MPr_bounded_by_2rho_Z2_action
	

		# if end1!=1:
		# 	for P in range(X._num_singularities):
		# 		length_flag=len(MPr_bounded_by_2rho_alterable)
		# 		for Q in range(len(MPr_bounded_by_2rho_alterable)):
		# 			if F_M_MPr_bounded_by_2rho[P]==MPr_bounded_by_2rho_alterable[Q]:
		# 				MPr_bounded_by_2rho_alterable.remove(MPr_bounded_by_2rho_alterable[Q])
		# 				break
		# 			elif len(F_M_MPr_bounded_by_2rho[P])==len(MPr_bounded_by_2rho_alterable[Q]):
		# 				for copy_of_plane in range(len(MPr_bounded_by_2rho_alterable[Q])-1):
		# 					hold=MPr_bounded_by_2rho_alterable[Q][0]
		# 					MPr_bounded_by_2rho_alterable[Q].remove(hold)
		# 					MPr_bounded_by_2rho_alterable[Q].append(hold)
		# 					if F_M_MPr_bounded_by_2rho[P]==MPr_bounded_by_2rho_alterable[Q0]:
		# 						MPr_bounded_by_2rho_alterable.remove(MPr_bounded_by_2rho_alterable[Q0])
		# 						break
		# 			if len(MPr_bounded_by_2rho_alterable)==length_flag:
		# 				end2=1
		# 		if end2==1:
		# 			break
		# 	if len(MPr_bounded_by_2rho_alterable)==0:
		# 		matrices_such_that_points_match.append(M)
	# Of the matrices in matrices_such_that_points_match, A_r will keep only those matrices that preserve the $\mathbb{Z}_2$-action on MPr_bounded_by_2rho
	# A_r=[]
	# MPr_bounded_by_2rho_Z2_action=[[[[MPr_bounded_by_2rho[P][copy_of_plane][point],a,b] for point in range(len(MPr_bounded_by_2rho[P][copy_of_plane]))] for copy_of_plane in range(len(MPr_bounded_by_2rho[P]))] for P in range(len(MPr_bounded_by_2rho))]
	return(matrices_such_that_points_match)


def A_r4(radius,norm_bound):
	# matrices_to_check will be a list comprised of matrices in $\text{SL}_2\mathbb{Z}$ whose Frobenius norms are bounded above by norm_bound and which take a basis of vectors in marked_periods(radius) to another basis of vectors in marked_periods(radius)
	#### CHECK THIS TO MAKE SURE matrices_to_check IS EXHAUSTIVE; IT DOESN'T SEEM TO FIND [[1,0],[2,1]] FOR THE L-SHAPED SURFACE
	matrices_to_check=[]
	MPr=marked_periods2(radius)
	MPr_bounded_by_2rho=[]
	# preemptive_test=[]
	for P in range(X._num_singularities):
		MPr_P_bounded_by_2rho=[]
		for copy_of_plane in range(len(MPr[P])):
			v0=MPr[P][copy_of_plane][0][0]
			v1=MPr[P][copy_of_plane][1][0]
			MPr_P_copy_of_plane_bounded_by_2rho=[]
			# preemptive_test_copy_of_plane=[]
			for i in range(len(MPr[P][copy_of_plane])):
				mp0=MPr[P][copy_of_plane][i][0]
				if mp0[0]^2+mp0[1]^2<=(2*rho)^2:
					# MPr_P_copy_of_plane_bounded_by_2rho.append(tuple(mp0))
					MPr_P_copy_of_plane_bounded_by_2rho.append(MPr[P][copy_of_plane][i])
					# preemptive_test_copy_of_plane.append(tuple(mp0))
				for k in range(len(MPr[P][copy_of_plane])):
					mp1=MPr[P][copy_of_plane][k][0]
					M=matrix([[mp0[0],mp1[0]],[mp0[1],mp1[1]]])
					if M.det()!=0:
						w0=M.inverse()*vector(v0)
						w1=M.inverse()*vector(v1)
						N=matrix([[w0[0],w1[0]],[w0[1],w1[1]]])
						frobenius_norm_squared=w0[0]^2+w1[0]^2+w0[1]^2+w1[1]^2
						if N not in matrices_to_check and N.det()==1 and frobenius_norm_squared<=norm_bound^2:
							matrices_to_check.append(N)
			MPr_P_bounded_by_2rho.append(sorted((MPr_P_copy_of_plane_bounded_by_2rho)))
			# preemptive_test.append(set(preemptive_test_copy_of_plane))
		MPr_bounded_by_2rho.append(MPr_P_bounded_by_2rho)
		MPr_bounded_by_2rho_sorted=sorted(MPr_bounded_by_2rho,key=len)

	# The matrices_such_that_points_match list will record matrices from matrices_to_check such that points in MPr_bounded_by_2rho match points in F_M_MPr_bounded_by_2rho (up to permutation of singularities and cyclic permutations of copys of planes for a particular singularity); this does not yet check the preservation of the $\mathbb{Z}_2$-action
	A_r=[]		
	for M in matrices_to_check:
		# MPr_bounded_by_2rho_alterable=deepcopy(MPr_bounded_by_2rho)
		F_M_MPr_bounded_by_2rho=[]
		end0=0
		end1=0
		end2=0

		for P in range(X._num_singularities):
			F_M_MPr_P_bounded_by_2rho=[]
			for copy_of_plane in range(len(MPr[P])):
				F_M_MPr_P_copy_of_plane_bounded_by_2rho=[]
				# F_M_preemptive_test=[]
				for i in range(len(MPr[P][copy_of_plane])):
					mp=MPr[P][copy_of_plane][i][0]
					point=M*vector([mp[0],mp[1]])
					if point[0]^2+point[1]^2<=(2*rho)^2:
						# F_M_MPr_P_copy_of_plane_bounded_by_2rho.append(tuple(point))
						F_M_MPr_P_copy_of_plane_bounded_by_2rho.append([point,MPr[P][copy_of_plane][i][1]])
						# F_M_preemptive_test.append(tuple(point))
				# if set(F_M_preemptive_test) not in preemptive_test:
					# end0=1
					# break
				# F_M_MPr_P_bounded_by_2rho.append(set(F_M_MPr_P_copy_of_plane_bounded_by_2rho))
				F_M_MPr_P_bounded_by_2rho.append(sorted(F_M_MPr_P_copy_of_plane_bounded_by_2rho))
			# if end0==1:
				# end1=1
				# break
			F_M_MPr_bounded_by_2rho.append(F_M_MPr_P_bounded_by_2rho)
			F_M_MPr_bounded_by_2rho_sorted=sorted(F_M_MPr_bounded_by_2rho)

### DO PERMUTATIONS AND SEE IF HOLONOMIES MATCH; IF SO, TEST Z2-ACTION ###

		# List of sizes of cone angles after sorted
		cone_angles=[len(MPr_bounded_by_2rho_sorted[P]) for P in range(X._num_singularities)]
		# Records number of singularities of a particular order; will create symmetric groups on these number of elements and later take cartesian product of these groups for various permutations of singularities
		permutation_group_orders=[]
		angle_previous=None
		for i in range(len(cone_angles)):
			angle=cone_angles[i]
			if angle!=angle_previous:
				permutation_group_orders.append(cone_angles.count(angle))
			angle_previous=angle

		# Product of symmetric groups of orders from permutation_group_orders; will be used to permute singularities of the same order
		permutations_of_singularities=SymmetricGroup(permutation_group_orders[0])
		for i in permutation_group_orders[1:]:
			permutations_of_singularities=permutations_of_singularities.cartesian_product(SymmetricGroup(i))

		# Product of cyclic groups from orders of cone angles; will be used to check permutations of copies of the plane associated to a particular singularity
		permutations_of_copies_of_plane=CyclicPermutationGroup(cone_angles[0])
		for j in cone_angles[1:]:
			permutations_of_copies_of_plane=permutations_of_copies_of_plane.cartesian_product(CyclicPermutationGroup(j))

		# Permute singularities and copies of planes of F_M_MPr_bounded_by_2rho_sorted and see if the sets of holomonies (together with a singularity of specific order and ordering of copies of planes) matches the set of holonomies of MPr_bounded_by_2rho_sorted
		for permute_singularites in permutations_of_singularities.list():
			F_M_copy=deepcopy(F_M_MPr_bounded_by_2rho_sorted)
			place=0
			for n in range(len(permutation_group_orders)):
				num_sings_of_same_angle=permutation_group_orders[n]
				F_M_copy[place:place+num_sings_of_same_angle]=permute_singularites(F_M_copy[place:place+num_sings_of_same_angle])
				place+=n 
			### SOMETHING DOESN'T SEEM TO BE WORKING HERE
			for permute_planes in permutations_of_copies_of_plane.list():
				if X._num_singularities==1:
					F_M_copy[0]=permute_planes(F_M_copy[0])
				else:
					for P in range(X._num_singularities):
						permute_planes_P=permute_planes[P]
						F_M_copy[P]=permute_planes_P(F_M_copy[P])

			### CHECK IF HOLONOMIES MATCH





		# MPr_cone_angles_list=[len(MPr_bounded_by_2rho_Z2_action[P]) for P in range(X._num_singularities)]
		# F_M_MPr_cone_angles_list=[len(F_M_MPr_bounded_by_2rho_Z2_action[P]) for P in range(X._num_singularities)]



		# if end1!=1:
		# 	for P0 in range(X._num_singularities):
		# 		length_flag=len(MPr_bounded_by_2rho_Z2_action)
		# 		for P1 in range(len(MPr_bounded_by_2rho_Z2_action)):
		# 			if F_M_MPr_bounded_by_2rho_Z2_action[P0]==MPr_bounded_by_2rho_Z2_action[P1]:
		# 				MPr_bounded_by_2rho_Z2_action.remove(MPr_bounded_by_2rho_Z2_action[P1])
		# 				break
		# 			elif len(F_M_MPr_bounded_by_2rho_Z2_action[P0])==len(MPr_bounded_by_2rho_Z2_action[P1]):
		# 				for copy_of_plane in range(len(MPr_bounded_by_2rho_Z2_action[P1])-1):
		# 					hold=MPr_bounded_by_2rho_Z2_action[P1][0]
		# 					MPr_bounded_by_2rho_Z2_action[P1].remove(hold)
		# 					MPr_bounded_by_2rho_Z2_action[P1].append(hold)
		# 					num_F_M_holonomies=len(F_M_MPr_bounded_by_2rho_Z2_action[P0][copy_of_plane])
		# 					if len(MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane])==num_F_M_holonomies:
		# 						for i in range(num_F_M_holonomies):
		# 							if [F_M_MPr_bounded_by_2rho_Z2_action[P0][copy_of_plane][sc][0] for sc in range(len(F_M_MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane]))]==[MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane][(sc+i)%num_F_M_holonomies][0] for sc in range(len(MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane]))]:
		# 								if [F_M_MPr_bounded_by_2rho_Z2_action[P0][copy_of_plane][sc][1:] for sc in range(len(F_M_MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane]))]==[MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane][sc+i][1:] for sc in range(len(MPr_bounded_by_2rho_Z2_action[P1][copy_of_plane]))]:
		# 									MPr_bounded_by_2rho_Z2_action.remove(MPr_bounded_by_2rho_Z2_action[Q0])
		# 									break
		# 								# else:
		# 									# for P in range(X._num_singularities-1)
		# 			if len(MPr_bounded_by_2rho_Z2_action)==length_flag:
		# 				end2=1
		# 		if end2==1:
		# 			break
		# 	if len(MPr_bounded_by_2rho_Z2_action)==0:
		# 		matrices_such_that_points_match.append(M)

	### As a first check, just keep track of various matrices_such_that_points_match; for each keep track of corresponding F_M_MPr_bounded_by_2rho_Z2_action, and afterwards check if Z2-action mathces MPr_bounded_by_2rho_Z2_action
	

		# if end1!=1:
		# 	for P in range(X._num_singularities):
		# 		length_flag=len(MPr_bounded_by_2rho_alterable)
		# 		for Q in range(len(MPr_bounded_by_2rho_alterable)):
		# 			if F_M_MPr_bounded_by_2rho[P]==MPr_bounded_by_2rho_alterable[Q]:
		# 				MPr_bounded_by_2rho_alterable.remove(MPr_bounded_by_2rho_alterable[Q])
		# 				break
		# 			elif len(F_M_MPr_bounded_by_2rho[P])==len(MPr_bounded_by_2rho_alterable[Q]):
		# 				for copy_of_plane in range(len(MPr_bounded_by_2rho_alterable[Q])-1):
		# 					hold=MPr_bounded_by_2rho_alterable[Q][0]
		# 					MPr_bounded_by_2rho_alterable[Q].remove(hold)
		# 					MPr_bounded_by_2rho_alterable[Q].append(hold)
		# 					if F_M_MPr_bounded_by_2rho[P]==MPr_bounded_by_2rho_alterable[Q0]:
		# 						MPr_bounded_by_2rho_alterable.remove(MPr_bounded_by_2rho_alterable[Q0])
		# 						break
		# 			if len(MPr_bounded_by_2rho_alterable)==length_flag:
		# 				end2=1
		# 		if end2==1:
		# 			break
		# 	if len(MPr_bounded_by_2rho_alterable)==0:
		# 		matrices_such_that_points_match.append(M)
	# Of the matrices in matrices_such_that_points_match, A_r will keep only those matrices that preserve the $\mathbb{Z}_2$-action on MPr_bounded_by_2rho
	# A_r=[]
	# MPr_bounded_by_2rho_Z2_action=[[[[MPr_bounded_by_2rho[P][copy_of_plane][point],a,b] for point in range(len(MPr_bounded_by_2rho[P][copy_of_plane]))] for copy_of_plane in range(len(MPr_bounded_by_2rho[P]))] for P in range(len(MPr_bounded_by_2rho))]
	return(matrices_such_that_points_match)
