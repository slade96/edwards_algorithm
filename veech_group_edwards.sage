# ALGORITHM 7.1

# Takes as input translation surface (X,omega) such that $\Gamma(X,omega)\cap \text{SO}(2,\mathbb{R})\subseteq \{\pm \text{Id}\}$

def veech_group(X):
	from flatsurf.geometry.polygon import polygons,PolygonPosition
	assert isinstance(X,TranslationSurface)


#0 Compute the Voronoi decomposition of X
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
				if already_accounted==0:
					vertex_equivalence_classes_sets.append(vertex_is_equivalent_to)
					vertex_equivalence_classes.append([v for v in vertex_is_equivalent_to])

	# Gives the vertices of each polygon of X in AA so that the Voronoi Diagram can be computed with these points
	vertices=[[[AA(c) for c in vert] for vert in X.polygon(poly).vertices()] for poly in range(X._num_polygons)]
	# Gives a list whose entries are the Voronoi Diagrams subordinate to the vertices of each polygon
	VD=[VoronoiDiagram(vertices[poly]) for poly in range(X._num_polygons)]
	# A list of lists; each entry of the outer list corresponds to a polygon of X; each entry of the inner lists is the (potentially) unbounded polyhedron in the Voronoi Diagram associated to a particular vertex of said polygon
	VD_polyhedra_unbdd=[[VD[poly].regions()[VD[poly].points()[v]] for v in range(X.polygon(poly).num_edges())] for poly in range(X._num_polygons)]
	# Similar to the list above, but now the polyhedra are intersected with the polygon of interest, so they are bounded
	VD_polyhedra=[[VD_polyhedra_unbdd[poly][vert].intersection(Polyhedron(vertices=vertices[poly])) for vert in range(X.polygon(poly).num_edges())] for poly in range(X._num_polygons)]



#1 Calculate $\rho=\rho(X,\omega)$ using the Voronoi 2-cells
	# Creates a list of the max distinces in each polygon of the Voronoi Diagram from the singularity to the boundary of the polygon
	distances_from_vertices=[]
	for poly in range(X._num_polygons):
		for vert in range(X.polygon(poly).num_edges()):
			# Gives the polygon of the Voronoi Diagram associated to polygon 'poly' and vertex 'vert' of X
			Voro_poly=VD_polyhedra[poly][vert]
			# List of vertices of aforementioned polygon
			Voro_poly_verts=Voro_poly.vertices_list()
			Voro_poly_num_edges=len(Voro_poly_verts)

			# Gives main vertex, 'vert,' of X that this section is based on
			vertices_tuple=[tuple(v) for v in vertices[poly]]
			Voro_poly_verts_tuple=[tuple(v) for v in Voro_poly_verts]
			main_vertex=[list(set(vertices_tuple).intersection(Voro_poly_verts_tuple))[0][0],list(set(vertices_tuple).intersection(Voro_poly_verts_tuple))[0][1]]

			distances_from_main_vertex=[sqrt((Voro_poly_verts[v][0]-main_vertex[0])^2+(Voro_poly_verts[v][1]-main_vertex[1])^2) for v in range(Voro_poly_num_edges)]
			distances_from_vertices.append(max(distances_from_main_vertex))
	rho=max(distances_from_vertices)



#2 Let $r=2\rho$ and let $b=\chi_1^{-1}(2\rho/r)=\sqrt{2}$
	def chi_1_inv(x):
		return sqrt(x^2+x^(-2))
	r=2*rho
	b=sqrt(2)

#3 Calculate $MP_P^r(X,\omega)$ for all $P\in\Sigma$
### SEEMS TO WORK IF THE SINGULARITY IS REMOVABLE; NEED TO CONSIDER SINGULARITIES WITH ANGLE $>2\pi$ ###	


	# This records information to determine which disjoint copy of the plane an element of $MP_P^r$ should lie in
	# vertices_in_disjoint_planes[P][n] gives a list of vertex data such that separatrices leaving these vertices should develop in the n^th copy of the plane corresponding to singularity P
	# There is some ambiguity for some vertices, such as vertex 5 of the regular hexagon...these ambiguous vertices where separatrices may develop in either the n^th or (n+1)^st plane lead the list corresponding to the (n+1)^st plane
	vertices_in_disjoint_planes=[]
	for P in range(X._num_singularities):
		copy_of_plane=0
		vertices_in_disjoint_planes_P=[]
		singularity_P=sorted(vertex_equivalence_classes[P])
		vertex_data=singularity_P[0]
		P_angle=X.polygon(vertex_data[0]).angle(vertex_data[1])
		# direction=X.polygon(vertex_data[0]).edge(vertex_data[1])
		vertices_in_disjoint_planes_P_copy_of_plane=[vertex_data]
		# P_angle+=X.polygon(vertex_data[0]).angle(vertex_data[1])

		# for i in range(X.angles()[P]):
		while copy_of_plane<X.angles()[P]:
			vertex_data=X.opposite_edge(vertex_data[0],(vertex_data[1]-1)%X.polygon(vertex_data[0]).num_edges())
			P_angle+=X.polygon(vertex_data[0]).angle(vertex_data[1])
			if P_angle>copy_of_plane+1:
				vertices_in_disjoint_planes_P.append(vertices_in_disjoint_planes_P_copy_of_plane)
				vertices_in_disjoint_planes_P_copy_of_plane=[]
				copy_of_plane+=1
			vertices_in_disjoint_planes_P_copy_of_plane.append(vertex_data)

		vertices_in_disjoint_planes.append(vertices_in_disjoint_planes_P)



	### ALSO NEED TO KEEP TRACK OF $\mathbb{Z}_2$-ACTION ###
	### NEED TO CONSIDER WHEN A VERTEX'S MEMBERSHIP OF THE n^th COPY OF THE PLANE IS AMBIGUOUS ###
	# marked_periods(radius)[P][n] is the set of marked periods from singularity P of X bounded by radius corresponding to the nth copy of the plane
	def marked_periods(radius):
		MPr=[]
		for P in range(X._num_singularities):
			MPr_P=[[] for i in X.angles()[P]]
			for vertex_data in vertex_equivalence_classes[P]:
				copy_of_plane=[n for n in range(len(vertices_in_disjoint_planes[P])) if vertex_data in vertices_in_disjoint_planes[P][n]][0]
				polygon=vertex_data[0]
				vertex=vertex_data[1]
				saddle_connections=X.saddle_connections(radius^2,polygon,vertex)
				saddle_connections_directions=[sc.direction() for sc in saddle_connections]
				poly=X.polygon(polygon)
				point=poly.vertex(vertex)
				V=poly.parent().vector_space()
				saddle_connection_endpoints=[]
				for direction in saddle_connections_directions:
					is_vertex=0
					pnt=point
					polygon2=polygon
					poly2=poly
					holonomy=(0,0)
					while is_vertex==0:
						flow_to_bdry=poly2.flow_to_exit(V(pnt),V(direction))
						holonomy=V(holonomy)+V(flow_to_bdry[0])-V(pnt)
						if flow_to_bdry[1].is_vertex()==True:
							is_vertex=1
							saddle_connection_endpoints.append(holonomy)
						else:
							continue_flow_from=X.opposite_edge(polygon2,flow_to_bdry[1].get_edge())
							pnt=X.edge_transformation(polygon2,flow_to_bdry[1].get_edge())(flow_to_bdry[0])
							polygon2=continue_flow_from[0]
							poly2=X.polygon(polygon2)

				# saddle_connection_endpoints=[poly.flow_to_exit(V(point),V(direction))[0] for direction in saddle_connections_directions]
				# holonomies=[endpoint-point for endpoint in saddle_connection_endpoints]
				for endpoint in saddle_connection_endpoints:
					MPr_P[copy_of_plane].append(endpoint)
			MPr.append(MPr_P)		
		return MPr 
	



#4 Calculate $A_r=\{M\in\Gamma(X,\omega)\ |\ ||M||\le b\}=\text{SO}(2,\mathbb{R})\cap\Gamma(X,\omega)$ using Theorem 18
### LOOK OVER THIS ###

	def A_r(radius):
		# filler=[]
		matrices_to_check=[]
		MPr=marked_periods(radius)
		MPr_bounded_by_2rho=[]
		# MPr_bounded_by_2rho=set()
		for P in range(X._num_singularities):
			v0=MPr[P][0]
			for mp in MPr[P]:
				if v0[0]*mp[1]-mp[0]*v0[1]!=0:
					v1=mp
					break
			MPr_P_bounded_by_2rho=[]
			for mp0 in MPr[P]:
				if mp0[0]^2+mp0[1]^2<=(2*rho)^2:
					MPr_P_bounded_by_2rho.append(tuple(mp0))
				for mp1 in MPr[P]:
					M=matrix([[mp0[0],mp1[0]],[mp0[1],mp1[1]]])
					if M.det()!=0:
						w0=M.inverse()*vector(v0)
						w1=M.inverse()*vector(v1)
						N=matrix([[w0[0],w1[0]],[w0[1],w1[1]]])
						# filler.append(N)
						frobenius_norm_squared=w0[0]^2+w1[0]^2+w0[1]^2+w1[1]^2
						if N not in matrices_to_check and N.det()==1 and frobenius_norm_squared<=b^2:
							matrices_to_check.append(N)
			MPr_bounded_by_2rho.append(set(MPr_P_bounded_by_2rho))

		A_r=[]			
		for M in matrices_to_check:
			F_M_MPr_bounded_by_2rho=[]
			# F_M_MPr_bounded_by_2rho=set()
			for P in range(X._num_singularities):
				F_M_MPr_P_bounded_by_2rho=[]
				for mp in MPr[P]:
					point=M*vector([mp[0],mp[1]])
					if point[0]^2+point[1]^2<=(2*rho)^2:
						F_M_MPr_P_bounded_by_2rho.append(tuple(point))
				F_M_MPr_bounded_by_2rho.append(set(F_M_MPr_P_bounded_by_2rho))
				# F_M_MPr_bounded_by_2rho.add(set(F_M_MPr_P_bounded_by_2rho))
			if F_M_MPr_bounded_by_2rho==MPr_bounded_by_2rho:
				A_r.append(M)
		return(A_r)
			# [M*vector([mp[0],mp[1]]) for mp in MPr[0]]



	A_r(r)




		








