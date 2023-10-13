
#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it cannot perform an operation (i.e., because
    the resulting mesh does not have a valid representation).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
 * add_face: add a standalone face to the mesh
 *  sides: number of sides
 *  radius: distance from vertices to origin
 *
 * We provide this method as an example of how to make new halfedge mesh geometry.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::add_face(uint32_t sides, float radius) {
	//faces with fewer than three sides are invalid, so abort the operation:
	if (sides < 3) return std::nullopt;


	std::vector< VertexRef > face_vertices;
	//In order to make the first edge point in the +x direction, first vertex should
	// be at -90.0f - 0.5f * 360.0f / float(sides) degrees, so:
	float const start_angle = (-0.25f - 0.5f / float(sides)) * 2.0f * PI_F;
	for (uint32_t s = 0; s < sides; ++s) {
		float angle = float(s) / float(sides) * 2.0f * PI_F + start_angle;
		VertexRef v = emplace_vertex();
		v->position = radius * Vec3(std::cos(angle), std::sin(angle), 0.0f);
		face_vertices.emplace_back(v);
	}

	assert(face_vertices.size() == sides);

	//assemble the rest of the mesh parts:
	FaceRef face = emplace_face(false); //the face to return
	FaceRef boundary = emplace_face(true); //the boundary loop around the face

	std::vector< HalfedgeRef > face_halfedges; //will use later to set ->next pointers

	for (uint32_t s = 0; s < sides; ++s) {
		//will create elements for edge from a->b:
		VertexRef a = face_vertices[s];
		VertexRef b = face_vertices[(s+1)%sides];

		//h is the edge on face:
		HalfedgeRef h = emplace_halfedge();
		//t is the twin, lies on boundary:
		HalfedgeRef t = emplace_halfedge();
		//e is the edge corresponding to h,t:
		EdgeRef e = emplace_edge(false); //false: non-sharp

		//set element data to something reasonable:
		//(most ops will do this with interpolate_data(), but no data to interpolate here)
		h->corner_uv = a->position.xy() / (2.0f * radius) + 0.5f;
		h->corner_normal = Vec3(0.0f, 0.0f, 1.0f);
		t->corner_uv = b->position.xy() / (2.0f * radius) + 0.5f;
		t->corner_normal = Vec3(0.0f, 0.0f,-1.0f);

		//thing -> halfedge pointers:
		e->halfedge = h;
		a->halfedge = h;
		if (s == 0) face->halfedge = h;
		if (s + 1 == sides) boundary->halfedge = t;

		//halfedge -> thing pointers (except 'next' -- will set that later)
		h->twin = t;
		h->vertex = a;
		h->edge = e;
		h->face = face;

		t->twin = h;
		t->vertex = b;
		t->edge = e;
		t->face = boundary;

		face_halfedges.emplace_back(h);
	}

	assert(face_halfedges.size() == sides);

	for (uint32_t s = 0; s < sides; ++s) {
		face_halfedges[s]->next = face_halfedges[(s+1)%sides];
		face_halfedges[(s+1)%sides]->twin->next = face_halfedges[s]->twin;
	}

	return face;
}


/*
 * bisect_edge: split an edge without splitting the adjacent faces
 *  e: edge to split
 *
 * returns: added vertex
 *
 * We provide this as an example for how to implement local operations.
 * (and as a useful subroutine!)
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {
	// Phase 0: draw a picture
	//
	// before:
	//    ----h--->
	// v1 ----e--- v2
	//   <----t---
	//
	// after:
	//    --h->    --h2->
	// v1 --e-- vm --e2-- v2
	//    <-t2-    <--t--
	//

	// Phase 1: collect existing elements
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->vertex;
	VertexRef v2 = t->vertex;

	// Phase 2: Allocate new elements, set data
	VertexRef vm = emplace_vertex();
	vm->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, vm); //set bone_weights

	EdgeRef e2 = emplace_edge();
	e2->sharp = e->sharp; //copy sharpness flag

	HalfedgeRef h2 = emplace_halfedge();
	interpolate_data({h, h->next}, h2); //set corner_uv, corner_normal

	HalfedgeRef t2 = emplace_halfedge();
	interpolate_data({t, t->next}, t2); //set corner_uv, corner_normal

	// The following elements aren't necessary for the bisect_edge, but they are here to demonstrate phase 4
    FaceRef f_not_used = emplace_face();
    HalfedgeRef h_not_used = emplace_halfedge();

	// Phase 3: Reassign connectivity (careful about ordering so you don't overwrite values you may need later!)

	vm->halfedge = h2;

	e2->halfedge = h2;

	assert(e->halfedge == h); //unchanged

	//n.b. h remains on the same face so even if h->face->halfedge == h, no fixup needed (t, similarly)

	h2->twin = t;
	h2->next = h->next;
	h2->vertex = vm;
	h2->edge = e2;
	h2->face = h->face;

	t2->twin = h;
	t2->next = t->next;
	t2->vertex = vm;
	t2->edge = e;
	t2->face = t->face;
	
	h->twin = t2;
	h->next = h2;
	assert(h->vertex == v1); // unchanged
	assert(h->edge == e); // unchanged
	//h->face unchanged

	t->twin = h2;
	t->next = t2;
	assert(t->vertex == v2); // unchanged
	t->edge = e2;
	//t->face unchanged


	// Phase 4: Delete unused elements
    erase_face(f_not_used);
    erase_halfedge(h_not_used);

	// Phase 5: Return the correct iterator
	return vm;
}


/*
 * split_edge: split an edge and adjacent (non-boundary) faces
 *  e: edge to split
 *
 * returns: added vertex. vertex->halfedge should lie along e
 *
 * Note that when splitting the adjacent faces, the new edge
 * should connect to the vertex ccw from the ccw-most end of e
 * within the face.
 *
 * Do not split adjacent boundary faces.
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(EdgeRef e) {
	// A2L2 (REQUIRED): split_edge
	
	//(void)e; //this line avoids 'unused parameter' warnings. You can delete it as you fill in the function.

	HalfedgeRef he = e->halfedge;
	HalfedgeRef twin = he->twin;
	bool he_boundary = he->face->boundary;
	bool twin_boundary = twin->face->boundary;
	if(he_boundary && twin_boundary) return std::nullopt;

	//3 cases total (2 single boundary face cases, 1 no boundary face case)
	
	if(he_boundary){ //twin half-edge half-subidivision case

		//Reduce to next case
		e->halfedge = twin;
		return split_edge(e);

	}else if(twin_boundary){ //he half-edge half-subidivision case

		//Phase 1: Collect all needed elements

		//Halfedges
		HalfedgeRef he_n = he->next;
		HalfedgeRef he_n2 = he_n->next;
		HalfedgeRef he_p = he;
		do{
			he_p = he_p->next;
		}while(he_p->next != he);

		//Vertices
		VertexRef old_v_left = he_n2->vertex;
		VertexRef old_v_right = twin->vertex;
		VertexRef old_v_bot_right = he->vertex;


		//Faces
		FaceRef old_f = he->face;
		FaceRef bound_f = twin->face;

		//Phase 2: Create any new elements

		//Halfedges
		HalfedgeRef new_out = emplace_halfedge();
		HalfedgeRef new_in = emplace_halfedge();
		HalfedgeRef split_in = emplace_halfedge();
		HalfedgeRef split_out = emplace_halfedge();

		//Edges
		EdgeRef new_e = emplace_edge();
		EdgeRef split_e = emplace_edge();

		//Vertices
		VertexRef new_v = emplace_vertex();

		//Faces
		FaceRef new_f = emplace_face();

		//Set data for new things

		//Halfedges
		new_out->set_tnvef(new_in,he_n2,new_v,new_e,old_f);
		new_in->set_tnvef(new_out,he,old_v_left,new_e,new_f);
		split_out->set_tnvef(split_in,twin->next,new_v,split_e,bound_f);
		split_in->set_tnvef(split_out,new_out,old_v_bot_right,split_e,old_f);

		//Edges
		new_e->halfedge = new_out;
		new_e->sharp = false; //inside edge
		split_e->halfedge = split_in;
		split_e->sharp = true; //boundary edge
	

		//Vertices
		new_v->position = e->center();
		new_v->halfedge = he;

		//Faces
		new_f->boundary = false;
		new_f->halfedge = he;

		//Phase 3: Connect up old elements

		//Halfedges
		he->face = new_f;
		he->vertex = new_v;
		twin->next = split_out;
		he_n->next = new_in;
		he_n->face = new_f;
		he_p->next = split_in;

		//Vertices
		old_v_left->halfedge = new_in;
		old_v_right->halfedge = he_n;
		old_v_bot_right->halfedge = split_in;

		//Faces
		old_f->halfedge = split_in;
		bound_f->halfedge = twin;

		//Debugging


		return new_v;

	}else{ //Normal case

		//Phase 1: Collect all needed elements

		//Get half-edges
		HalfedgeRef he_n = he->next;
		HalfedgeRef he_n2 = he_n->next;
		HalfedgeRef he_p = he;
		do{
			he_p = he_p->next;
		}while(he_p->next != he);

		HalfedgeRef twin_n = twin->next;
		HalfedgeRef twin_n2 = twin_n->next;
		HalfedgeRef twin_p = twin;
		do{
			twin_p = twin_p->next;
		}while(twin_p->next != twin);

		//Get vertices
		VertexRef old_v_left = he_n2->vertex;
		VertexRef old_v_right = twin_n2->vertex;
		VertexRef old_top = twin->vertex;	
		VertexRef old_bot = he->vertex;	

		//Get faces
		FaceRef old_f_left = he->face;
		FaceRef old_f_right = twin->face;


		//Phase 2: Create new elements

		//Half-edges first
		HalfedgeRef new_l_out = emplace_halfedge();
		HalfedgeRef new_l_in = emplace_halfedge();
		HalfedgeRef new_r_out = emplace_halfedge();
		HalfedgeRef new_r_in = emplace_halfedge();
		HalfedgeRef new_out = emplace_halfedge();
		HalfedgeRef new_in = emplace_halfedge();

		//Edges next
		EdgeRef new_l = emplace_edge();
		EdgeRef new_r = emplace_edge();
		EdgeRef new_split = emplace_edge();

		//Vertices
		VertexRef new_v = emplace_vertex();

		//Faces
		FaceRef new_f_left = emplace_face();
		FaceRef new_f_right = emplace_face();

		//Set data for new elements

		//Half edges
		new_l_out->set_tnvef(new_l_in,he_n2,new_v,new_l,old_f_left);
		new_l_in->set_tnvef(new_l_out,new_out,old_v_left,new_l,new_f_left);
		new_r_out->set_tnvef(new_r_in,twin_n2,new_v,new_r,old_f_right);
		new_r_in->set_tnvef(new_r_out,twin,old_v_right,new_r,new_f_right);
		new_out->set_tnvef(new_in,he_n,new_v,new_split,new_f_left);
		new_in->set_tnvef(new_out,new_r_out,old_top,new_split,old_f_right);

		//Edges
		new_l->halfedge = new_l_out;
		new_l->sharp = false; //inner edge 
		new_r->halfedge = new_r_out;
		new_r->sharp = false; //also inner edge
		new_split->halfedge = new_out;
		new_split->sharp = false; //also inner edge

		//Vertices
		new_v->position = e->center(); //actually update position to midpt of edge
		new_v->halfedge = new_l_out;

		//Faces
		new_f_left->halfedge = new_out;
		new_f_left->boundary = false;
		new_f_right->halfedge = twin;
		new_f_right->boundary = false;

		//Phase 3: Connect up old elements

		//Half edges first
		he->set_tnvef(twin,new_l_out,old_bot,e,old_f_left);
		he_n->next = new_l_in;
		he_n->face = new_f_left;
		he_p->next = he;
		
		twin->set_tnvef(he,twin_n,new_v,e,new_f_right);
		twin_n->next = new_r_in;
		twin_n->face = new_f_right;
		twin_p->next = new_in;

		//Vertices
		old_v_left->halfedge = new_l_in;
		old_v_right->halfedge = new_r_in;
		old_top->halfedge = new_in;
		old_bot->halfedge = he;

		//Faces
		old_f_left->halfedge = he;
		old_f_right->halfedge = new_in;
		new_f_left->halfedge = new_out;
		new_f_right->halfedge = twin;

		return new_v;
	}
}



/*
 * inset_vertex: divide a face into triangles by placing a vertex at f->center()
 *  f: the face to add the vertex to
 *
 * returns:
 *  std::nullopt if insetting a vertex would make mesh invalid
 *  the inset vertex otherwise
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
	// A2Lx4 (OPTIONAL): inset vertex
	
	(void)f;
    return std::nullopt;
}


/* [BEVEL NOTE] Note on the beveling process:

	Each of the bevel_vertex, bevel_edge, and extrude_face functions do not represent
	a full bevel/extrude operation. Instead, they should update the _connectivity_ of
	the mesh, _not_ the positions of newly created vertices. In fact, you should set
	the positions of new vertices to be exactly the same as wherever they "started from."

	When you click on a mesh element while in bevel mode, one of those three functions
	is called. But, because you may then adjust the distance/offset of the newly
	beveled face, we need another method of updating the positions of the new vertices.

	This is where bevel_positions and extrude_positions come in: these functions are
	called repeatedly as you move your mouse, the position of which determines the
	amount / shrink parameters. These functions are also passed an array of the original
	vertex positions, stored just after the bevel/extrude call, in order starting at
	face->halfedge->vertex, and the original element normal, computed just *before* the
	bevel/extrude call.

	Finally, note that the amount, extrude, and/or shrink parameters are not relative
	values -- you should compute a particular new position from them, not a delta to
	apply.
*/

/*
 * bevel_vertex: creates a face in place of a vertex
 *  v: the vertex to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(VertexRef v) {
	//A2Lx5 (OPTIONAL): Bevel Vertex
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in bevel_vertex_helper (A2Lx5h)

	(void)v;
    return std::nullopt;
}

/*
 * bevel_edge: creates a face in place of an edge
 *  e: the edge to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(EdgeRef e) {
	//A2Lx6 (OPTIONAL): Bevel Edge
	// Reminder: This function does not update the vertex positions.
	// remember to also fill in bevel_edge_helper (A2Lx6h)

	(void)e;
    return std::nullopt;
}

/*
 * extrude_face: creates a face inset into a face
 *  f: the face to inset
 *
 * returns: reference to the inner face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::extrude_face(FaceRef f) {
	//A2L4: Extrude Face
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in extrude_helper (A2L4h)

	(void)f;
    return std::nullopt;
}

/*
 * flip_edge: rotate non-boundary edge ccw inside its containing faces
 *  e: edge to flip
 *
 * if e is a boundary edge, does nothing and returns std::nullopt
 * if flipping e would create an invalid mesh, does nothing and returns std::nullopt
 *
 * otherwise returns the edge, post-rotation
 *
 * does not create or destroy mesh elements.
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(EdgeRef e) {
	//A2L1: Flip Edge

	//Can't flip boundary edges
	if(e->on_boundary()) return std::nullopt;

	//Otherwise, there should be at least two triangles on each side

	//First get everything we need access to

	//Half edges first
	HalfedgeRef he = e->halfedge;
	HalfedgeRef twin = he->twin;

	HalfedgeRef he_back = he;
	do{
		he_back = he_back->next;
	}while(he_back->next != he);

	HalfedgeRef he_front = twin->next;
	HalfedgeRef he_front2 = he_front->next;

	HalfedgeRef twin_back = twin;
	do{
		twin_back = twin_back->next;
	}while(twin_back->next != twin);
	
	HalfedgeRef twin_front = he->next;
	HalfedgeRef twin_front2 = twin_front->next;

	//Get vertices
	VertexRef v_old_start = he->vertex;
	VertexRef v_old_end = twin->vertex;
	VertexRef v_new_start = he_front2->vertex;
	VertexRef v_new_end = twin_front2->vertex;

	//Get faces
	FaceRef fl = he->face;
	FaceRef fr = twin->face;

	//Now remove edge e
	he_back->next = he_front;
	twin_back->next = twin_front;

	he_front->vertex = v_old_start;
	v_old_start->halfedge = he_front;

	twin_front->vertex = v_old_end;
	v_old_end->halfedge = twin_front;

	//Re-attach edge e in rotated place
	v_new_start->halfedge = he_front2;
	v_new_end->halfedge = twin_front2;

	he_front->next = he;
	he->next = twin_front2;

	twin_front->next = twin;
	twin->next = he_front2;

	//Re-attach face info
	fl->halfedge = he;
	fr->halfedge = twin;

	he_front->face = fl;
	twin_front->face = fr;

	he->vertex = v_new_start;
	twin->vertex = v_new_end;

	return e;
}


/*
 * make_boundary: add non-boundary face to boundary
 *  face: the face to make part of the boundary
 *
 * if face ends up adjacent to other boundary faces, merge them into face
 *
 * if resulting mesh would be invalid, does nothing and returns std::nullopt
 * otherwise returns face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::make_boundary(FaceRef face) {
	//A2Lx7: (OPTIONAL) make_boundary

	return std::nullopt; //TODO: actually write this code!
}

/*
 * dissolve_vertex: merge non-boundary faces adjacent to vertex, removing vertex
 *  v: vertex to merge around
 *
 * if merging would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_vertex(VertexRef v) {
	// A2Lx1 (OPTIONAL): Dissolve Vertex

    return std::nullopt;
}

/*
 * dissolve_edge: merge the two faces on either side of an edge
 *  e: the edge to dissolve
 *
 * merging a boundary and non-boundary face produces a boundary face.
 *
 * if the result of the merge would be an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_edge(EdgeRef e) {
	// A2Lx2 (OPTIONAL): dissolve_edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data
	
    return std::nullopt;
}

/* collapse_edge: collapse edge to a vertex at its middle
 *  e: the edge to collapse
 *
 * if collapsing the edge would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(EdgeRef e) {
	//A2L3: Collapse Edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

	HalfedgeRef he = e->halfedge;
	HalfedgeRef twin = he->twin;

	//First check if edge has a boundary face on some side
	bool he_boundary = he->face->boundary;
	bool twin_boundary = twin->face->boundary;
	if(twin_boundary){

		//Handle tri case differently than v>=4 case

		//Phase 1: Collect all elements

		//Half edges
		HalfedgeRef he_n = he->next;
		HalfedgeRef he_n_twin = he_n->twin;
		HalfedgeRef he_p = he;
		do{
			he_p = he_p->next;
		}while(he_p->next != he);
		HalfedgeRef he_p_twin = he_p->twin;

		HalfedgeRef twin_n = twin->next;
		HalfedgeRef twin_p = he;
		do{
			twin_p = twin_p->next->twin;
		}while(twin_p->next != twin);

		HalfedgeRef he_n_twin_n = he_n_twin->next;
		HalfedgeRef he_n_twin_p = he_n;
		do{
			he_n_twin_p = he_n_twin_p->next->twin;
		}while(he_n_twin_p->next != he_n_twin);

		//HalfedgeRef he_p_twin_n = he_p_twin->next;
		HalfedgeRef he_p_twin_p = twin;
		do{
			he_p_twin_p = he_p_twin_p->next->twin;
		}while(he_p_twin_p->next != he_p_twin);

		//Edges
		EdgeRef old_e_left = he_n->edge;
		EdgeRef old_e_top = he->edge;
		EdgeRef old_e_right = he_p->edge;

		//Vertices
		VertexRef old_v_left = he_n->vertex;
		VertexRef old_v_bot = he_p->vertex;
		VertexRef old_v_right = he->vertex;

		//Faces
		FaceRef old_face = he->face;
		//FaceRef old_face_left = he_n_twin->face;
		//FaceRef old_face_right = he_p_twin->face;
		FaceRef bound_face = twin->face;

		//Data
		Vec3 e_center = e->center();

		assert(old_face->degree() >= 3);

		std::cout << "wuts it now: " << he_n_twin_n->id << std::endl;


		if(old_face->degree() == 3){ //Tri case


			//First check if collapse is valid: is the face a "thin face",creates invalid edge on collapse?
		
			if(he_n->twin->face->boundary && he_p->twin->face->boundary) return std::nullopt;

			//Special case: change halfedges in nbhd of old_v_right to old_v_left
			HalfedgeRef start = old_v_right->halfedge;
			HalfedgeRef cur = start;
			do{	
				cur->vertex = old_v_left;
				cur = cur->twin->next;
			}while(cur != start);
			
			//Phase 2: Delete anything we need to

			//Half edges
			erase_halfedge(he);
			erase_halfedge(twin);
			erase_halfedge(he_n);
			erase_halfedge(he_p);

			//Edges
			erase_edge(old_e_right);
			erase_edge(old_e_top);

			//Vertices
			erase_vertex(old_v_right);

			//Faces
			erase_face(old_face);

			//Phase 3: Connect up what's left

			//Half edges
			twin_p->next = twin_n;
			he_n_twin->twin = he_p_twin;
			he_p_twin->twin = he_n_twin;
			he_p_twin->edge = old_e_left;
	
			//Edges
			old_e_left->halfedge = he_n_twin;

			//Vertices
			old_v_left->halfedge = twin_n;
			old_v_bot->halfedge = he_n_twin;
			old_v_left->position = e_center;

			std::cout << "old_v_left: " << old_v_left->id << std::endl;
			std::cout << "old_v_left_he: " << old_v_left->halfedge->id << std::endl;

			//Faces
			bound_face->halfedge = twin_p;


			return old_v_left;

		}else{ //Unconditional collapse in this case (v>=4 case)

			//Phase 2: Delete anything

			//Special case: change halfedges in nbhd of old_v_right to old_v_left
			HalfedgeRef start = old_v_right->halfedge;
			HalfedgeRef cur = start;
			do{	
				cur->vertex = old_v_left;
				cur = cur->twin->next;
			}while(cur != start);
			
			//Half edges
			erase_halfedge(he);
			erase_halfedge(twin);

			//Edges
			erase_edge(e);

			//Vertices
			erase_vertex(old_v_right);

			//Faces

			//Phase 3: Connect up what's left

			//Half edges
			twin_p->next = twin_n;
			he_p->next = he_n;

			//Edges

			//Vertices
			old_v_left->halfedge = he_n;
			old_v_left->position = e_center;

			//Faces
			old_face->halfedge = he_n;
			bound_face->halfedge = twin_p;

			return old_v_left;
		}
		
	}else if(he_boundary){
		std::cout << "HE BOUNDARY\n";

		e->halfedge = twin;
		return collapse_edge(e);
	}else{

		


		//Phase 1: Collect elements

		//Half edges
		HalfedgeRef he_n = he->next;
		HalfedgeRef he_p = he;
		do{
			he_p = he_p->next;
		}while(he_p->next != he);
		
		HalfedgeRef twin_n = twin->next;
		HalfedgeRef twin_p = twin;
		do{
			twin_p = twin_p->next;
		}while(twin_p->next != twin);
		HalfedgeRef twin_p_twin = twin_p->twin;

		HalfedgeRef he_n_twin = he_n->twin;

		HalfedgeRef he_p_twin = he_p->twin;
		//HalfedgeRef he_p_twin_n = he_p_twin->next;
		HalfedgeRef he_p_twin_p = he_p;
		do{
			he_p_twin_p = he_p_twin_p->next->twin;
		}while(he_p_twin_p->next != he_p_twin);
		
		HalfedgeRef twin_n_twin = twin_n->twin;
		//HalfedgeRef twin_n_twin_n = twin_n_twin->next;
		HalfedgeRef twin_n_twin_p = twin_n;
		do{
			twin_n_twin_p = twin_n_twin_p->next->twin;
		}while(twin_n_twin_p->next != twin_n_twin);

		

		//Edges
		//EdgeRef e_collapse = e;
		EdgeRef top_left = he_n->edge;
		EdgeRef bot_left = he_p->edge;
		EdgeRef top_right = twin_p->edge;
		EdgeRef bot_right = twin_n->edge;

		//Vertices
		//VertexRef v_top_left = he_n->next->vertex;
		VertexRef v_bot_left = he_p->vertex;
		//VertexRef v_top_right = twin_p->vertex;
		VertexRef v_bot_right = twin_n->next->vertex;
		VertexRef v_top = he_n->vertex;
		VertexRef v_bot = he->vertex;

		//Faces
		FaceRef f_left = he->face;
		FaceRef f_right = twin->face;
		//FaceRef f_he_n = he_n_twin->face;
		//FaceRef f_twin_n = twin_n_twin->face;

		//Data
		Vec3 e_center = e->center();

		//First check for hourglass shape case

		bool top_bound = false;
		bool bot_bound = false;
		HalfedgeRef cur = v_top->halfedge;
		do{
			if(cur->face->boundary) top_bound = true;
			cur = cur->twin->next;
		}while(cur != v_top->halfedge);
		cur = v_bot->halfedge;
		do{
			if(cur->face->boundary) bot_bound = true;
			cur = cur->twin->next;
		}while(cur != v_bot->halfedge);

		//Hourglass means each endpoint of edge has a boundary face in its orbit
		if(top_bound && bot_bound) return std::nullopt;

		//3 cases: tri & tri, tri & not tri,not tri & not tri
		uint32_t size_left = f_left->degree();
		uint32_t size_right = f_right->degree();

		if((size_left == 3 ) && (size_right == 3)){ //Tri | Tri 

			//Phase 2: Delete stuff
			
			//Before deleting, update half edges to pt to collapsed vertex v_top
			HalfedgeRef cur = v_bot->halfedge;
			do{
				cur->vertex = v_top;
				cur = cur->twin->next;
			}while(cur != v_bot->halfedge);

			//Half edges
			erase_halfedge(he);
			erase_halfedge(he_n);
			erase_halfedge(he_p);

			erase_halfedge(twin);
			erase_halfedge(twin_p);
			erase_halfedge(twin_n);			

			//Edges
			erase_edge(e);
			erase_edge(bot_left);
			erase_edge(bot_right);

			//Vertices
			erase_vertex(v_bot);

			//Faces
			erase_face(f_left);
			erase_face(f_right);

			//Phase 3: Connect up stuff

			//Half edges
			he_n_twin->twin = he_p_twin;
			he_n_twin->edge = top_left;
			he_p_twin->twin = he_n_twin;
			he_p_twin->edge = top_left;

			twin_n_twin->twin = twin_p_twin;
			twin_n_twin->edge = top_right;
			twin_p_twin->twin = twin_n_twin;
			twin_p_twin->edge = top_right;

			//Edges
			top_left->halfedge = he_n_twin;
			top_right->halfedge = twin_p_twin;

			//Vertices
			v_top->halfedge = twin_p_twin;
			v_top->position = e_center;
			v_bot_left->halfedge = he_n_twin;
			v_bot_right->halfedge = twin_n_twin;

			//Faces

			return v_top;

		}else if((size_left == 3) && (size_right != 3)){ //Tri | Face

			//Phase 2: delete stuff

			//First move all bot halfedges to have top vert due to collapse
			HalfedgeRef start = v_bot->halfedge;
			HalfedgeRef cur = start;
			do{
				cur->vertex = v_top;
				cur = cur->twin->next;
			}while(cur != start);

			//Half edges
			erase_halfedge(he);
			erase_halfedge(he_n);
			erase_halfedge(he_p);
			erase_halfedge(twin);

			//Edges
			erase_edge(e);
			erase_edge(bot_left);

			//Vertices
			erase_vertex(v_bot);

			//Faces
			erase_face(f_left);

			//Data
			v_top->position = e_center;

			//Phase 3: connect up stuff

			//Half edges
			he_n_twin->twin = he_p_twin;
			he_n_twin->edge = top_left;
			he_p_twin->twin = he_n_twin;
			he_p_twin->edge = top_left;
			twin_p->next = twin_n;

			//Edges
			top_left->halfedge = he_n_twin;

			//Vertices
			v_bot_left->halfedge = he_n_twin;
			v_top->halfedge = he_p_twin;

			//Faces
			f_right->halfedge = twin_p;

			return v_top;
		}else if((size_left != 3) && (size_right == 3)){ // Face | Face
			e->halfedge = twin;
			return collapse_edge(e);
		}else{ //General non-tri case

			//Part 2: Delete stuff

			//First move all bot halfedges to have top vert due to collapse
			HalfedgeRef start = v_bot->halfedge;
			HalfedgeRef cur = start;
			do{
				cur->vertex = v_top;
				cur = cur->twin->next;
			}while(cur != start);

			//Half edges
			erase_halfedge(he);
			erase_halfedge(twin);

			//Edges
			erase_edge(e);

			//Vertices
			erase_vertex(v_bot);

			//Faces

			//Data
			v_top->position = e_center;

			//Part 3: connect stuff

			//Half edges
			he_p->next = he_n;
			twin_p->next = twin_n;

			//Edges

			//Vertices
			v_top->halfedge = he_n;

			//Faces
			f_left->halfedge = he_n;
			f_right->halfedge = twin_n;

			return v_top;
		}
	}

	
}

/*
 * collapse_face: collapse a face to a single vertex at its center
 *  f: the face to collapse
 *
 * if collapsing the face would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(FaceRef f) {
	//A2Lx3 (OPTIONAL): Collapse Face

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

    return std::nullopt;
}

/*
 * weld_edges: glue two boundary edges together to make one non-boundary edge
 *  e, e2: the edges to weld
 *
 * if welding the edges would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns e, updated to represent the newly-welded edge
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::weld_edges(EdgeRef e, EdgeRef e2) {
	//A2Lx8: Weld Edges

	//Reminder: use interpolate_data() to merge bone_weights data on vertices!

    return std::nullopt;
}



/*
 * bevel_positions: compute new positions for the vertices of a beveled vertex/edge
 *  face: the face that was created by the bevel operation
 *  start_positions: the starting positions of the vertices
 *     start_positions[i] is the starting position of face->halfedge(->next)^i
 *  direction: direction to bevel in (unit vector)
 *  distance: how far to bevel
 *
 * push each vertex from its starting position along its outgoing edge until it has
 *  moved distance `distance` in direction `direction`. If it runs out of edge to
 *  move along, you may choose to extrapolate, clamp the distance, or do something
 *  else reasonable.
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after bevel_vertex or bevel_edge.
 * (So you can assume the local topology is set up however your bevel_* functions do it.)
 *
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::bevel_positions(FaceRef face, std::vector<Vec3> const &start_positions, Vec3 direction, float distance) {
	//A2Lx5h / A2Lx6h (OPTIONAL): Bevel Positions Helper
	
	// The basic strategy here is to loop over the list of outgoing halfedges,
	// and use the preceding and next vertex position from the original mesh
	// (in the start_positions array) to compute an new vertex position.
	
}

/*
 * extrude_positions: compute new positions for the vertices of an extruded face
 *  face: the face that was created by the extrude operation
 *  move: how much to translate the face
 *  shrink: amount to linearly interpolate vertices in the face toward the face's centroid
 *    shrink of zero leaves the face where it is
 *    positive shrink makes the face smaller (at shrink of 1, face is a point)
 *    negative shrink makes the face larger
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after extrude_face.
 * (So you can assume the local topology is set up however your extrude_face function does it.)
 *
 * Using extrude face in the GUI will assume a shrink of 0 to only extrude the selected face
 * Using bevel face in the GUI will allow you to shrink and increase the size of the selected face
 * 
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::extrude_positions(FaceRef face, Vec3 move, float shrink) {
	//A2L4h: Extrude Positions Helper

	//General strategy:
	// use mesh navigation to get starting positions from the surrounding faces,
	// compute the centroid from these positions + use to shrink,
	// offset by move
	
}

