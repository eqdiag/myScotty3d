
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <stack>

#include <iostream>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims; ///< number of primitives in the bucket
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh

	// Keep these
    nodes.clear();
    primitives = std::move(prims);

    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.

	//TODO
	//std::cout << "\nBULDING BVH... LEAF_SIZE: " << max_leaf_size << std::endl;
	root_idx = buildHelper(0,n_primitives(),max_leaf_size);
}

template<typename Primitive>
size_t BVH<Primitive>::buildHelper(size_t start,size_t size,size_t max_leaf_size){

	size_t nodeId;

	if(size <= max_leaf_size){ //Leaf node case
		//std::cout << "\n LEAF CASE: [" << start << "," << start + size << ")\n";
		BBox leafBox;
		for(size_t i = start; i < start + size;i++){
			auto& primitive = primitives.at(i);
			leafBox.enclose(primitive.bbox());
		}
		//Note: l = r in this case
		nodeId = new_node(leafBox,start,size);
	}else{ 
		//Internal node case, so we need to
		//1) Find a good way to cut primitives in sublist
		//2) Recurse on children
		//3) Create parent node and return


		float smallestCost = FLT_MAX;
		int smallestAxis = 0;
		float sliceValue = 0;

		//For each axis
		for(int i = 0;i<3;i++){
			//Init buckets
			std::vector<SAHBucketData> buckets(NUM_BUCKETS,SAHBucketData{});


			float minI = FLT_MAX;
			float maxI = -FLT_MAX;
			//Find min/max extent along axis i
			for(int j = start; j < start + size;j++){
				auto& primitive = primitives.at(j);
				Vec3 center = primitive.bbox().center();
				if(center[i] < minI){
					minI = center[i];
				}
				if(center[i] > maxI){
					maxI = center[i];
				}
			}

			float bucketSize = (maxI - minI) / static_cast<float>(NUM_BUCKETS);
			//Make buckets slightly bigger to resolve edge cases (when centers lie on boundaries)
			bucketSize += 0.0001;
			
			//Drop primitives in buckets
			for(int j = start; j < start + size;j++){
				auto& primitive = primitives.at(j);
				auto box = primitive.bbox();
				Vec3 center = box.center();
				//Find bucket the primitive lands in
				size_t bucketId = static_cast<size_t>(fabs(center[i] - minI)/bucketSize);
				//std::cout << "DIST: " << fabs(center[i] - minI) << std::endl;
				//std::cout << "RANGE: " << maxI - minI << std::endl;
				//std::cout << "bucket id: " << bucketId << std::endl;
				auto& bucket = buckets.at(bucketId);
				//Update bucket bbox and counts
				bucket.bb.enclose(box);
				bucket.num_prims++;
			}



			//Check all NUM_BUCKETS - 1 partitions on axis i
			
			for(int k = 1;k < NUM_BUCKETS;k++){
				BBox left;
				BBox right;
				size_t numLeft = 0;
				size_t numRight = 0;

				//Compute left bbox
				for(int j = 0;j < k;j++){
					left.enclose(buckets.at(j).bb);
					numLeft += buckets.at(j).num_prims;
				}

				//Compute right bbox
				for(int j = k;j < NUM_BUCKETS;j++){
					right.enclose(buckets.at(j).bb);
					numRight += buckets.at(j).num_prims;
				}

				BBox parent;
				parent.enclose(left);
				parent.enclose(right);
				float parentSA = parent.surface_area();

				//Partition cost
				float pCost = (left.surface_area() / parentSA)*numLeft;
				pCost += (right.surface_area() / parentSA)*numRight;
				//std::cout << "\nLEFT SA: " << left.surface_area() / parentSA << std::endl;
				//std::cout << "\nRIGHT SA: " << right.surface_area() / parentSA << std::endl;

				//Update if we found a better partition
				if(pCost < smallestCost){
					smallestCost = pCost;
					smallestAxis = i;
					sliceValue = minI + k*bucketSize;
				}
			}
		}

		//Now for the splitting...
		auto first = primitives.begin() + start;
		auto last = first + size; //Because partition is exclusive on right end
		auto mid = std::partition(first,last,[smallestAxis,sliceValue](auto& primitive){
			return primitive.bbox().center()[smallestAxis] <= sliceValue;
		});

		auto newFirst = primitives.begin() + start;
		size_t leftSize = std::distance(newFirst,mid);
		size_t rightSize = size - leftSize;
		
		//Create children
		size_t leftNode = buildHelper(start,leftSize,max_leaf_size);
		size_t rightNode = buildHelper(start + leftSize,rightSize,max_leaf_size);

		BBox parentBox;
		parentBox.enclose(nodes.at(leftNode).bbox);
		parentBox.enclose(nodes.at(rightNode).bbox);

		nodeId = new_node(parentBox,start,size,leftNode,rightNode);
	}

	return nodeId;
}


template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.
	Trace ret;
	ret.hit = false;
	//TODO: replace this code with a more efficient traversal:
    /*for(const Primitive& prim : primitives) {
        Trace hit = prim.hit(ray);
        ret = Trace::min(ret, hit);
    }*/
	if(nodes.size() > 0){
		auto root = nodes.at(root_idx);
		ret = hitHelper(ray,root);
	}

	return ret;	
}

template<typename Primitive> Trace BVH<Primitive>::hitHelper(const Ray& ray,const Node& node) const{

	//First check if ray hits bbox
	Vec2 times(0.0,FLT_MAX);
	Trace ret;
	ret.hit = false;
	//If hits node bbox, look inside node
	if(node.bbox.hit(ray,times)){

		//In leaf case, just traverse primitives inside
		if(node.is_leaf()){
			for(int i = node.start;i < node.start + node.size;i++){
				const auto& primitive = primitives.at(i);
				ret = Trace::min(ret,primitive.hit(ray));
			}
		}else{ //Interior node case

			//Check left and right nodes, take min
			Trace leftRet = hitHelper(ray,nodes.at(node.l));
			Trace rightRet = hitHelper(ray,nodes.at(node.r));

			ret = Trace::min(ret,leftRet);
			ret = Trace::min(ret,rightRet);
		}
	}

	return ret;
}


template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
