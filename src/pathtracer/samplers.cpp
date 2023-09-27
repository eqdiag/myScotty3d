
#include "samplers.h"
#include "../scene/shape.h"
#include "../util/rand.h"

#include <iostream>

constexpr bool IMPORTANCE_SAMPLING = true;

namespace Samplers {

Vec2 Rect::sample(RNG &rng) const {
	//A3T1 - step 2 - supersampling

    // Return a point selected uniformly at random from the rectangle [0,size.x)x[0,size.y)
    // Useful function: rng.unit()
    return Vec2{rng.unit()*size.x,rng.unit()*size.y};
}

float Rect::pdf(Vec2 at) const {
	if (at.x < 0.0f || at.x > size.x || at.y < 0.0f || at.y > size.y) return 0.0f;
	return 1.0f / (size.x * size.y);
}

Vec3 Point::sample(RNG &rng) const {
	return point;
}

float Point::pdf(Vec3 at) const {
	return at == point ? 1.0f : 0.0f;
}

Vec3 Triangle::sample(RNG &rng) const {
	float u = std::sqrt(rng.unit());
	float v = rng.unit();
	float a = u * (1.0f - v);
	float b = u * v;
	return a * v0 + b * v1 + (1.0f - a - b) * v2;
}

float Triangle::pdf(Vec3 at) const {
	float a = 0.5f * cross(v1 - v0, v2 - v0).norm();
	float u = 0.5f * cross(at - v1, at - v2).norm() / a;
	float v = 0.5f * cross(at - v2, at - v0).norm() / a;
	float w = 1.0f - u - v;
	if (u < 0.0f || v < 0.0f || w < 0.0f) return 0.0f;
	if (u > 1.0f || v > 1.0f || w > 1.0f) return 0.0f;
	return 1.0f / a;
}

Vec3 Hemisphere::Uniform::sample(RNG &rng) const {

	float Xi1 = rng.unit();
	float Xi2 = rng.unit();

	float theta = std::acos(Xi1);
	float phi = 2.0f * PI_F * Xi2;

	float xs = std::sin(theta) * std::cos(phi);
	float ys = std::cos(theta);
	float zs = std::sin(theta) * std::sin(phi);

	return Vec3(xs, ys, zs);
}

float Hemisphere::Uniform::pdf(Vec3 dir) const {
	if (dir.y < 0.0f) return 0.0f;
	return 1.0f / (2.0f * PI_F);
}

Vec3 Hemisphere::Cosine::sample(RNG &rng) const {

	float phi = rng.unit() * 2.0f * PI_F;
	float cos_t = std::sqrt(rng.unit());

	float sin_t = std::sqrt(1 - cos_t * cos_t);
	float x = std::cos(phi) * sin_t;
	float z = std::sin(phi) * sin_t;
	float y = cos_t;

	return Vec3(x, y, z);
}

float Hemisphere::Cosine::pdf(Vec3 dir) const {
	if (dir.y < 0.0f) return 0.0f;
	return dir.y / PI_F;
}

Vec3 Sphere::Uniform::sample(RNG &rng) const {
	//A3T7 - sphere sampler

    // Generate a uniformly random point on the unit sphere.
    // Tip: start with Hemisphere::Uniform
	Vec3 dir = hemi.sample(rng);
	//With 0.5 prob flip the direction
	return rng.coin_flip(0.5) ? -dir : dir;
}

float Sphere::Uniform::pdf(Vec3 dir) const {
	return 1.0f / (4.0f * PI_F);
}

Sphere::Image::Image(const HDR_Image& image) {
    //A3T7 - image sampler init

    // Set up importance sampling data structures for a spherical environment map image.
    // You may make use of the _pdf, _cdf, and total members, or create your own.

    const auto [_w, _h] = image.dimension();
    w = _w;
    h = _h;

	//First compute the unnormalized pdf at each pixel
	for(int row = 0; row < h;row++){
		for(int col = 0;col < w;col++){
			float L = image.at(col,row).luma();
			float theta = (1.0f - (static_cast<float>(row)/static_cast<float>(h - 1)))*PI_F;
			_pdf.push_back(std::sin(theta)*L);
		}
	}

	//Normalize
	float norm = 0.0f;
	for(auto v: _pdf) norm += v;
	for(auto& v: _pdf) v /= norm;

	//Compute cdf
	float total = 0.0f;
	for(auto v: _pdf){
		total += v;
		_cdf.push_back(total);
	}
}

Vec3 Sphere::Image::sample(RNG &rng) const {
	if(!IMPORTANCE_SAMPLING) {
		// Step 1: Uniform sampling
		// Declare a uniform sampler and return its sample
		Sphere::Uniform sampler;
    	return sampler.sample(rng);
	} else {
		// Step 2: Importance sampling
		// Use your importance sampling data structure to generate a sample direction.
		// Tip: std::upper_bound
		float r = rng.unit();

		//Sphere::Uniform sampler;
    	//return sampler.sample(rng);

		auto cut = std::upper_bound(_cdf.begin(),_cdf.end(),r);
		uint32_t index = std::distance(_cdf.begin(),cut);

		uint32_t x =  index % w;
		uint32_t y = index / w;
		float theta = (1.0f - (static_cast<float>(y)/static_cast<float>(h - 1)))*PI_F;
		float phi = (static_cast<float>(x)/static_cast<float>(w - 1)) * 2.0f*PI_F;
    	return Vec3{std::sin(theta)*cos(phi),std::cos(theta),std::sin(theta)*std::sin(phi)};
	}
}

float Sphere::Image::pdf(Vec3 dir) const {
    if(!IMPORTANCE_SAMPLING) {
		// Step 1: Uniform sampling
		// Declare a uniform sampler and return its pdf
    	Sphere::Uniform sampler;
		return sampler.pdf(dir);
	} else {
		// A3T7 - image sampler importance sampling pdf
		// What is the PDF of this distribution at a particular direction?
		//Sphere::Uniform sampler;
		//return sampler.pdf(dir);
		Vec2 uv = Shapes::Sphere::uv(dir);

		float theta = (1.0f - uv.y)*PI_F;
		float sinTheta = std::sin(theta);
		float jacobian = w*h/(2*PI_F*PI_F*sinTheta);

		uint32_t x = static_cast<uint32_t>(uv.x*w);
		uint32_t y = static_cast<uint32_t>(uv.y*h);

		uint32_t index = y*w + x;
		assert(index <= w*h);
		float prob = _pdf[index];
		return prob*jacobian;
	}
}

} // namespace Samplers
