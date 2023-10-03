# pTrace

A physically-based path tracer based based on party of [15-462/662 Computer Graphics](http://15462.courses.cs.cmu.edu).

User documentation and example works are available at https://cmu-graphics.github.io/Scotty3D-docs/ .

Setup is described [here](https://github.com/eqdiag/myScotty3d/tree/main#readme).


# Features
- [x] Triangle Mesh and Sphere Intersections
- [x] Bounding Volume Hierachy (BVH) for Ray Acceleration
- [x] Recursive Path Tracing solution to the [Rendering Equation](https://en.wikipedia.org/wiki/Rendering_equation).
- [x] Support for multiple material BSDFs: Lambertian (diffuse), Perfect mirror, Dielectric (glass)
- [x] Direct light sampling to speed up convergence
- [x] Single-sample multiple-importance sampling between BSDF and direct light samples
- [x] Support for custom environment maps as light sources

# Screenshots

## Comparison (BRDF sampling vs Direct Light sampling w/ 32 samples per pixel)
### 32 samples per pixel (No direct lighting)
![compare_no_direct_light](/screenshots/cbox_compare_no_direct.png "No Direct lighting cbox")
### 32 samples per pixel (With direct lighting)
![compare_direct_light](/screenshots/cbox_compare_direct.png "Direct lighting cbox")


## Comparison (BRDF sampling vs Direct Light sampling w/ 32 samples per pixel for Environment Map)
![env_no_direct_light](/screenshots/env_compare_no_direct.png "No Direct lighting with Environment Map")
![env_direct_light](/screenshots/env_compare_direct.png "Direct lighting with Environment Map")

# Showcase

## Cornell Box (All Diffuse Materials,512 spp)
![cbox_all_diffuse](/screenshots/cbox_diffuse.png "Cbox all diffuse")

## Cornell Box (Perfect Reflective & Glass Materials,512 spp)
![cbox_reflect_glass](/screenshots/cbox_reflect_glass.png "Cbox reflect and glass")

## Cornell Box (Bunny Triangle Mesh, 512 spp)
![cbox_reflect_glass](/screenshots/cbox_bunny.png "Cbox bunny")

## Environment Map with Mesh (Dragon Triangle Mesh, 512 spp)
![env_dragon](/screenshots/env_ennis_dragon.png "Environment dragon")