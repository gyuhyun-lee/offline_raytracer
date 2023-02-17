# Offline Raytracer
Simple offline raytracer with full BSDF, written in C++ & Objective-C(for MacOS platform code) from scratch.

# Features
- Multiple object shapes, including arbitrary meshes
- HDR support, depth of field & shaped lights
- Custom OBJ, PLY loader and bounding box hierarchy with octree
- Optimized using SIMD and multi-threading

# Future improvements
- Profile cache usage for maximum performace
- Use better random ray generator instead of XOR Shift(https://en.wikipedia.org/wiki/Xorshift)
- Add x86/64 architecture support
- Add support for more interesting geometries

# Showcase
Showcase #1
![Sample 1](showcase/1.hdr)

Showcase #2(2022/11/16)
![Sample 2](showcase/2.hdr)








