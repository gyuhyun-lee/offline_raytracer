#include "ray.h"

#define Very_Small_Number 0.000001f
#define euler_number 2.71828182845904523536028747135266249f
#define Hit_t_Threshold 0.000001f

// creates a row_major matrix that rotates certain vector along Z
internal m3
rotation_matrix_along_z(v3 source)
{
    m3 result = m3_();

    if(!compare_0(cross(source, v3_(0, 0, 1))))
    {
        // Not parallel
        v3 a = normalize(source);
        v3 b = normalize(cross(v3_(0, 0, 1), a));

        if(compare_0(b))
        {
            // try x
            b = normalize(cross(v3_(1, 0, 0), a));
        }

        v3 c = cross(a, b);

        result.rows[0] = b;
        result.rows[1] = c;
        result.rows[2] = a;
    }

    return result;
}

#if 0
internal void
push_sphere_into_csg(MemoryArena *arena, CSG *csg, v3 sphere_center, f32 sphere_r)
{
    if(!csg.base)
    {
        assert(csg.used == 0);
        csg.base = ((u8 *)arena->base + used);
    }

    Sphere *sphere = ;
}

internal void
push_intersection_into_csg(memoryarena *arena, CSG *csg)
{
}
#endif

struct IntersectionTestResult
{
    f32 hit_t;
    v3 hit_normal;
    b32 inner_hit;
};

// TODO(joon) : This works for both inward & outward normals, 
// but we might only want to test it against the outward normal
internal IntersectionTestResult
ray_intersect_with_triangle(v3 v0, v3 v1, v3 v2, v3 ray_origin, v3 ray_dir)
{
    /*
       NOTE(joon) : 
       Moller-Trumbore line triangle intersection argorithm

       |t| =       1           |(T x E1) * E2|
       |u| = -------------  x  |(ray_dir x E2) * T |
       |v| = (ray_dir x E2)    |(T x E1) * ray_dir |

       where T = ray_origin - v0, E1 = v1 - v0, E2 = v2 - v0,
       ray  = ray_origin + t * ray_dir;

       u & v = barycentric coordinates of the triangle, as a triangle can be represented in a form of 
       (1 - u - v)*v0 + u*v1 + v*v2;

       Note that there are a lot of same cross products, which we can calculate just once and reuse
       Also, v0, v1, v2 can be in any order
    */

    IntersectionTestResult result = {};
    result.hit_t = -1.0f;

    v3 e1 = v1 - v0;
    v3 e2 = v2 - v0;
    v3 cross_ray_e2 = cross(ray_dir, e2);

    f32 det = dot(cross_ray_e2, e1);
    v3 T = ray_origin - v0;

    // TODO(joon) : completely made up number
    // TODO(joon) : Why bumping this up cause werid artifact?
    f32 tolerance = 0.000001f;
    if(det <= -tolerance || det >= tolerance) // if the determinant is 0, it means the ray is parallel to the triangle
    {
        v3 a = cross(T, e1);

        f32 t = dot(a, e2) / det;
        f32 u = dot(cross_ray_e2, T) / det;
        f32 v = dot(a, ray_dir) / det;

        if(t >= Hit_t_Threshold && u >= 0.0f && v >= 0.0f && u+v <= 1.0f)
        {
            result.hit_t = t;
            // TODO(joon): calculate normal based on the ray dir, so that the next normal is facing the incoming ray dir
            // otherwise, the reflection vector will be totally busted?
            result.hit_normal = cross(e1, e2);
        }
    }

    return result;
}

internal b32 
is_inside_sphere(v3 center, f32 r, v3 p)
{
    b32 result = false;

    v3 rel_p = p - center;
    // TODO(joon) do we also need to include =?
    if(length_square(rel_p) < square(r))
    {
        result = true;
    }

    return result;
}

internal IntersectionTestResult
ray_intersect_with_sphere(v3 center, f32 r, v3 ray_origin, v3 ray_dir)
{
    IntersectionTestResult result = {};
    result.hit_t = -1.0f;
    
    v3 rel_ray_origin = ray_origin - center;
    f32 a = dot(ray_dir, ray_dir);
    f32 b = dot(ray_dir, rel_ray_origin);
    f32 c = dot(rel_ray_origin, rel_ray_origin) - r*r;

    f32 root_term = b*b - a*c;

    f32 tolerance = 0.00001f;
    if(root_term >= tolerance)
    {
        f32 sqrt_root_term = sqrtf(root_term);

        // two intersection points
        f32 tn = (-b - sqrt_root_term)/a;
        f32 tp = (-b + sqrt_root_term)/a;

        f32 hit_normal_c = 1.0f;
        f32 t = 0.0f;
        if(tn < 0.0f)
        {
            //ray started inside the sphere
            t = tp;
            //hit_normal_c = -1.0f;
            result.inner_hit = true;
        }
        else
        {
            t = tn;
        }

        if(t > Hit_t_Threshold)
        {
            result.hit_t = t;
            result.hit_normal = hit_normal_c*(ray_origin + result.hit_t * ray_dir - center);
        }
    }
    else if(root_term < tolerance && root_term > -tolerance)
    {
        // one intersection point
        f32 t = (-b)/(2*a);
        if(t > Hit_t_Threshold)
        {
            result.hit_t = t;
            result.hit_normal = ray_origin + result.hit_t * ray_dir - center;
        }
    }
    else
    {
        // no intersection
    }

    return result;
}

/*
    Normal plane can be expressed in : N * P = d, where N being a normal and P being any point in plane.

    Because the planes are axis aligned, we can also simplify the typical ray - plane intersection code.
    Also, we can get all 6 planes(3 slabs) just by using the min & max points of the box.
    For example, when we solve the equation above with a ray, we get t = (d - dot(N, O)) / dot(N, V).
    Because the planes are axis aligned, the normals should be (1, 0, 0), (0, 1, 0), (0, 0, 1).
    So tx = (Px - Ox) / Vx;
    So ty = (Py - Oy) / Vy;
    So tz = (Pz - Oz) / Vz;

    If there is an intersection, a overlapping interval should exist.

*/
internal IntersectionTestResult
ray_intersect_with_aab(v3 min, v3 max, v3 ray_origin, v3 ray_dir)
{
    // TODO(joon) This value can be stored
    v3 inv_ray_dir = v3_(1.0f/ray_dir.x, 1.0f/ray_dir.y, 1.0f/ray_dir.z);

    IntersectionTestResult result = {};
    result.hit_t = -1.0f;

    v3 t0 = hadamard((min - ray_origin), inv_ray_dir); // normals are (-1, 0, 0), (0, -1, 0), (0, 0, -1)
    v3 t1 = hadamard((max - ray_origin), inv_ray_dir);// normals are (1, 0, 0), (0, 1, 0), (0, 0, 1)

    v3 t_min = gather_min_elements(t0, t1);
    v3 t_max = gather_max_elements(t0, t1);

    f32 max_component_of_t_min = max_component(t_min);
    f32 min_component_of_t_max = min_component(t_max);

    /*
       Three lines : goes from tmin to tmax for the slabs(two parallel slabs give us min & max t)
       To intersect, those three lines should be intersecting to each other, which means
       the max of the min should be smaller than the min of the max

        Intersect:
        min------max
            min------max
        Don't intersect:
        min------max
                    min------max
    */
    if(min_component_of_t_max >= max_component_of_t_min)
    {
        f32 t_min_x_t = t0.x;
        v3 t_min_x_normal = v3_(-1, 0, 0);
        if(t_min_x_t > t1.x)
        {
            t_min_x_t = t1.x;
            t_min_x_normal = v3_(1, 0, 0);
        }

        f32 t_min_y_t = t0.y;
        v3 t_min_y_normal = v3_(0, -1, 0);
        if(t_min_y_t > t1.y)
        {
            t_min_y_t = t1.y;
            t_min_y_normal = v3_(0, 1, 0);
        }

        f32 t_min_z_t = t0.z;
        v3 t_min_z_normal = v3_(0, 0, -1);
        if(t_min_z_t > t1.z)
        {
            t_min_z_t = t1.z;
            t_min_z_normal = v3_(0, 0, 1);
        }

        f32 max_of_min_t = t_min_x_t;
        v3 max_of_min_t_normal = t_min_x_normal;

        if(max_of_min_t < t_min_y_t)
        {
            max_of_min_t = t_min_y_t;
            max_of_min_t_normal = t_min_y_normal;
        }

        if(max_of_min_t < t_min_z_t)
        {
            max_of_min_t = t_min_z_t;
            max_of_min_t_normal = t_min_z_normal;
        }

        result.hit_t = max_component_of_t_min;
        result.hit_normal = max_of_min_t_normal;
        //result.hit_normal = v3_(0
    }

    return result;
}


internal IntersectionTestResult
ray_intersect_with_cylinder(v3 base, v3 axis, f32 radius, 
                            v3 ray_origin, v3 ray_dir)
{

    IntersectionTestResult result = {};
    result.hit_t = -1.0f;

    // rotate axis to (0, 0, 1)
    m3 rotation = rotation_matrix_along_z(axis);
    ray_origin = rotation * (ray_origin - base);
    ray_dir = rotation * ray_dir;

    // intersection test with two slabs(top & bottom)
    // as the slabs have normal of (0, 0, 1), we only need to calculate just using the z value
    f32 t_slab_bottom = (-ray_origin.z) / ray_dir.z; 
    f32 t_slab_top = (length(axis) - ray_origin.z) / ray_dir.z;

    f32 t_slab_min = minimum(t_slab_bottom, t_slab_top);
    f32 t_slab_max = maximum(t_slab_bottom, t_slab_top);

    // intersection test with a cylinder that is inf large in z dir
    // put ray equation inside lx^2 + y^2 = r^2
    f32 a = dot(ray_dir.xy, ray_dir.xy);
    f32 b = dot(ray_dir.xy, ray_origin.xy);
    f32 c = dot(ray_origin.xy, ray_origin.xy) - radius * radius;

    f32 det = b*b - a*c;
    if(det >= 0.0f)
    {
        f32 sqrt_det = sqrtf(det);
        f32 t_cylinder_min = (-b - sqrt_det) / a;
        f32 t_cylinder_max = (-b + sqrt_det) / a;

        f32 t_max_in_min = maximum(t_slab_min, t_cylinder_min);
        f32 t_min_in_max = minimum(t_slab_max, t_cylinder_max);

#if 1
        if(t_max_in_min <= t_min_in_max)
        {


            result.hit_t = t_max_in_min;

            result.hit_normal = v3_(0, 1, 0);
            if(t_slab_min < t_cylinder_min)
            {
                if(axis.x != 0.0f)
                {
                    i32 breakerer = 0;
                }
                // update hit normal
                result.hit_normal = v3_((ray_origin + result.hit_t * ray_dir).xy, 0);
                //result.hit_normal = v3_(1, 0, 0);
            }
            m3 inv_rotation = transpose(rotation);
            result.hit_normal = inv_rotation *  result.hit_normal;
        }
#endif
    }
    else
    {
        // no intersection
    }

    return result;
}

#if 0
internal void
tiled_raytrace(World *world, Camera *camera, 
                v3 *output_buffer, i32 output_width, i32 output_height, 
                i32 tile_min_x, i32 tile_min_y,
                i32 tile_one_past_max_x, i32 tile_one_past_max_y)
{
    f32 width_over_height = output_width/(f32)output_height;

    v3 *row = output_buffer + tile_min_y * output_width + tile_min_x;
    for(i32 y = tile_min_y;
            y < tile_one_past_max_y;
            ++y)
    {
        f32 normalized_y = 2.0f * ((y + 0.5f) / (f32)output_height) - 1.0f;

        v3 *pixel = row;
        for(i32 x = tile_min_x;
                x < tile_one_past_max_x;
                ++x)
        {
            f32 normalized_x = 2.0f * ((x + 0.5f) / (f32)output_width) - 1.0f;

            // TODO(joon) What's the math here? 
            //v3 p = width_over_height * normalized_x * camera->x_axis + normalized_y * camera->y_axis - camera->z_axis;
            //v3 ray_dir = p - camera->p; // TODO(joon) don't need to normalize this, as long as we use the same vector for all of our intersection tests?
            v3 ray_dir = normalize(normalized_x * camera->x_axis + 
                                    normalized_y * camera->y_axis - 
                                    camera->z_axis);
            v3 ray_origin = camera->p;

            f32 min_t = Flt_Max;
            u32 hit_mat_index = 0;

            f32 hit_t_threshold = 0.0001f;
#if 1
            for(u32 sphere_index = 0;
                    sphere_index < world->sphere_count;
                    ++sphere_index)
            {
                Sphere *sphere = world->spheres + sphere_index;
                f32 hit_t = ray_intersect_with_sphere(sphere->center, sphere->r, ray_origin, ray_dir);

                if(hit_t >= hit_t_threshold && hit_t < min_t)
                {
                    min_t = hit_t;
                    hit_mat_index = sphere->mat_index;
                }
            }
#endif

#if 1
            for(u32 box_index = 0;
                    box_index < world->box_count;
                    ++box_index)
            {
                AAB *box = world->boxes + box_index;

                f32 hit_t = ray_intersect_with_aab(aab->min, aab->max, ray_origin, ray_dir);
                if(hit_t >= hit_t_threshold && 
                    hit_t < min_t)
                {
                    min_t = hit_t;
                    hit_mat_index = box->mat_index;
                }
            }
#endif

#if 1
            for(u32 cylinder_index = 0;
                    cylinder_index < world->cylinder_count;
                    ++cylinder_index)
            {
                Cylinder *cylinder = world->cylinders + cylinder_index;
                f32 hit_t = ray_intersect_with_cylinder(cylinder->base, cylinder->axis, cylinder->r, ray_origin, ray_dir);

                if(hit_t >= hit_t_threshold && 
                    hit_t < min_t)
                {
                    min_t = hit_t;
                    hit_mat_index = cylinder->mat_index;
                }
            }
#endif

#if 1
            for(u32 mesh_index = 0;
                    mesh_index < world->mesh_count;
                    ++mesh_index)
            {
                Mesh *mesh = world->meshes + mesh_index;

                f32 aabb_hit_t = ray_intersect_with_aab(mesh->aabb_min, mesh->aabb_max, 
                                                   ray_origin, ray_dir);
                if(aabb_hit_t >= hit_t_threshold && 
                    aabb_hit_t < min_t)
                {
                    for(u32 index_index = 0;
                            index_index < mesh->index_count;
                            index_index += 3)
                    {
                        u32 i_0 = mesh->indices[index_index];
                        u32 i_1 = mesh->indices[index_index + 1];
                        u32 i_2 = mesh->indices[index_index + 2];

                        v3 v_0 = mesh->vertices[i_0];
                        v3 v_1 = mesh->vertices[i_1];
                        v3 v_2 = mesh->vertices[i_2];

                        f32 hit_t = ray_intersect_with_triangle(v_0, v_1, v_2, 
                                                                ray_origin, ray_dir);
                        if(hit_t >= hit_t_threshold && 
                            hit_t < min_t)
                        {
                            min_t = hit_t;
                            hit_mat_index = mesh->mat_index;
                        }
                    }
                }
            }
#endif

            if(hit_mat_index)
            {
                *pixel = world->materials[hit_mat_index].diffuse;
                //*pixel = v3_(1, 0, 0);
            }

            pixel++;
        }
        
        row += output_width;
    }
}
#endif

// TODO(joon) Diffuse reflection, specular reflection is coming soon
// NOTE(joon) N should be normalized!!!
// NOTE(joon) Diffuse means uniform scatter along the normal, so we just use the spherical coordinate,
internal v3
choose_random_direction(RandomSeries *series, v3 N)
{
    // result = (cos(phi) * cos(theta), cos(phi) * sin(theta), sin(phi))
    v3 result = {};

    result = random_spherical_coordinate(series, 0, pi_32/2.0f, 0, 2.0f * pi_32);

    if(compare_equal_f32(N.z, 1.0f))
    {
        // NOTE(joon) because N is already normalized, 
        // this means that N is already (0, 0, 1), so no rotation required
    }
    else if(compare_equal_f32(N.z, 1.0f))
    {
        // N = (0, 0, -1), so just flip
        result.y *= -1.0f;
        result.z *= -1.0f;
    }
    else
    {
        // NOTE(joon) N is nor (0, 0, 1) or (0, 0, -1), so we need to rotate
        // This is basically same as using the inversed rotation matrix along Z that we built,
        // but without meaningless checks
        v3 a = N;
        v3 b = v3_(-N.y, N.x, 0); // cross(a, Z)
        v3 c = cross(a, b);

        // inversed(which means transposed) of the rotation above
        result = result.x * b + result.y * c + result.z * a;
    }
    
    return result;
}

struct SampleRandomLightResult
{
    u32 mat_index;
    v3 p; // random sampled point on the light shape
    v3 normal;

    f32 pdf;
};

internal SampleRandomLightResult
sample_random_lights(TempMemory *light_push_buffer, u32 light_count, RandomSeries *series)
{
    SampleRandomLightResult result = {};

    u32 random_light_index = random_between_u32(series, 0, light_count);
    
    b32 found = false;
    u32 light_index = 0;
    for(u32 consumed = 0;
            !found && consumed < light_push_buffer->used;
            )
    {
        ShapeType *type = (ShapeType *)((u8 *)light_push_buffer->base + consumed);
        consumed += sizeof(*type);
        switch(*type)
        {
            case Shape_Type_Sphere:
            {
                Sphere **sphere_ptr = (Sphere **)((u8 *)light_push_buffer->base + consumed);
                Sphere *sphere = *sphere_ptr;
                consumed += sizeof(sphere);

                if(light_index == random_light_index)
                {
                    result.p = sphere->r * random_spherical_coordinate(series, -pi_32/2.0f, pi_32/2.0f, 0, 2.0f * pi_32) + 
                                      sphere->center;
                    result.normal = normalize(result.p - sphere->center);
                    result.mat_index = sphere->mat_index;

                    // NOTE(joon) pdf if 1/(area of the sphere * light count)
                    result.pdf = 1.0f/(4.0f * pi_32 * sphere->r * sphere->r * light_count);

                    found = true;
                }
            }break;

            case Shape_Type_Cylinder:
            {
                Cylinder **cylinder_ptr = (Cylinder **)((u8 *)light_push_buffer->base + consumed);
                Cylinder *cylinder = *cylinder_ptr;
                consumed += sizeof(*cylinder_ptr);
            }break;

            case Shape_Type_Mesh:
            {
                Mesh **mesh_ptr = (Mesh **)((u8 *)light_push_buffer->base + consumed);
                Mesh *sphere = *mesh_ptr;
                consumed += sizeof(*mesh_ptr);

                // TODO(joon) support for random sampling the mesh lights?? 
            }break;

            default:
            {
                // NOTE(joon) I don't support other shapes for the lights other than the sphere
                invalid_code_path;
            }break;
        }

        light_index++;
    }

    return result;
}

internal void
push_node_to_queue(BVHQueue *queue, BVHOctreeNode *node)
{
    BVHOctreeNode **n = (BVHOctreeNode **)(queue->base + queue->used);
    *n = node;
    queue->used += sizeof(BVHOctreeNode *);

    assert(queue->used <= queue->size);
}

struct RaycastBVHResult
{
    f32 hit_t;
    v3 hit_normal;
    u32 hit_mat_index;

    u32 test_shape_count;

    b32 inner_hit;
};

internal void
raycast_bvh(BVHQueue *queue, v3 ray_origin, v3 ray_dir, RaycastBVHResult *result)
{
    result->hit_t = Flt_Max;

    while(queue->read_cursor < queue->used)
    {
        BVHOctreeNode **node_ptr = (BVHOctreeNode **)(queue->base + queue->read_cursor);
        BVHOctreeNode *node = *node_ptr;

        // NOTE(joon) not doing any aabb intersection testing here,
        // as node will be pushed only if they pass the intersection test

        for(u32 consumed = 0;
                consumed < node->push_buffer.used;
           )
        {
            BVHShapeHeader *header = (BVHShapeHeader *)((u8 *)node->push_buffer.base + consumed);
            consumed += sizeof(*header);

            switch(header->type)
            {
#if 1
                case Shape_Type_Sphere:
                {
                    Sphere *sphere = (Sphere *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*sphere);

                    IntersectionTestResult intersection = ray_intersect_with_sphere(sphere->center, sphere->r, ray_origin, ray_dir);
                    if(intersection.hit_t >= Hit_t_Threshold && intersection.hit_t < result->hit_t)
                    {
                        result->hit_t = intersection.hit_t;
                        result->hit_normal = intersection.hit_normal;
                        result->hit_mat_index = sphere->mat_index;
                        result->inner_hit = intersection.inner_hit;
                    }

                    result->test_shape_count++;
                }break;

                case Shape_Type_AAB:
                {
                    AAB *aab = (AAB *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*aab);

                    IntersectionTestResult intersection = ray_intersect_with_aab(aab->min, aab->max, ray_origin, ray_dir);
                    if(intersection.hit_t >= Hit_t_Threshold && intersection.hit_t < result->hit_t)
                    {
                        result->hit_t = intersection.hit_t;
                        result->hit_normal = intersection.hit_normal;
                        result->hit_mat_index = aab->mat_index;
                        result->inner_hit = intersection.inner_hit;
                    }
                    result->test_shape_count++;
                }break;

                case Shape_Type_Cylinder:
                {
                    Cylinder *cylinder = (Cylinder *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*cylinder);

                    IntersectionTestResult intersection = ray_intersect_with_cylinder(cylinder->base, cylinder->axis, cylinder->r, ray_origin, ray_dir);
                    if(intersection.hit_t >= Hit_t_Threshold && intersection.hit_t < result->hit_t)
                    {
                        result->hit_t = intersection.hit_t;
                        result->hit_normal = intersection.hit_normal;
                        result->hit_mat_index = cylinder->mat_index;
                        result->inner_hit = intersection.inner_hit;
                    }
                    result->test_shape_count++;
                }break;
#endif

                case Shape_Type_Triangle:
                {
                    Triangle *triangle = (Triangle *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*triangle);

                    v3 v_0 = triangle->mesh->vertices[triangle->i_0];
                    v3 v_1 = triangle->mesh->vertices[triangle->i_1];
                    v3 v_2 = triangle->mesh->vertices[triangle->i_2];

                    IntersectionTestResult intersection = ray_intersect_with_triangle(v_0, v_1, v_2, ray_origin, ray_dir);

                    if(intersection.hit_t >= Hit_t_Threshold && intersection.hit_t < result->hit_t)
                    {
                        result->hit_t = intersection.hit_t;
                        result->hit_normal = intersection.hit_normal;
                        result->hit_mat_index = triangle->mesh->mat_index;
                        result->inner_hit = intersection.inner_hit;
                    }
                    result->test_shape_count++;
                }break;

                case Shape_Type_CSG:
                {
                    CSG *csg = (CSG *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*csg);

#if 0
                    IntersectionTestResult *total_intersection = 0;

                    IntersectionTestResult aab_intersection = 
                        ray_intersect_with_aab(csg->aab.min, csg->aab.max, ray_origin, ray_dir);

                    IntersectionTestResult sphere_intersection = 
                        ray_intersect_with_sphere(csg->sphere.center, csg->sphere.r, ray_origin,  ray_dir);

                    if(aab_intersection.hit_t < Hit_t_Threshold)
                    {
                        total_intersection = &sphere_intersection;
                    }
                    else if(sphere_intersection.hit_t < Hit_t_Threshold)
                    {
                        total_intersection = &aab_intersection;
                    }
                    else
                    {
                        total_intersection = &aab_intersection;
                        ray_intersect_with_sphere(csg->sphere.center, csg->sphere.r, ray_origin, ray_dir)

                        // NOTE(joon) both intersections are valid
                        if(aab_intersection.hit_t > sphere_intersection.hit_t)
                        {
                            total_intersection =  &aab_intersection;
                        }
                        else
                        {
                            total_intersection = &sphere_intersection;
                        }
                    }

                    if(total_intersection)
                    {
                        if(total_intersection->hit_t >= Hit_t_Threshold && total_intersection->hit_t < result->hit_t)
                        {
                            result->hit_t = total_intersection->hit_t;
                            result->hit_normal = total_intersection->hit_normal;
                            result->hit_mat_index = csg->mat_index;
                            result->inner_hit = total_intersection->inner_hit;
                        }
                    }
#endif
                }break;

                default:
                {
                    invalid_code_path;
                }
            }// switch
        }//for

        if(node->first_child)
        {
            // NOTE(joon) check if any of the child intersects with the ray,
            // and add them to the queue
            for(u32 child_index = 0;
                    child_index < 8;
                    ++child_index)
            {
                BVHOctreeNode *child = node->first_child + child_index;

                // NOTE(joon) do the testing only if the child a leaf node with shapes inside it
                // or the child is a parent node
                if((child->is_leaf && child->push_buffer.used) || child->first_child)
                {
                    b32 should_add = false;
#if 1
                    if(in_rect(ray_origin, child->aabb_min, child->aabb_max))
                    {
                        should_add = true;
                    }
                    else
#endif
                    {
                        IntersectionTestResult intersection = ray_intersect_with_aab(child->aabb_min, child->aabb_max, ray_origin, ray_dir);
                        if(intersection.hit_t >= Hit_t_Threshold && intersection.hit_t < result->hit_t)
                        {
                            should_add = true;
                        }
                    }

                    if(should_add)
                    {
                        push_node_to_queue(queue, child);
                    }
                }
            }
        }

        queue->read_cursor += sizeof(BVHOctreeNode *);
    }

    result->hit_normal = normalize(result->hit_normal);

    // NOTE(joon) clear the queue, this should always happen here!!!
    queue->read_cursor = 0;
    queue->used = 0;
}

// NOTE(joon) d is a dot()
internal v3
fresnel(v3 Ks, f32 l_dot_h)
{
    v3 result = {};
    result = Ks + (1 - powf(1.0f - absolute(l_dot_h), 5.0f)) * (v3_(1, 1, 1) - Ks);
    return result;
}

// NOTE(joon) This is the D part of the brdf equation
internal f32
ggx_distribution(v3 N, v3 H, f32 roughness)
{
    f32 result = 0.0f;

    f32 n_dot_h = dot(N, H);

    if(n_dot_h > 1.0f || n_dot_h < -1.0f)
    {
        int a = 1;
        f32 n_length = length(N);
        f32 H_length = length(H);
    }

    if(n_dot_h > 0.0f)
    {
        f32 roughness_square = square(roughness);

        // NOTE(joon) math behind this is that tan = sin/cos,
        // cos = dot, sin = sqrt(1 - cos^2)
        f32 tan_theta = sqrt(1.0f - square(n_dot_h)) / n_dot_h;

        f32 nom = roughness_square;
        f32 denom = pi_32 * powf(n_dot_h, 4.0f) * square(roughness_square + square(tan_theta));
        if(!compare_equal_f32(denom, 0.0f))
        {
            result = nom / denom;
        }
    }

    return result;
}

// NOTE(joon) w can be either wo or wi
internal f32
geometry(v3 w, v3 N, v3 m, f32 roughness)
{
    f32 result = 0.0f;

    f32 w_dot_n = dot(w, N);
    f32 w_dot_m = dot(w, m);

    if(!compare_equal_f32(w_dot_m, 0.0f) && (w_dot_n / w_dot_m) > 0)
    {
        if(w_dot_m > 1.0f)
        {
            result = 1.0f;
        }
        else
        {
            // NOTE(joon) math behind this is that tan = sin/cos,
            // cos = dot, sin = sqrt(1 - cos^2)
            f32 tan_theta = sqrt(1.0f - square(w_dot_n)) / w_dot_n;
            if(!compare_equal_f32(tan_theta, 0.0f))
            {
                f32 roughness_square = square(roughness);
                f32 n_dot_m =dot(N, m);
                result = 2.0f / (1.0f + sqrt(1 + roughness_square * square(tan_theta)));
            }
        }
    }

    return result;
}

internal f32
get_radicand(v3 m, v3 wo, f32 n)
{
    // NOTE(joon) radicant r=1–η2(1−(ωo⋅m)2).
    return 1 - square(n)*(1-square(dot(wo, m))); 
}

struct Beern
{
    f32 ni;
    f32 no;
    f32 n;
};

// NOTE(joon) beer's law, n = ni/no
internal Beern
get_beer_n(v3 N, v3 wo, f32 ior)
{
    Beern result = {};
    // NOTE(joon) handle dot(N, wo) == 0 as special case?
    if(dot(N, wo) >= 0.0f)
    {
        result.ni = 1.0f;
        result.no = ior;
    }
    else
    {
        result.ni = ior;
        result.no = 1.0f;
    }

    result.n = result.ni / result.no;

    return result;
}

// NOTE(joon) distance is used to attenuate the light 
internal v3
eval_scattering(v3 N, v3 wi, v3 wo, v3 Kd, v3 Ks, v3 Kt, f32 ior, f32 roughness, f32 distance)
{
    v3 Ed = Kd / pi_32;
    
    v3 H = sign(dot(wi, N)) * normalize((wo + wi));
    f32 wi_dot_h = dot(wi, H);
    
    v3 Es = v3_(0, 0, 0);

    f32 wi_dot_n = dot(wi, N);
    f32 wo_dot_n = dot(wo, N);

    if(wi_dot_h > 0.0f && length_square(Ks) > 0.0f)
    {
        // NOTE(joon) fresnel already contains the Ks(specular term)
        v3 F = fresnel(Ks, wi_dot_h);
        f32 D = ggx_distribution(N, H, roughness);
        f32 G = geometry(wi, N, H, roughness) * geometry(wo, N, H, roughness);
        Es = ((D*G) / (4.0f * absolute(wi_dot_n) * absolute(wo_dot_n))) * F;
    }

    v3 Et = v3_(0, 0, 0);
    if(length_square(Kt) > 0.0f)
    {
        v3 At = v3_(1, 1, 1);
        if(wo_dot_n < 0)
        {
            At.r = powf(euler_number, distance*logf(Kt.r));
            At.g = powf(euler_number, distance*logf(Kt.g));
            At.b = powf(euler_number, distance*logf(Kt.b));
        }
        Beern beer_n = get_beer_n(N, wo, ior);

        v3 m = normalize(-(beer_n.ni*wi + beer_n.no*wo));
        f32 r = get_radicand(m, wo, beer_n.n);

        if(r < 0.0f)
        {
            if(length_square(Ks) > 0.0f)
            {
                // NOTE(joon) total internal reflection
                Et = hadamard(At, Es);
            }
        }
        else
        {
            f32 wi_dot_m = dot(wi, m);
            f32 wo_dot_m = dot(wo, m);

            v3 F = v3_(1, 1, 1) - fresnel(Ks, wi_dot_m);
            f32 D = ggx_distribution(N, m, roughness);
            f32 G = geometry(wi, N, m, roughness) * geometry(wo, N, m, roughness);

            f32 denom = (absolute(wi_dot_n)*absolute(wo_dot_n)*square(beer_n.ni*wi_dot_m + beer_n.no*wo_dot_m));
            if(!compare_equal_f32(denom, 0.0f))
            {
                v3 nom = (D * G * absolute(wi_dot_m)*absolute(wo_dot_m)*square(beer_n.no)) * F;
                Et = hadamard(At, nom / denom);
            }
        }
    }

    v3 result = absolute(wi_dot_n) * (Ed + Es + Et);
    if(is_nan(result))
    {
        int br = 0;
    }
    return result;
}

internal f32
pdf_brdf(RandomSeries *series, v3 N, v3 wi, v3 wo, f32 roughness, v3 Kd, v3 Ks, v3 Kt, f32 ior)
{
    f32 Kd_length = length(Kd);
    f32 Ks_length = length(Ks);
    f32 Kt_length = length(Kt);

    f32 s = Kd_length + Ks_length + Kt_length;

    f32 pd_c = Kd_length / s;
    f32 ps_c = Ks_length / s;
    f32 pt_c = Kt_length / s;

    //NOTE(joon) diffuse term : |N⋅ωi| /π
    f32 pd = absolute(dot(wi, N)) / pi_32;

    v3 H = sign(dot(N, wi)) * normalize((wo + wi));
    f32 n_dot_h = dot(N, H);
    f32 wi_dot_h = dot(wi, H);
    f32 wo_dot_h = dot(wo, H);

    //NOTE(joon) specular term
    // TODO(joon) This is not properly tested!
    f32 ps = 0.0f;
    if(ps_c > 0.0f)
    {
        f32 denom = (4.0f * absolute(wi_dot_h));
        if(!compare_equal_f32(denom, 0.0f))
        {
            f32 D = ggx_distribution(N, H, roughness);
            ps = D * absolute(n_dot_h) / denom;
        }
    }

    Beern beer_n = get_beer_n(N, wo, ior);
    v3 m = normalize(-(beer_n.ni*wi + beer_n.no*wo));
    f32 r = get_radicand(m, wo, beer_n.n); 

    // NOTE(joon) transmission
    f32 pt = ps;
    if(pt_c > 0.0f && r >= 0.0f)
    {
        // TODO(joon) Double check whether we should use H or m here
        f32 n_dot_m = dot(N, m);
        f32 wi_dot_m = dot(wi, m);
        f32 wo_dot_m = dot(wo, m);

        f32 denom = square(beer_n.no*wo_dot_m + beer_n.no*wo_dot_m);
        if(!compare_equal_f32(denom, 0.0f))
        {
            f32 D = ggx_distribution(N, m, roughness);
            pt = D * absolute(n_dot_m) * square(beer_n.no) * absolute(wi_dot_m) / denom;
        }
    }

    return pd_c*pd + ps_c*ps + pt_c*pt;
}

internal v3
sample_lobe(v3 N, f32 c, f32 phi)
{
    // TODO(joon) Make sure that N is already normalized, and remove this ?
    N = normalize(N);
    v3 result = {};
    f32 s = sqrtf(1.0f - c*c);
    v3 K = v3_(s*cos(phi), s * sin(phi), c);

    if(absolute(N.z - 1.0f) < 0.0001f)
    {
        result = K;
    }
    else if(absolute(N.z + 1.0f) < 0.0001f)
    {
        result = v3_(K.x, -K.y, -K.z);
    }
    else
    {
        v3 B = normalize(v3_(-N.y, N.x, 0)); // Z x A
        v3 C = cross(N, B);

        result = K.x*B + K.y*C + K.z*N;
    }

    return result;
}

struct SampleBRDFResult
{
    v3 wi;
    b32 is_transmission;
};

// NOTE(joon) returns wi
internal SampleBRDFResult
sample_brdf(RandomSeries *series, v3 N, v3 wo, f32 roughness, v3 Kd, v3 Ks, v3 Kt, f32 ior)
{
    SampleBRDFResult result = {};

    f32 Kd_length = length(Kd);
    f32 Ks_length = length(Ks);
    f32 Kt_length = length(Kt);

    f32 s = Kd_length + Ks_length + Kt_length;

    f32 pd_c = Kd_length / s;
    f32 ps_c = Ks_length / s;
    f32 pt_c = Kt_length / s;

    f32 e0 = random_between_0_1(series);
    f32 e1 = random_between_0_1(series);

    f32 choice = random_between_0_1(series);

    if(choice < pd_c)
    {
        // diffuse
        result.wi = sample_lobe(N, sqrtf(e0), 2.0f * pi_32 * e1);
    }
    else if(choice >= pd_c && choice < pd_c + ps_c)
    {
        // specular
        f32 ggx_cos_theta = cos(atan2f(roughness * sqrtf(e0), sqrtf(1.0f - e0)));
        v3 m = sample_lobe(N, ggx_cos_theta, 2.0f * pi_32 * e1);

        result.wi = 2.0f * absolute(dot(wo, m)) * m - wo;
    }
    else
    {
        // transmission

        // NOTE(joon) this part lis same as specular 
        f32 ggx_cos_theta = cos(atan2f(roughness * sqrtf(e0), sqrtf(1.0f - e0)));
        v3 m = sample_lobe(N, ggx_cos_theta, 2.0f * pi_32 * e1);

        Beern beer_n = get_beer_n(N, wo, ior);

        f32 r = get_radicand(m, wo, beer_n.n); 

        if(r < 0.0f)
        {
            // NOTE(joon) Again, same as the specular
            result.wi = 2.0f * absolute(dot(wo, m)) * m - wo;
        }
        else
        {
            // TODO(joon) make sure that I got the equation right
            result.wi = (beer_n.n * dot(wo, m) - sign(dot(wo, N)) * sqrtf(r)) * m - beer_n.n*wo;
            result.is_transmission = true;
        }
    }

    result.wi = normalize(result.wi);
    
    return result;
}

// NOTE(joon) Push the top most node to the queue,
// this has a same functionality of raytracing the whole scene.
internal RaycastBVHResult
raycast_top_most_node(BVHQueue *queue, BVHOctreeNode *top_most_node, u64 *test_shape_count, v3 ray_origin, v3 ray_dir)
{
    RaycastBVHResult result = {};

    push_node_to_queue(queue, top_most_node);
    raycast_bvh(queue, ray_origin, ray_dir, &result);

    *test_shape_count += result.test_shape_count;

    return result;
}

internal u64
tiled_raytrace_bvh(World *world, Camera *camera, BVHOctreeNode *top_most_node,
                v3 *output_buffer, i32 output_width, i32 output_height, 
                i32 tile_min_x, i32 tile_min_y,
                i32 tile_one_past_max_x, i32 tile_one_past_max_y, 
                RandomSeries *series, u32 ray_per_pixel_count, f32 russian_roulette_value)
{
    u64 test_shape_count = 0;
    // TODO(joon) : priority sorting!
    BVHQueue queue = {};

    // TODO(joon) malloc can be pontentially bad. Replace this with the memory arena 
    queue.size = megabytes(4);
    queue.base = (u8 *)malloc(queue.size);
    queue.used = 0;

    f32 roughness = 0.01f;
    // TODO(joon) check whether this is working fine with the refraction
    f32 dont_get_too_close_epsilon = 0.0001f;

    f32 focal_length = length(camera->p - v3_(0, 0, 0.2f));
    f32 aperture_radius = 0.1f;

    v3 *row = output_buffer + tile_min_y * output_width + tile_min_x;
    for(i32 y = tile_min_y;
            y < tile_one_past_max_y;
            ++y)
    {
        v3 *pixel = row;
        for(i32 x = tile_min_x;
                x < tile_one_past_max_x;
                ++x)
        {
            v3 color = {};

#if 1
            // NOTE(joon) get the point to the center of the pixel from the eye
            f32 pixel_x = (2.0f * x / (f32)output_width) - 1.0f;
            f32 pixel_y = (2.0f * y / (f32)output_height) - 1.0f;

            // NOTE(joon) math here is that camera_z is facing the opposite direction of where the camera is facing,
            // and x and y are responsible for where the pixel is in camera plane
            v3 camera_to_pixel = normalize(pixel_x*camera->x_axis + pixel_y*camera->y_axis - camera->z_axis);
            v3 focal_point = camera->p + focal_length*camera_to_pixel;
                
                // NOTE(joon) This is the first ray that goes through the screen pixel
                // v3 initial_ray_dir = normalize(focal_point - random_p_in_aperture);
#endif

            for(u32 ray_index = 0;
                    ray_index < ray_per_pixel_count;
                    ++ray_index)
            {

                f32 random_rad = random_between(series, 0.0f, 2 * pi_32);
                v3 random_p_in_aperture = camera->p + aperture_radius*cosf(random_rad) * camera->x_axis +
                                          aperture_radius*sinf(random_rad) * camera->y_axis - 0.1f * camera->z_axis;
                
                // NOTE(joon) This is the first ray that goes through the screen pixel
                v3 initial_ray_dir = normalize(focal_point - random_p_in_aperture);

                // NOTE(joon) outgoing light direction
                v3 wo = -normalize(initial_ray_dir);

                v3 ray_origin = random_p_in_aperture;
                b32 ray_alive = true;
                v3 previous_ray_dir = {};
                v3 hit_normal = v3_(0, 0, 0);
                v3 weight = v3_(1, 1, 1);
                Material *hit_mat = 0;

                RaycastBVHResult initial_raycast_result = 
                        raycast_top_most_node(&queue, top_most_node, &test_shape_count, ray_origin, initial_ray_dir);
                if(initial_raycast_result.hit_mat_index)
                {
                    Material *mat = world->materials + initial_raycast_result.hit_mat_index;
                    if(mat->is_light)
                    {
                        // NOTE(joon) This is an initial ray, we don't need to multiply by the weight
                        color += mat->emit_color; 
                        ray_alive = false;
                    }
                    else
                    {
                        ray_origin = ray_origin + (initial_raycast_result.hit_t - dont_get_too_close_epsilon) * initial_ray_dir;
                        hit_normal = initial_raycast_result.hit_normal;
                        hit_mat = mat;
                        previous_ray_dir = initial_ray_dir;

                        if(length_square(hit_mat->diffuse) > 0.0f)
                        {
                            weight = hadamard(weight, hit_mat->diffuse);
                        }

                        if(hit_mat->transmission.x > 0.0f)
                        {
                            int a = 1;
                        }
                    }
                }

#if 1
                while(ray_alive && random_between_0_1(series) < russian_roulette_value)
                {
                    // NOTE(joon): Explicit light connection
                    SampleRandomLightResult sampled_light = sample_random_lights(&world->light_push_buffer, world->light_count, series);
                    v3 explicit_wi = normalize(sampled_light.p - ray_origin);
#if 0
                    f32 geometry_factor = absolute(dot(hit_normal, -explicit_wi) * dot(sampled_light.normal, -explicit_wi) / square(dot(-explicit_wi, -explicit_wi)));
                    f32 probability = sampled_light.pdf / geometry_factor;
                    //f32 q = pdf_brdf(series, hit_normal, wi, wo, roughness, hit_mat->diffuse, hit_mat->specular.xyz, hit_mat->transmission, hit_mat->ior) * russian_roulette_value; // Prob the explicit light could be chosen implicitly w =p2/(p2+q2)
                    //f32 w_miss = square(probability) / (square(q) * square(probability));

                    if(probability > 0.00001f)
                    {
                        RaycastBVHResult shadow_raycast_result = 
                            raycast_top_most_node(&queue, top_most_node, &test_shape_count, ray_origin, explicit_wi);

                        if(shadow_raycast_result.hit_t >= Hit_t_Threshold)
                        {
                            v3 shadow_ray_hit_p = ray_origin + shadow_raycast_result.hit_t * explicit_wi;
                            if(compare_equal(shadow_ray_hit_p, sampled_light.p))
                            {
                                // NOTE(joon) : There is no obstacle between the hit point & the light p that we sampled, 
                                // thus two points are connected
                                Material *light_mat = world->materials + sampled_light.mat_index;

                                // TODO(joon) something is def wrong with this function. check wi and wo!
                                // TODO(joon) take account of q, too!
                                //v3 f = eval_scattering(hit_normal, wi, wo, hit_mat->diffuse, hit_mat->specular.xyz, roughness);

                                // NOTE(joon) this already take account of the diffuse color of the hit_mat
                                v3 c = 0.2f * hadamard(hit_mat->diffuse, hadamard(weight, light_mat->emit_color));
                                if(!is_nan(c) && !is_infinite(c))
                                {
                                    color += c;
                                }
                                else
                                {
                                    int br = 0;
                                }
                            }
                        }
                    }

                    if(hit_mat->transmission.x > 0.0f)
                    {
                        int b = 0;
                    }
#endif

                    if(hit_mat->transmission.x > 0.0f)
                    {
                        int b = 0;
                    }

                    // NOTE(joon) Implicit light connection 
                    SampleBRDFResult sampled_brdf = sample_brdf(series, hit_normal, wo, roughness, hit_mat->diffuse, hit_mat->specular.xyz, hit_mat->transmission, hit_mat->ior); // NOTE(joon) random direction from P
                    if(is_nan(sampled_brdf.wi))
                    {
                        int b = 0;
                    }
                    v3 wi = sampled_brdf.wi;

#if 1
                    // NOTE(joon) If the ray is going to refract, we move the ray a little bit more so that the ray 
                    // can be inside the object that is going to refract.
                    if(sampled_brdf.is_transmission)
                    {
                        ray_origin += 2.0f * dont_get_too_close_epsilon * previous_ray_dir;
                    }
#endif
#if 1
                    // TODO(joon) double check that the P is next_ray_origin
                    RaycastBVHResult extended_path_raycast_result = 
                        raycast_top_most_node(&queue, top_most_node, &test_shape_count, ray_origin, wi);

                    if(extended_path_raycast_result.hit_mat_index)
                    {
                        Material *extended_path_mat = world->materials + extended_path_raycast_result.hit_mat_index;
                        if(extended_path_mat->is_light)
                        {
                            // TODO(joon) take account of q, too!
                            v3 c = hadamard(weight, extended_path_mat->emit_color);
                            if(!is_nan(c) && !is_infinite(c))
                            {
                                color += c;
                            }
                            else 
                            {
                                int br = 0;
                            }
                            ray_alive = false;
                        }
                        else
                        {
                            f32 a = pdf_brdf(series, extended_path_raycast_result.hit_normal, wi, wo, roughness, extended_path_mat->diffuse, extended_path_mat->specular.xyz, extended_path_mat->transmission, extended_path_mat->ior) * russian_roulette_value;
                            if(isnan(a))
                            {
                                int b = 0;
                            }

                            f32 p = pdf_brdf(series, extended_path_raycast_result.hit_normal, wi, wo, roughness, extended_path_mat->diffuse, extended_path_mat->specular.xyz, extended_path_mat->transmission, extended_path_mat->ior) * russian_roulette_value;

                            // NOTE(joon) avoid division by 0
                            if(p > 0.000001f)
                            {
                                // NOTE(joon) wi is always unit length
                                f32 distance = extended_path_raycast_result.hit_t;
                                v3 f0 = eval_scattering(extended_path_raycast_result.hit_normal, wi, wo, 
                                                        extended_path_mat->diffuse, 
                                                        extended_path_mat->specular.xyz, 
                                                        extended_path_mat->transmission, extended_path_mat->ior, 
                                                        roughness, distance);
                                if(is_nan(f0))
                                {
                                    int break_here = 0;
                                }

                                v3 f = eval_scattering(extended_path_raycast_result.hit_normal, wi, wo, 
                                                        extended_path_mat->diffuse, 
                                                        extended_path_mat->specular.xyz, 
                                                        extended_path_mat->transmission, extended_path_mat->ior, 
                                                        roughness, distance);

                                weight = hadamard(f/p, weight);

                            }
                            else
                            {
                                //ray_alive = false;
                            }

                            ray_origin = ray_origin + (extended_path_raycast_result.hit_t - dont_get_too_close_epsilon) * wi;
                            hit_normal = extended_path_raycast_result.hit_normal;
                            hit_mat = extended_path_mat;
                            previous_ray_dir = wi;
                            wo = -wi;
                        }
                    }
                    else 
                    {
                        ray_alive = false;
                    }
#endif
                }
#endif

            }

            *pixel = color / ray_per_pixel_count;
            if(is_nan(*pixel) || compare_0(*pixel))
            {
                int a = 1;
            }
            pixel++;
        }
        
        row += output_width;
    }

    // NOTE(joon) thread tile region test
#if 0
    row = output_buffer + tile_min_y * output_width + tile_min_x;
    for(i32 y = tile_min_y;
            y < tile_one_past_max_y;
            ++y)
    {
        v3 *pixel = row;
        for(i32 x = tile_min_x;
                x < tile_one_past_max_x;
                ++x)
        {
            if(x == 0 || x == tile_one_past_max_x-1 ||
                y == 0 || y == tile_one_past_max_y-1)
            {
                *pixel = v3_(1, 0, 0);
            }
            pixel++;
        }

        row += output_width;
    }
#endif

    free(queue.base);

    return test_shape_count;
}

// NOTE(joon) This is not a great name...
struct bvh_octree_node_child_info
{
    u8 bitmask;// how should we get the next child?
    v3 center;
    v3 half_dim;
};

internal bvh_octree_node_child_info
get_bvh_octree_node_child_info(v3 node_center, v3 node_half_dim, v3 shape_center)
{
    bvh_octree_node_child_info result = {};
    result.bitmask = 0b11111111;
    result.center = node_center;
    result.half_dim = 0.5f*node_half_dim;

    if(shape_center.x >= node_center.x)
    {
        result.bitmask = result.bitmask & BVH_Node_Pos_X_Mask;
        result.center.x += result.half_dim.x;
    }
    else
    {
        result.bitmask = result.bitmask & BVH_Node_Neg_X_Mask;
        result.center.x -= result.half_dim.x;
    }

    if(shape_center.y >= node_center.y)
    {
        result.bitmask = result.bitmask & BVH_Node_Pos_Y_Mask;
        result.center.y += result.half_dim.y;
    }
    else
    {
        result.bitmask = result.bitmask & BVH_Node_Neg_Y_Mask;
        result.center.y -= result.half_dim.y;
    }

    if(shape_center.z >= node_center.z)
    {
        result.bitmask = result.bitmask & BVH_Node_Pos_Z_Mask;
        result.center.z += result.half_dim.z;
    }
    else
    {
        result.bitmask = result.bitmask & BVH_Node_Neg_Z_Mask;
        result.center.z -= result.half_dim.z;
    }
    
    // TODO(joon) If we go too deep, this will fire because of the floating point precision
    //assert(in_rect(shape_center, result.center - result.half_dim, result.center + result.half_dim));
    assert(result.bitmask != 0b11111111);

    return result;
}

internal void
push_shape(TempMemory *push_buffer, void *shape, ShapeType type)
{
    if(push_buffer->base == 0)
    {
        push_buffer->total_size = kilobytes(256);
        push_buffer->base = (u8 *)malloc(push_buffer->total_size);
    }
    else
    {
        if(push_buffer->used >= 0.9f*push_buffer->total_size)
        {
            push_buffer->total_size *= 2;

            u8 *new_memory = (u8 *)malloc(push_buffer->total_size);
            memcpy(new_memory, push_buffer->base, push_buffer->used);
            free(push_buffer->base);

            push_buffer->base = new_memory;
        }
    }

#if 0
    switch(type)
    {
        case Shape_Type_Sphere:
        {
            BVHShapeHeader *header = push_struct(push_buffer, BVHShapeHeader);
            header->type = Shape_Type_Sphere;
            Sphere **sphere_ptr = (Sphere **)push_size(push_buffer, sizeof(Sphere *));
            *sphere_ptr = (Sphere *)shape;
        }break;
        case Shape_Type_Cylinder:
        {
            BVHShapeHeader *header = push_struct(push_buffer, BVHShapeHeader);
            header->type = Shape_Type_Cylinder;
            Cylinder **cylinder_ptr = (Cylinder **)push_size(push_buffer, sizeof(Cylinder *));
            *cylinder_ptr = (Cylinder *)shape;
        }break;
        case Shape_Type_AAB:
        {
            BVHShapeHeader *header = push_struct(push_buffer, BVHShapeHeader);
            header->type = Shape_Type_AAB;
            AAB **aab_ptr = (AAB **)push_size(push_buffer, sizeof(AAB *));
            *aab_ptr = (AAB *)shape;
        }break;
        case Shape_Type_Mesh:
        {
            BVHShapeHeader *header = push_struct(push_buffer, BVHShapeHeader);
            header->type = Shape_Type_Mesh;
            Mesh **mesh_ptr = (Mesh **)push_size(push_buffer, sizeof(Mesh *));
            *mesh_ptr = (Mesh *)shape;
        }break;
    }
#endif

#if 1
    switch(type)
    {
        case Shape_Type_Sphere:
        {
            BVHShapeHeader *header = push_struct(push_buffer, BVHShapeHeader);
            header->type = Shape_Type_Sphere;
            Sphere *sphere = (Sphere *)push_size(push_buffer, sizeof(Sphere));
            *sphere = *(Sphere *)shape;
        }break;

        case Shape_Type_Cylinder:
        {
            BVHShapeHeader *header = push_struct(push_buffer, BVHShapeHeader);
            header->type = Shape_Type_Cylinder;
            Cylinder *cylinder = (Cylinder *)push_size(push_buffer, sizeof(Cylinder));
            *cylinder = *(Cylinder *)shape;
        }break;

        case Shape_Type_AAB:
        {
            BVHShapeHeader *header = push_struct(push_buffer, BVHShapeHeader);
            header->type = Shape_Type_AAB;
            AAB *aab = (AAB *)push_size(push_buffer, sizeof(AAB));
            *aab = *(AAB *)shape;
        }break;

        case Shape_Type_Triangle:
        {
            BVHShapeHeader *header = push_struct(push_buffer, BVHShapeHeader);
            header->type = Shape_Type_Triangle;
            Triangle *triangle = (Triangle *)push_size(push_buffer, sizeof(Triangle));
            *triangle = *(Triangle *)shape;
        }break;

        case Shape_Type_CSG:
        {
            BVHShapeHeader *header = push_struct(push_buffer, BVHShapeHeader);
            header->type = Shape_Type_CSG;
            CSG *csg = (CSG *)push_size(push_buffer, sizeof(CSG));
            *csg = *(CSG *)shape;
        }break;

        default:
        {
            invalid_code_path;
        }
    }
#endif
}

internal bvh_octree_node_child_info
get_bvh_octree_node_child_info(v3 node_center, v3 node_half_dim, u8 bitmask)
{
    bvh_octree_node_child_info result = {};
    result.bitmask = bitmask;
    result.center = node_center;
    result.half_dim = 0.5f*node_half_dim;

    if(bitmask & BVH_Node_Pos_X_Mask)
    {
        result.center.x += result.half_dim.x;
    }
    else
    {
        result.center.x -= result.half_dim.x;
    }

    if(bitmask & BVH_Node_Pos_Y_Mask)
    {
        result.center.y += result.half_dim.y;
    }
    else
    {
        result.center.y -= result.half_dim.y;
    }

    if(bitmask & BVH_Node_Pos_Z_Mask)
    {
        result.center.z += result.half_dim.z;
    }
    else
    {
        result.center.z -= result.half_dim.z;
    }

    return result;
}

struct ShapeAABB
{
    v3 center;
    v3 half_dim;
};

internal ShapeAABB
get_shape_aabb(void *shape, ShapeType type)
{
    ShapeAABB result = {};

    switch(type)
    {
        case Shape_Type_Sphere:
        {
            Sphere *sphere = (Sphere *)shape;
            result.center = sphere->center;
            result.half_dim = sphere->r * v3_(1, 1, 1);
        }break;

        case Shape_Type_Cylinder:
        {
            Cylinder *cylinder = (Cylinder *)shape;

            v3 cylinder_another_base = cylinder->base + cylinder->axis;
            v3 e = cylinder->r * (v3_(1, 1, 1) - sqrt(hadamard(cylinder->axis, cylinder->axis) / dot(cylinder->axis, cylinder->axis)));
            v3 min = gather_min_elements(cylinder->base - e, cylinder_another_base - e);
            v3 max = gather_max_elements(cylinder->base + e, cylinder_another_base + e);

            result.center = 0.5f * (min + max);
            result.half_dim = max - result.center;
        }break;

        case Shape_Type_AAB:
        {
            AAB *aab = (AAB *)shape;

            result.center = 0.5f * (aab->min + aab->max);
            result.half_dim = aab->max - result.center;
        }break;

        case Shape_Type_Mesh:
        {
            Mesh *mesh = (Mesh *)shape;
            result.center = 0.5f * (mesh->aabb_min + mesh->aabb_max);
            result.half_dim = mesh->aabb_max - result.center;
        }break;

        case Shape_Type_Triangle:
        {
            Triangle *triangle = (Triangle *)shape;

            v3 v_0 = triangle->mesh->vertices[triangle->i_0];
            v3 v_1 = triangle->mesh->vertices[triangle->i_1];
            v3 v_2 = triangle->mesh->vertices[triangle->i_2];

            v3 aabb_min = gather_min_elements(gather_min_elements(v_0, v_1), v_2); 
            v3 aabb_max = gather_max_elements(gather_max_elements(v_0, v_1), v_2); 

            result.center = 0.5f*(aabb_min + aabb_max);
            result.half_dim = aabb_max - result.center;
        }break;

        case Shape_Type_CSG:
        {
            CSG *csg = (CSG *)shape;
            result.center = 0.5f*(csg->aabb_min + csg->aabb_max);
            result.half_dim = csg->aabb_max - result.center;
        }break;

        default:
        {
            invalid_code_path;
        }
    }

    return result;
}

internal void
initialize_all_childs(BVHOctreeNode *node)
{
    assert(node->first_child);

    for(u32 child_index = 0;
            child_index < 8;
            ++child_index)
    {
        BVHOctreeNode *child = node->first_child + child_index;
        zero(child);
        child->is_leaf = true;
        child->aabb_min = v3_(Flt_Max, Flt_Max, Flt_Max);
        child->aabb_max = v3_(Flt_Min, Flt_Min, Flt_Min);
    }
}

internal void
update_aabb_min_max(v3 *dest_min, v3 *dest_max, void *shape, ShapeType type)
{
    ShapeAABB aabb = get_shape_aabb(shape, type);
    
    v3 aabb_min = aabb.center - aabb.half_dim;
    v3 aabb_max = aabb.center + aabb.half_dim;

    *dest_min = gather_min_elements(*dest_min, aabb_min);
    *dest_max = gather_max_elements(*dest_max, aabb_max);

    assert(dest_min->x <= dest_max->x && dest_min->y <= dest_max->y && dest_min->z <= dest_max->z);
}

internal b32
test_aabb_aabb(v3 a_center, v3 a_half_dim,
                v3 b_center, v3 b_half_dim)
{
    b32 result = false;

    v3 half_dim_sum = a_half_dim + b_half_dim; 
    v3 rel_a_center = a_center - b_center;

    // TODO(joon) early return is faster?
    if((rel_a_center.x >= -half_dim_sum.x && rel_a_center.x < half_dim_sum.x) && 
       (rel_a_center.y >= -half_dim_sum.y && rel_a_center.y < half_dim_sum.y) && 
       (rel_a_center.z >= -half_dim_sum.z && rel_a_center.z < half_dim_sum.z))
    {
        result = true;
    }

    return result;
}

internal void
push_shape_inside_node(MemoryArena *node_arena, 
                    BVHOctreeNode *node, 
                    v3 node_center, v3 node_half_dim,
                    u32 current_depth, u32 desired_depth, 
                    void *shape, ShapeType type)
{
    update_aabb_min_max(&node->aabb_min, &node->aabb_max, shape, type);

    if(current_depth < desired_depth)
    {
        if(node->first_child)
        {
            // NOTE(joon) move to next child
            ShapeAABB shape_aabb = get_shape_aabb(shape, type);

            bvh_octree_node_child_info child_info = 
                get_bvh_octree_node_child_info(node_center, node_half_dim, shape_aabb.center);

            u32 child_index = find_least_significant_bit(child_info.bitmask);

            BVHOctreeNode *child = node->first_child + child_index;

            push_shape_inside_node(node_arena,  
                    child, 
                    child_info.center, child_info.half_dim,
                    current_depth + 1, desired_depth, 
                    shape, type);
        }
        else
        {
            if(node->push_buffer.used == 0)
            {
                // NOTE(joon) : node is empty and there is no child node, we can just push it here
                push_shape(&node->push_buffer, shape, type);
            }
            else
            {
                // NOTE(joon) : node was not empty, relocate both this shape & whatever shapes that were inside this node
                // to the triangle->child nodes
                node->first_child = push_array(node_arena, BVHOctreeNode, 8);
                initialize_all_childs(node);

                for(u32 consumed = 0;
                        consumed < node->push_buffer.used;
                        )
                {
                    BVHShapeHeader *header = (BVHShapeHeader *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*header);

                    void *existing_shape = 0;
                    switch(header->type)
                    {
                        case Shape_Type_Sphere:
                        {
                            Sphere *sphere = (Sphere *)((u8 *)node->push_buffer.base + consumed);
                            consumed += sizeof(*sphere);

                            existing_shape = sphere;
                        }break;

                        case Shape_Type_AAB:
                        {
                            AAB *aab = (AAB *)((u8 *)node->push_buffer.base + consumed);
                            consumed += sizeof(*aab);

                            existing_shape = aab;
                        }break;

                        case Shape_Type_Cylinder:
                        {
                            Cylinder *cylinder = (Cylinder *)((u8 *)node->push_buffer.base + consumed);
                            consumed += sizeof(*cylinder);

                            existing_shape = cylinder;
                        }break;

                        case Shape_Type_Mesh:
                        {
                            Mesh *mesh = (Mesh *)((u8 *)node->push_buffer.base + consumed);
                            consumed += sizeof(*mesh);

                            existing_shape = mesh;
                        }break;

                        case Shape_Type_Triangle:
                        {
                            Triangle *triangle = (Triangle *)((u8 *)node->push_buffer.base + consumed);
                            consumed += sizeof(*triangle);

                            existing_shape = triangle;
                        }break;

                        case Shape_Type_CSG:
                        {
                            CSG *csg = (CSG *)((u8 *)node->push_buffer.base + consumed);
                            consumed += sizeof(*csg);

                            existing_shape = csg;
                        }break;

                        default:
                        {
                            invalid_code_path;
                        }
                    }

                    ShapeAABB shape_aabb = get_shape_aabb(existing_shape, header->type);

                    bvh_octree_node_child_info child_info = 
                        get_bvh_octree_node_child_info(node_center, node_half_dim, shape_aabb.center);

                    u32 child_index = find_least_significant_bit(child_info.bitmask);

                    BVHOctreeNode *child = node->first_child + child_index;

                    push_shape_inside_node(node_arena,  
                            child, 
                            child_info.center, child_info.half_dim,
                            current_depth + 1, desired_depth, 
                            existing_shape, header->type);
                }
                node->push_buffer.total_size = 0; 
                node->push_buffer.used = 0;
                node->push_buffer.base = 0;
                node->is_leaf = false;

                ShapeAABB shape_aabb = get_shape_aabb(shape, type);

                bvh_octree_node_child_info child_info = 
                    get_bvh_octree_node_child_info(node_center, node_half_dim, shape_aabb.center);

                u32 child_index = find_least_significant_bit(child_info.bitmask);
                BVHOctreeNode *child = node->first_child + child_index;

                push_shape_inside_node(node_arena,  
                        child, 
                        child_info.center, child_info.half_dim,
                        current_depth + 1, desired_depth, 
                        shape, type);
            }
        }
    }
    else
    {
        assert(node->is_leaf == true);
        // NOTE(joon) we cannot go deeper, we just put the shape here
        push_shape(&node->push_buffer, shape, type);
    }
}

struct ValidateNodesResult
{
    u32 sphere_count; 
    u32 cylinder_count; 
    u32 box_count; 
    u32 mesh_count; 
    u32 csg_count;
    u32 memory_used;
};

internal void
validate_nodes_and_reallocate_shapes(MemoryArena *arena, BVHOctreeNode *node, ValidateNodesResult *result)
{
    assert((node->push_buffer.used == 0 && !node->is_leaf) || node->is_leaf);
    assert((node->is_leaf && !node->first_child) || !node->is_leaf);
    //assert(node->aabb_min.x <= node->aabb_max.x && node->aabb_min.y <= node->aabb_max.y && node->aabb_min.z <= node->aabb_max.z);

    if(node->push_buffer.used)
    {
        TempMemory new_push_buffer = start_temp_memory(arena, node->push_buffer.used);
        memcpy(new_push_buffer.base, node->push_buffer.base, node->push_buffer.used);
        new_push_buffer.used = node->push_buffer.used;
        free(node->push_buffer.base);
        node->push_buffer = new_push_buffer;

        result->memory_used += node->push_buffer.used;

        for(u32 consumed = 0;
                consumed < node->push_buffer.used;
           )
        {
            BVHShapeHeader *header = (BVHShapeHeader *)((u8 *)node->push_buffer.base + consumed);
            consumed += sizeof(*header);

            switch(header->type)
            {
                case Shape_Type_Sphere:
                {
                    Sphere *sphere = (Sphere *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*sphere);
                    result->sphere_count++;
                }break;

                case Shape_Type_AAB:
                {
                    AAB *aab = (AAB *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*aab);
                    result->box_count++;
                }break;

                case Shape_Type_Cylinder:
                {
                    Cylinder *cylinder = (Cylinder *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*cylinder);
                    result->cylinder_count++;
                }break;

                case Shape_Type_Mesh:
                {
                    Mesh *mesh = (Mesh *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*mesh);
                    result->mesh_count++;
                }break;

                case Shape_Type_CSG:
                {
                    CSG *csg = (CSG *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*csg);
                    result->csg_count++;
                }break;

                case Shape_Type_Triangle:
                {
                    Triangle *triangle = (Triangle *)((u8 *)node->push_buffer.base + consumed);
                    consumed += sizeof(*triangle);
                }break;

                default:
                {
                    invalid_code_path;
                }
            }// switch
        }//for
    }
     
    if(node->first_child)
    {
        for(u32 child_index = 0;
                child_index < 8;
                ++child_index)
        {
            BVHOctreeNode *child = node->first_child + child_index;
            validate_nodes_and_reallocate_shapes(arena, child, result);
        }
    }
}










