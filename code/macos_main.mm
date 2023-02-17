#include <mach/mach_time.h> // mach_absolute_time
#include <stdio.h> // printf for debugging purpose
#include <sys/stat.h>
#include <libkern/OSAtomic.h>
#include <pthread.h>
#include <semaphore.h>
#include <Carbon/Carbon.h>
#include <dlfcn.h> // dlsym
#include <arm_neon.h>

#undef internal
#undef assert

#include "types.h"
#include "math.h"
#include "platform.h"
#include "intrinsic.h"
#include "random.h"

// NOTE(joon) unity build
#include "ray.cpp"
#include "parser.cpp"

global dispatch_semaphore_t semaphore;

internal u64 
mach_time_diff_in_nano_seconds(u64 begin, u64 end, r32 nano_seconds_per_tick)
{
    return (u64)(((end - begin)*nano_seconds_per_tick));
}


PLATFORM_GET_FILE_SIZE(macos_get_file_size) 
{
    u64 result = 0;

    int File = open(filename, O_RDONLY);
    struct stat FileStat;
    fstat(File , &FileStat); 
    result = FileStat.st_size;
    close(File);

    return result;
}

PLATFORM_READ_FILE(debug_macos_read_file)
{
    PlatformReadFileResult Result = {};

    int File = open(filename, O_RDONLY);
    int Error = errno;
    if(File >= 0) // NOTE : If the open() succeded, the return value is non-negative value.
    {
        struct stat FileStat;
        fstat(File , &FileStat); 
        off_t fileSize = FileStat.st_size;

        if(fileSize > 0)
        {
            // TODO/Joon : NO MORE OS LEVEL ALLOCATION!
            Result.size = fileSize;
            Result.memory = (u8 *)malloc(Result.size);
            if(read(File, Result.memory, fileSize) == -1)
            {
                free(Result.memory);
                Result.size = 0;
            }
        }

        close(File);
    }

    return Result;
}

PLATFORM_WRITE_ENTIRE_FILE(debug_macos_write_entire_file)
{
    int file = open(file_name, O_WRONLY|O_CREAT|O_TRUNC, S_IRWXU);

    if(file >= 0) 
    {
        if(write(file, memory_to_write, size) == -1)
        {
            // TODO(joon) : LOG here
        }

        close(file);
    }
    else
    {
        // TODO(joon) :LOG
        printf("Failed to create file\n");
    }
}

PLATFORM_FREE_FILE_MEMORY(debug_macos_free_file_memory)
{
    free(memory);
}

// TODO(joon) : It seems like this combines read & write barrier, but make sure
// TODO(joon) : mfence?(DSB)
#define write_barrier OSMemoryBarrier(); asm volatile("": : :"memory");
#define read_barrier OSMemoryBarrier(); asm volatile("": : :"memory");

struct MacOSThread
{
    u32 ID;
    ThreadWorkQueue *queue;

    // TODO(joon): I like the idea of each thread having a random number generator that they can use throughout the whole process
    // though what should happen to the 0th thread(which does not have this structure)?
    //simd_random_series series;
};

// NOTE(joon) : use this to add what thread should do
internal 
THREAD_WORK_CALLBACK(print_string)
{
    char *stringToPrint = (char *)data;
    printf("%s\n", stringToPrint);
}

struct ThreadWorkRaytraceTileData
{
    World *world;
    Camera *camera;

    BVHOctreeNode *top_most_node;

    RandomSeries random_series;

    v3 *output_buffer;
    i32 output_width;
    i32 output_height;

    // NOTE(joon) only for printing out at the end
    u32 tile_x;
    u32 tile_y;

    i32 tile_min_x;
    i32 tile_min_y;
    i32 tile_one_past_max_x;
    i32 tile_one_past_max_y;

    u32 rays_per_pixel_count;
    f32 russian_roulette_value;
};

internal 
THREAD_WORK_CALLBACK(thread_callback_tiled_raytrace)
{
    ThreadWorkRaytraceTileData *input_data = (ThreadWorkRaytraceTileData *)data;
    u64 shape_count = tiled_raytrace_bvh(input_data->world, input_data->camera, input_data->top_most_node, 
                      input_data->output_buffer, 
                    input_data->output_width, input_data->output_height,
                    input_data->tile_min_x, input_data->tile_min_y,
                    input_data->tile_one_past_max_x, input_data->tile_one_past_max_y, 
                    &input_data->random_series, input_data->rays_per_pixel_count, input_data->russian_roulette_value);

    printf("Thread %d has rendered tile(%d, %d) with %llu tested shapes\n", thread_index, input_data->tile_x, input_data->tile_y, shape_count);
    OSAtomicIncrement32Barrier((volatile i32 *)&input_data->world->rendered_tile_count);
}

// NOTE(joon): This is single producer multiple consumer - 
// meaning, it _does not_ provide any thread safety
// For example, if the two threads try to add the work item,
// one item might end up over-writing the other one
internal void
macos_add_thread_work_item(ThreadWorkQueue *queue,
                            thread_work_callback *work_callback,
                            void *data)
{
    assert(data); // TODO(joon) : There might be a work that does not need any data?
    ThreadWorkItem *item = queue->items + queue->add_index;
    item->callback = work_callback;
    item->data = data;
    item->written = true;

    // TODO(joon) prevent compiler reordering!!!
    write_barrier;
    queue->add_index++;
    assert(queue->add_index <= queue->item_count);

    // increment the semaphore value by 1
    dispatch_semaphore_signal(semaphore);
}

internal b32
macos_do_thread_work_item(ThreadWorkQueue *queue, u32 thread_index)
{
    b32 did_work = false;
    if(queue->work_index != queue->add_index)
    {
        int original_work_index = queue->work_index;
        int expected_work_index = original_work_index + 1;

        if(OSAtomicCompareAndSwapIntBarrier(original_work_index, expected_work_index, &queue->work_index))
        {
            ThreadWorkItem *item = queue->items + original_work_index;
            read_barrier;
            item->callback(item->data, thread_index);

            did_work = true;
        }
    }

    return did_work;
}

internal 
PLATFORM_COMPLETE_ALL_THREAD_WORK_QUEUE_ITEMS(macos_complete_all_thread_work_queue_items)
{
    // TODO(joon): If there was a last thread that was working on the item,
    // this does not guarantee that the last work will be finished.
    // Maybe add some flag inside the thread? (sleep / working / ...)
    while(queue->work_index != queue->add_index) 
    {
        macos_do_thread_work_item(queue, 0);
    }
}

internal void*
thread_proc(void *data)
{
    MacOSThread *thread = (MacOSThread *)data;
    while(1)
    {
        if(macos_do_thread_work_item(thread->queue, thread->ID))
        {
        }
        else
        {
            // dispatch semaphore puts the thread into sleep until the semaphore is signaled
            dispatch_semaphore_wait(semaphore, DISPATCH_TIME_FOREVER);
        }
    }

    return 0;
}

internal u32
v3_to_rgbe(v3 color)
{
    u32 result = 0;

    f32 max_component = maximum(maximum(color.r, color.g), color.b);
    i32 e;

    f32 color_threshold = 1e-32f;
    if (max_component >= color_threshold) 
    {
        f32 denom = frexp(max_component, &e) * 255.0f / max_component;
        result = ((round_r32_u32(color.r * denom) << 0) |
                (round_r32_u32(color.g * denom) << 8) |
                (round_r32_u32(color.b * denom) << 16) |
                ((e + 128) << 24));
    }

    return result;
}

internal void 
write_hdr_header(FILE *file, i32 width, i32 height)
{
    assert(width > 0 && height > 0);
    fprintf(file, "#?RADIANCE\n");
    fprintf(file, "FORMAT=32-bit_rle_rgbe\n\n");
    fprintf(file, "+Y %d +X %d\n", height, width);

#if 0
    if (info && (info->valid & RGBE_VALID_GAMMA)) {
        if (fprintf(fp,"GAMMA=%g\n",info->gamma) < 0)
            return rgbe_error(rgbe_write_error,NULL, errbuf);
    }

    if (info && (info->valid & RGBE_VALID_EXPOSURE)) {
        if (fprintf(fp,"EXPOSURE=%g\n",info->exposure) < 0)
            return rgbe_error(rgbe_write_error,NULL, errbuf);
    }

    if (fprintf(fp,"FORMAT=32-bit_rle_rgbe\n\n") < 0)
        return rgbe_error(rgbe_write_error,NULL, errbuf);
    if (fprintf(fp, "-Y %d +X %d\n", height, width) < 0)
        return rgbe_error(rgbe_write_error,NULL, errbuf);
#endif
}

int 
main(void)
{ 
    struct mach_timebase_info mach_time_info;
    mach_timebase_info(&mach_time_info);
    r32 nano_seconds_per_tick = ((r32)mach_time_info.numer/(r32)mach_time_info.denom);

    // NOTE(joon) Only using this to initialize the first series
    srand(time(NULL));
    RandomSeries series = start_random_series(rand()); 

    PlatformMemory platform_memory = {};

    platform_memory.permanent_memory_size = gigabytes(4);
    //platform_memory.transient_memory_size = gigabytes(1);
    u64 total_size = platform_memory.permanent_memory_size + platform_memory.transient_memory_size;
    vm_allocate(mach_task_self(), 
                (vm_address_t *)&platform_memory.permanent_memory,
                total_size, 
                VM_FLAGS_ANYWHERE);
    // TODO(joon) we dont need transient memory for now...
    //platform_memory.transient_memory = (u8 *)platform_memory.permanent_memory + platform_memory.permanent_memory_size;

    MemoryArena light_memory_arena = start_memory_arena(platform_memory.permanent_memory, megabytes(2));
    ParseSceneResult scene = {};
    // NOTE(joon) every memory allocation is explicit
    scene.light_push_buffer = start_temp_memory(&light_memory_arena, light_memory_arena.total_size);
    // TODO(joon) executable relative path?
    PlatformReadFileResult scene_file = debug_macos_read_file("/Volumes/digipen/offline_raytracer/data/testscene.scn");
    parse_scene(&scene, scene_file.memory, scene_file.size, "/Volumes/digipen/offline_raytracer/data/");
    scene.output_width = 1280;
    scene.output_height = 720;

    // NOTE(joon) manually create CSG
    MemoryArena csg_arena = start_memory_arena(malloc(megabytes(8)), megabytes(8));
    // TODO(joon) formalize csg population
    u32 csg_count = 0;
    CSG csgs[10] = {};
    csgs[0].sphere.center = v3_(0, 0, 0.8f);
    csgs[0].sphere.r = 0.35f;
    csgs[0].aab.min = v3_(0, 0, 0.8f) - v3_(0.3f, 0.3f, 0.3f);
    csgs[0].aab.max = v3_(0, 0, 0.8f) + v3_(0.3f, 0.3f, 0.3f);
    csgs[0].mat_index = 5;
    csg_count++;

    World world = {};
    world.ambient = scene.ambient;
    world.materials = scene.materials;
    world.mat_count = scene.mat_count;
    world.light_push_buffer = scene.light_push_buffer;
    world.light_count = scene.light_count;

    // TODO(joon) better make this better..
    Mesh meshes[100] = {};
    u32 mesh_count = 0;
    for(u32 mesh_index = 0;
            mesh_index < scene.mesh_count;
            mesh_index++)
    {
        MeshInfo *mesh_info = scene.mesh_infos + mesh_index;
        Mesh *mesh = meshes + mesh_count++;

        char extension_buffer[16] = {};
        get_extension(extension_buffer, mesh_info->file_path);
        if(string_compare(extension_buffer, "ply"))
        {
            PlatformReadFileResult mesh_file = debug_macos_read_file(mesh_info->file_path);
            ParsePlyHeaderResult ply_header = parse_ply_header(mesh_file.memory, mesh_file.size);

            mesh->vertex_count = ply_header.vertex_count;
            mesh->index_count = ply_header.index_count;
            mesh->vertices = (v3 *)malloc(sizeof(v3) * mesh->vertex_count);
            mesh->indices = (u32 *)malloc(sizeof(u32) * mesh->index_count);
            mesh->mat_index = mesh_info->mat_index;

            parse_ply(&ply_header, mesh_file.memory, mesh_file.size, mesh->vertices, mesh->indices);
        }
        else if(string_compare(extension_buffer, "obj"))
        {
            // replaced to obj, as I don't think we can benefit much from .x files
            PlatformReadFileResult mesh_file = debug_macos_read_file(mesh_info->file_path);

            PreParseObjResult pre_parse_obj_result = pre_parse_obj(mesh_file.memory, mesh_file.size);

            mesh->vertex_count = pre_parse_obj_result.position_count;
            mesh->index_count = pre_parse_obj_result.index_count;
            mesh->vertices = (v3 *)malloc(sizeof(v3) * mesh->vertex_count);
            mesh->indices = (u32 *)malloc(sizeof(u32) * mesh->index_count);
            mesh->mat_index = mesh_info->mat_index;

            parse_obj(&pre_parse_obj_result, mesh_file.memory, mesh_file.size, mesh->vertices, 0, 0, mesh->indices);
        }

        v3 min = v3_(Flt_Max, Flt_Max, Flt_Max);
        v3 max = v3_(Flt_Min, Flt_Min, Flt_Min);

        v4 axis_q = quaternion(v3_(0, 0, 1), 0.0174533f * mesh_info->axis.degree);
        v4 q = quaternion_mul(axis_q, mesh_info->quaternion); 

        // adjust vertices, based on rotation & translation & scale
        // also, build a very simple bounding box
        for(u32 vertex_index = 0;
                vertex_index < mesh->vertex_count;
                ++vertex_index)
        {
            v3 *v = mesh->vertices + vertex_index;

            // First rotate, scale, and then translate
            *v *= mesh_info->scale;  
            //*v = quaternion_rotation(mesh_info->quaternion, *v);
            *v = quaternion_rotation(mesh_info->quaternion, quaternion_rotation(v3_(0, 1, 0), 0.0174533f * mesh_info->axis.degree, *v));
            *v += mesh_info->translate;  

            // build bounding volume 
            min.x = minimum(min.x, v->x);
            min.y = minimum(min.y, v->y);
            min.z = minimum(min.z, v->z);

            max.x = maximum(max.x, v->x);
            max.y = maximum(max.y, v->y);
            max.z = maximum(max.z, v->z);

            mesh->aabb_min = min;
            mesh->aabb_max = max;
        }
    }
    
    IntersectionTestResult adfasdf = ray_intersect_with_aab(v3_(-1.5f, -1.5F, -1.5f), v3_(1.5f, 1.5f, 1.5f), v3_(0, 0, 0), v3_(1, 1, 1));

    u32 bvh_node_arena_size = megabytes(512);
    MemoryArena bvh_node_arena = start_memory_arena((u8 *)light_memory_arena.base + light_memory_arena.total_size, bvh_node_arena_size);

    BVHOctreeNode *top_most_node = push_struct(&bvh_node_arena, BVHOctreeNode);
    zero(top_most_node);
    top_most_node->aabb_min = v3_(Flt_Max, Flt_Max, Flt_Max);
    top_most_node->aabb_max = v3_(Flt_Min, Flt_Min, Flt_Min);

    for(u32 mesh_index = 0;
            mesh_index < scene.mesh_count;
            ++mesh_index)
    {
        Mesh *mesh = meshes + mesh_index;

        
        update_aabb_min_max(&top_most_node->aabb_min, &top_most_node->aabb_max, mesh, Shape_Type_Mesh);
    }

    for(u32 cylinder_index = 0;
            cylinder_index < scene.cylinder_count;
            ++cylinder_index)
    {
        Cylinder *cylinder = scene.cylinders + cylinder_index;
        update_aabb_min_max(&top_most_node->aabb_min, &top_most_node->aabb_max, cylinder, Shape_Type_Cylinder);
    }

    for(u32 box_index = 0;
            box_index < scene.box_count;
            ++box_index)
    {
        AAB *box = scene.boxes + box_index;
        update_aabb_min_max(&top_most_node->aabb_min, &top_most_node->aabb_max, box, Shape_Type_AAB);
    }

    for(u32 sphere_index = 0;
            sphere_index < scene.sphere_count;
            ++sphere_index)
    {
        Sphere *sphere = scene.spheres + sphere_index;
        update_aabb_min_max(&top_most_node->aabb_min, &top_most_node->aabb_max, sphere, Shape_Type_Sphere);
    }

    for(u32 csg_index = 0;
            csg_index < csg_count;
            ++csg_index)
    {
        CSG *csg = csgs + csg_index;
        csg->aabb_min = Flt_Max * v3_(1, 1, 1);
        csg->aabb_max = Flt_Min * v3_(1, 1, 1);
        update_aabb_min_max(&csg->aabb_min, &csg->aabb_max, &csg->sphere, Shape_Type_Sphere);
        update_aabb_min_max(&csg->aabb_min, &csg->aabb_max, &csg->aab, Shape_Type_AAB);
    }

    v3 top_most_node_center = 0.5f * (top_most_node->aabb_min + top_most_node->aabb_max);
    v3 top_most_node_half_dim = top_most_node->aabb_max - top_most_node_center;
    
    u32 desired_depth = 10;

    // NOTE(joon) get the top-most bounding volume node, which will contain every objects in the game
#if 1
    for(u32 mesh_index = 0;
            mesh_index < mesh_count;
            ++mesh_index)
    {
        Mesh *mesh = meshes + mesh_index;

        for(u32 index_index = 0;
                index_index < mesh->index_count;
                index_index += 3)
        {
            u32 i_0 = mesh->indices[index_index];
            u32 i_1 = mesh->indices[index_index + 1];
            u32 i_2 = mesh->indices[index_index + 2];

            assert(i_0 < 0x0000ffff);
            assert(i_1 < 0x0000ffff);
            assert(i_2 < 0x0000ffff);

            Triangle triangle = {};
            triangle.mesh = mesh;
            triangle.i_0 = i_0;
            triangle.i_1 = i_1;
            triangle.i_2 = i_2;

            push_shape_inside_node(&bvh_node_arena, top_most_node, top_most_node_center, top_most_node_half_dim, 0, desired_depth, &triangle, Shape_Type_Triangle);
        }
    }
#endif


    for(u32 cylinder_index = 0;
            cylinder_index < scene.cylinder_count;
            ++cylinder_index)
    {
        Cylinder *cylinder = scene.cylinders + cylinder_index;
        push_shape_inside_node(&bvh_node_arena, top_most_node, top_most_node_center, top_most_node_half_dim, 0, desired_depth, cylinder, Shape_Type_Cylinder);
    }

    for(u32 box_index = 0;
            box_index < scene.box_count;
            ++box_index)
    {
        AAB *aab = scene.boxes + box_index;
        push_shape_inside_node(&bvh_node_arena, top_most_node, top_most_node_center, top_most_node_half_dim, 0, desired_depth, aab, Shape_Type_AAB);
    }

    for(u32 sphere_index = 0;
            sphere_index < scene.sphere_count;
            ++sphere_index)
    {
        Sphere *sphere = scene.spheres + sphere_index;
        push_shape_inside_node(&bvh_node_arena, top_most_node, top_most_node_center, top_most_node_half_dim, 0, desired_depth, sphere, Shape_Type_Sphere);
    }

    for(u32 csg_index = 0;
            csg_index < csg_count;
            ++csg_index)
    {
        CSG *csg = csgs + csg_index;
        push_shape_inside_node(&bvh_node_arena, top_most_node, top_most_node_center, top_most_node_half_dim, 0, desired_depth, csg, Shape_Type_CSG);
    }

    u32 shape_arena_size = megabytes(8);
    MemoryArena shape_arena = start_memory_arena((u8 *)bvh_node_arena.base + bvh_node_arena.total_size, shape_arena_size);

    u32 total_shape_count = scene.sphere_count + scene.cylinder_count + + scene.box_count + mesh_count;
    ValidateNodesResult validate_nodes_result = {};
    validate_nodes_and_reallocate_shapes(&shape_arena, top_most_node, &validate_nodes_result);

    // TODO(joon) bunny!
    v3 *output_buffer = (v3 *)malloc(sizeof(v3) * scene.output_width * scene.output_height);

    Camera camera = {};
    camera.p = scene.camera_p;
    f32 rx = scene.camera_height_ratio * ((f32)scene.output_width / scene.output_height);

    camera.x_axis = rx * quaternion_rotation(scene.camera_quaternion, v3_(1, 0, 0));
    camera.y_axis = scene.camera_height_ratio * quaternion_rotation(scene.camera_quaternion, v3_(0, 1, 0));
    camera.z_axis = quaternion_rotation(scene.camera_quaternion, v3_(0, 0, 1));

    v3 z_axis = normalize(camera.p) - v3_(-0.3f, -0.5f, 2.8f); // - V3(0, 0, 0), which is the center of the world
    v3 x_axis = normalize(cross(v3_(0, 0, 1), camera.z_axis));
    v3 y_axis = normalize(cross(camera.z_axis, camera.x_axis));

    f32 ldfdfength = length(z_axis);

    // NOTE(joon) threads
    semaphore = dispatch_semaphore_create(0);

    pthread_attr_t thread_attribute = {};
    pthread_attr_init(&thread_attribute);

    ThreadWorkQueue work_queue = {};
    work_queue.item_count = 1024;
    work_queue.items = (ThreadWorkItem *)malloc(sizeof(ThreadWorkItem) * work_queue.item_count);

    u32 thread_spawn_count = 8; 
    //u32 thread_count = 1;
    // NOTE(joon): spawn threads
    if(thread_spawn_count > 0)
    {
        // TODO(joon) : spawn threads based on the core count
        MacOSThread *threads = (MacOSThread *)malloc(sizeof(MacOSThread) * thread_spawn_count);

        for(u32 thread_index = 0;
                thread_index < thread_spawn_count;
                ++thread_index)
        {
            MacOSThread *thread = threads + thread_index;
            thread->ID = thread_index + 1; // 0th thread is the main thread
            thread->queue = &work_queue;

            pthread_t thread_id; // we don't care about this one, we just generate our own id
            int result = pthread_create(&thread_id, &thread_attribute, &thread_proc, (void *)(thread));

            if(result != 0) // 0 means success
            {
                assert(0);
            }
        }
    }

    IntersectionTestResult a = ray_intersect_with_sphere(v3_(0, 0, 0), 10.0f, v3_(0, 0, 0), v3_(2, 2, 2));

    i32 tile_x_count = 32;
    i32 tile_y_count = 32;
    i32 pixel_count_per_tile_x = ceil_r32_i32(scene.output_width/(f32)tile_x_count);
    i32 pixel_count_per_tile_y = ceil_r32_i32(scene.output_height/(f32)tile_y_count);

    ThreadWorkRaytraceTileData *raytrace_tile_datas = (ThreadWorkRaytraceTileData *)malloc(sizeof(ThreadWorkRaytraceTileData) * tile_x_count * tile_y_count);
    u32 raytracer_data_index = 0;

    world.total_tile_to_render_count = tile_x_count * tile_y_count;

    u32 rays_per_pixel_count = 2048;
    for(i32 tile_y = 0;
            tile_y < tile_y_count;
            ++tile_y)
    {
        i32 min_y = tile_y * pixel_count_per_tile_y;
        i32 one_past_max_y = min_y + pixel_count_per_tile_y;

        if(one_past_max_y > scene.output_height)
        {
            one_past_max_y = scene.output_height;
        }

        for(i32 tile_x = 0;
                tile_x < tile_x_count;
                ++tile_x)
        {
            i32 min_x = tile_x * pixel_count_per_tile_x;
            i32 one_past_max_x = min_x + pixel_count_per_tile_x;

            if(one_past_max_x > scene.output_width)
            {
                one_past_max_x = scene.output_width;
            }

            ThreadWorkRaytraceTileData *data = raytrace_tile_datas + raytracer_data_index++;
            data->world = &world;
            data->camera = &camera;
            data->top_most_node = top_most_node;
            data->output_buffer = output_buffer;
            data->output_width = scene.output_width;
            data->output_height = scene.output_height;
            data->random_series = start_random_series(random_u32(&series));

            data->tile_x = tile_x;
            data->tile_y = tile_y;

            data->tile_min_x = min_x;
            data->tile_one_past_max_x = one_past_max_x;

            data->tile_min_y = min_y;
            data->tile_one_past_max_y = one_past_max_y;

            data->rays_per_pixel_count = rays_per_pixel_count;
            data->russian_roulette_value = 0.8f;

            macos_add_thread_work_item(&work_queue,
                                        thread_callback_tiled_raytrace,
                                        (void *)data);
        }
    }
    u64 start_time = mach_absolute_time();

    macos_complete_all_thread_work_queue_items(&work_queue);

    // NOTE(joon) This step is required as completing all thread work only gurantees that all the works are fired up, 
    // but does not gurantee that all of them are done
    while(world.total_tile_to_render_count != world.rendered_tile_count)
    {
    }

    u64 nano_seconds_elapsed = mach_time_diff_in_nano_seconds(start_time, mach_absolute_time(), nano_seconds_per_tick);
    u32 ms_elapsed = (u32)(nano_seconds_elapsed*nano_sec_to_milli_sec);
    u32 sec_elapsed = (u32)(ms_elapsed / sec_to_milli_sec);

    u32 total_bytes_used = light_memory_arena.used + bvh_node_arena.used + shape_arena.used;

    // TODO(joon) Also add time per bounce
    printf("Raytracing finished in %d seconds(%dms) : %dThreads, %dRays Per Pixel, %dMB used for the BVH", sec_elapsed, ms_elapsed, thread_spawn_count + 1, rays_per_pixel_count, (total_bytes_used / 1024) / 1024);

    // NOTE(joon) output to hdr file
    FILE *hdr_file = fopen("../output.hdr", "wb");
    write_hdr_header(hdr_file, scene.output_width, scene.output_height);

    v3 *row = output_buffer + (scene.output_height-1)*scene.output_width;
    for(i32 y = 0;
            y < scene.output_height;
            ++y)
    {
        v3 *pixel = row;
        for(i32 x = 0;
                x < scene.output_width;
                ++x)
        {

            u32 color_u32 = v3_to_rgbe(*pixel);
            fwrite(&color_u32, sizeof(u32), 1, hdr_file);

            pixel++;
        }

        // bottom up
        row -= scene.output_width;
    }

    fclose(hdr_file);

    return 0;
}










