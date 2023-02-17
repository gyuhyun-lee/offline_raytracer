#ifndef RAY_H
#define RAY_H

struct Sphere
{
    v3 center;
    f32 r;

    u32 mat_index;
};

struct AAB
{
    v3 min;
    v3 max;

    u32 mat_index;
};

struct Cylinder
{
    v3 base;
    v3 axis;
    f32 r;

    u32 mat_index;
};

// NOTE(joon) For any material, both diffuse and transmission should not be 0
struct Material
{
    v3 diffuse;
    v4 specular;

    v3 transmission;
    f32 ior;

    v3 emit_color;
    b32 is_light;
};

struct Camera
{
    v3 p;
    
    v3 x_axis;
    v3 y_axis;
    v3 z_axis;
};

struct Mesh
{
    // NOTE(joon) These are already rotated / scaled / translated
    v3 *vertices;
    u32 vertex_count;

    u32 *indices;
    u32 index_count;

    u32 mat_index;

    // aabb is defined at the start of the program
    v3 aabb_min;
    v3 aabb_max;
};

struct Triangle
{
    u32 i_0;
    u32 i_1;
    u32 i_2;

    Mesh *mesh;
};

struct World
{
    v3 ambient;

    u32 total_tile_to_render_count;
    volatile u32 rendered_tile_count;

    Material *materials;
    u32 mat_count;

    TempMemory light_push_buffer;
    u32 light_count;
};

#define BVH_Node_Pos_X_Mask 0b10101010
#define BVH_Node_Neg_X_Mask 0b01010101
#define BVH_Node_Pos_Y_Mask 0b11001100
#define BVH_Node_Neg_Y_Mask 0b00110011
#define BVH_Node_Pos_Z_Mask 0b11110000
#define BVH_Node_Neg_Z_Mask 0b00001111

enum  ShapeType
{
    Shape_Type_Null,
    Shape_Type_Sphere,
    Shape_Type_Cylinder,
    Shape_Type_AAB, // axis aligned box
    Shape_Type_Mesh, // axis aligned box
    Shape_Type_Triangle,
    Shape_Type_CSG,
};

struct BVHShapeHeader
{
    ShapeType type;
};

// NOTE(joon) depth 0 : only the top most node exists
// depth 1 : top most & 8 child nodes exist
struct BVHOctreeNode
{
    u8 child_bitmask;
    BVHOctreeNode *first_child; // 8 children nodes will be laid sequentially in memory

    // NOTE(joon) This normally contains only one shape,
    // but the very last leaf nodes might contain more than one shapes
    // TODO(joon) should we only save the pointer to the shape, or the shape itself?
    // shape itself makes a little bit more sense, as then we can divide the mesh into triangles and 
    // save the triangles?

    b32 is_leaf;
    TempMemory push_buffer;

    // NOTE(joon) aabb that contains all the shapes that this node has
    // This _should_ contain the volumes
    v3 aabb_min;
    v3 aabb_max;
};

struct BVHQueue
{
    u8 *base;

    u32 used;
    u32 size;

    u32 read_cursor;
};

enum CSGShapeIncludeType
{
    CSGShapeIncludeType_Null,
    CSGShapeIncludeType_Intersect,
    CSGShapeIncludeType_Union,
    CSGShapeIncludeType_Difference,
};

struct CSGShapeHeader
{
    ShapeType type;
    CSGShapeIncludeType include_type;
};

// NOTE(joon) structure of this is header/shape/header/shape...
struct CSG
{
#if 0
    // NOTE(joon) Shape push buffer for any shape 
    u8 *base;
    u32 size;
    u32 used;
#endif

    Sphere sphere;
    AAB aab;

    v3 aabb_min;
    v3 aabb_max;

    u32 mat_index;
};

#endif
