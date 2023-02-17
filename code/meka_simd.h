#ifndef MEKA_SIMD_H
#define MEKA_SIMD_H
// NOTE(joon): This usage of SIMD is a very simple & easy to use to make use of SIMD in a particular peace of code that you want to speed up.
// However, this also creates a bit of an overhead, as everything is wrapped into another layer of function call that might not get inlined or whatever weird thing that
// the compiler might do... So for a routine that needs an _extreme_ speed boot, use this as just a intermission layer, and eventually get rid of it 
// or maybe write an assmebly code

// TODO(joon): Make seperate types for each lane with different size, or a safegurad that we can insert inside the function that tells 
// that a certain function only works with specific size of lane
#define LANE_WIDTH 4 // NOTE(joon): Should be 1/4 for ARM, 1/4/8 for x86/64

#if MEKA_ARM
//NOTE(joon): Gets rid of bl instructions
#define force_inline static inline __attribute__((always_inline))

#if LANE_WIDTH == 4
#include "meka_simd_4x.h"

#elif LANE_WIDTH == 8

// TODO(joon): M1 pro does not support more that 128bit lane... 
//#include "meka_simd_8x.h"
#elif LANE_WIDTH == 1

typedef v3 LaneV3;
typedef r32 LaneF32;
typedef u32 LaneU32;
typedef i32 LaneI32;
#define LaneRandomSeries RandomSeries

//////////////////// LaneU32 //////////////////// 
force_inline LaneU32
LaneU32_(u32 value)
{
    return value;
}

force_inline LaneU32
LaneU32_load(u32 *ptr)
{
    LaneU32 result = *ptr;
    return result;
}

force_inline u32
get_lane(LaneU32 a, u32 lane)
{
    return a;
}

force_inline LaneF32
convert_f32_from_u32(LaneU32 a)
{
    // TODO(joon): This might have potential sign bug!!!
    return (LaneF32)a;
}

force_inline LaneU32
compare_not_equal(LaneU32 a, LaneU32 b)
{
    LaneU32 result = 0;
    if(a != b)
    {
        result = 0xffffffff;
    }

    return result;
}

force_inline LaneU32
overwrite(LaneU32 dest, LaneU32 mask, LaneU32 source)
{
    LaneU32 result = {};

    mask = (mask == 0) ? 0 : 0xffffffff;
    result = (dest & (~mask)) | (source & mask);

    return result;
}

force_inline LaneU32
is_lane_non_zero(LaneU32 a)
{
    LaneU32 result = 0;

    if(a != 0)
    {
        result = 0xffffffff;
    }

    return result;
}

force_inline u32
get_non_zero_lane_count_from_all_set_bit(LaneU32 a)
{
    u32 result = 0;

    if(a != 0)
    {
        result = 1;
    }

    return result;
}

force_inline u32
get_non_zero_lane_count(LaneU32 a)
{
    u32 result = 0;

    if(a != 0)
    {
        result = 1;
    }

    return result;
}

force_inline b32
all_lanes_zero(LaneU32 a)
{
    b32 result = false;

    if(a == 0)
    {
        result = true;
    }

    return result;
}

//////////////////// LaneF32 //////////////////// 

force_inline LaneF32
LaneF32_(f32 value)
{
    return value;
}

force_inline LaneF32
LaneF32_(f32 value0, f32 value1, f32 value2, f32 value3)
{
    return value0;
}

force_inline LaneF32
LaneF32_load(f32 *ptr)
{
    LaneF32 result = *ptr;
    return result;
}

force_inline f32
get_lane(LaneF32 a, u32 lane)
{
    return a;
}

force_inline LaneU32
compare_equal(LaneF32 a, LaneF32 b)
{
    LaneU32 result = 0;

    // TODO(joon): This is a totally made up number
    f32 epsilon = 0.00001f;
    f32 diff = a-b;

    if(diff > -epsilon && diff < epsilon)
    {
        result = 0xffffffff;
    }

    return result;
}

force_inline LaneU32
compare_not_equal(LaneF32 a, LaneF32 b)
{
    LaneU32 result = 0;

    result = (compare_equal(a, b) == 0) ? 0xffffffff : 0;

    return result;
}

force_inline LaneU32
compare_greater(LaneF32 a, LaneF32 b)
{
    LaneU32 result = 0;

    if(a > b)
    {
        result = 0xffffffff;
    }

    return result;
}

// TODO(joon): Also add epsilon?
force_inline LaneU32
compare_greater_equal(LaneF32 a, LaneF32 b)
{
    LaneU32 result = 0;

    if(a > b || compare_equal(a, b))
    {
        result = 0xffffffff;
    }

    return result;
}

force_inline LaneU32
compare_less(LaneF32 a, LaneF32 b)
{
    LaneU32 result = 0;
    if(a < b)
    {
        result = 0xffffffff;
    }

    return result;
}

force_inline LaneU32
compare_less_equal(LaneF32 a, LaneF32 b)
{
    LaneU32 result = 0;
    if(a <= b)
    {
        result = 0xffffffff;
    }

    return result;
}

force_inline LaneF32
overwrite(LaneF32 dest, LaneU32 mask, LaneF32 source)
{
    LaneF32 result = {};

    mask = (mask == 0) ? 0 : 0xffffffff;

    // TODO(joon): This is also prone to bug
    LaneU32 result_u32 = (((*(LaneU32 *)&dest) & (~mask)) | ((*(LaneU32 *)&source) & mask));
    result = *(LaneF32 *)&result_u32;

    return result;
}

force_inline LaneF32
add_all_lanes(LaneF32 a)
{
    return a;
}

//////////////////// LaneV3 //////////////////// 

force_inline LaneV3
LaneV3_()
{
    LaneV3 result = {};
    return result;
}

force_inline LaneV3
LaneV3_(v3 value)
{
    return value;
}

force_inline LaneV3
LaneV3_(v3 value0, v3 value1, v3 value2, v3 value3)
{
    LaneV3 result = value0;

    return result;
}

force_inline LaneV3
LaneV3_(LaneF32 value0, LaneF32 value1, LaneF32 value2)
{
    LaneV3 result = {};

    result.x = value0;
    result.y = value1;
    result.z = value2;

    return result;
}

force_inline LaneV3
operator*(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};

    result.x = a.x * b.x;
    result.y = a.y * b.y;
    result.z = a.z * b.z;

    return result;
}

force_inline LaneV3
operator/(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};

    result.x = a.x / b.x;
    result.y = a.y / b.y;
    result.z = a.z / b.z;

    return result;
}

force_inline LaneV3
min(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};
    result.x = minimum(a.x, b.x);
    result.y = minimum(a.y, b.y);
    result.z = minimum(a.z, b.z);

    return result;
}

force_inline LaneV3
max(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};

    result.x = maximum(a.x, b.x);
    result.y = maximum(a.y, b.y);
    result.z = maximum(a.z, b.z);

    return result;
}

force_inline LaneF32
min_component(LaneV3 a)
{
    LaneF32 result = minimum(minimum(a.x, a.y), a.z);

    return result;
}

force_inline LaneF32
max_component(LaneV3 a)
{
    LaneF32 result = maximum(maximum(a.x, a.y), a.z);

    return result;
}

force_inline LaneV3
overwrite(LaneV3 dest, LaneU32 mask, LaneV3 source)
{
    LaneV3 result = {};

    mask = (mask == 0) ? 0 : 0xffffffff;

    result.x = overwrite(dest.x, mask, source.x);
    result.y = overwrite(dest.y, mask, source.y);
    result.z = overwrite(dest.z, mask, source.z);

    return result;
}

force_inline LaneRandomSeries
start_random_series(u32 seed0, u32 seed1, u32 seed2, u32 seed3)
{
    LaneRandomSeries result = {};
    result.next_random = LaneU32_(seed0);

    return result;
}

#endif // #elif LANE_WIDTH == 1

///////////// codes that are common across different simd lane size should come here

#endif
