#ifndef MEKA_SIMD_4X_H
#define MEKA_SIMD_4X_H

#if MEKA_ARM
#include <arm_neon.h>

// TODO(note): Doing simd operations this way makes it easier for us to see and examine the code,
// but has a chance to make the performance worse as the compiler is not that smart :(
struct LaneF32
{
    float32x4_t v;
};

struct LaneU32
{
    uint32x4_t v;
};

struct LaneI32
{
    int32x4_t v;
};

struct LaneV3
{
    union
    {
        struct
        {
            float32x4_t x;
            float32x4_t y;
            float32x4_t z;
        };

        struct
        {
            float32x4_t r;
            float32x4_t g;
            float32x4_t b;
        };
        float32x4_t e[3];
    };
};

struct LaneV3i
{
    int32x4_t x;
    int32x4_t y;
    int32x4_t z;
};

struct LaneV3u
{
    uint32x4_t x;
    uint32x4_t y;
    uint32x4_t z;
};

force_inline u32
add_all_lanes(uint32x4_t a)
{
    u32 result = vaddvq_u32(a);

    return result;
}

force_inline f32
add_all_lanes(float32x4_t a)
{
    f32 result = vaddvq_f32(a);

    return result;
}

//////////////////// LaneU32 ////////////////////

force_inline LaneU32
LaneU32_(u32 dup)
{
    LaneU32 result = {};

    result.v = vdupq_n_u32(dup);

    return result;
}

force_inline LaneU32
LaneU32_(u32 value0, u32 value1, u32 value2, u32 value3)
{
    LaneU32 result = {};
    
    // TODO(joon): is using vld1q_u32 faster?
    result.v = {value0, value1, value2, value3};

    return result;
}

force_inline LaneU32
LaneU32_load(u32 *ptr)
{
    LaneU32 result = {};
    result.v = vld1q_u32(ptr);

    return result;
}

force_inline u32
get_lane(LaneU32 a, u32 lane)
{
    return (((u32 *)&a)[lane]);
}

// unary not opeartor
force_inline LaneU32
operator~(LaneU32 a)
{
    LaneU32 result = {};

    result.v = vmvnq_u32(a.v);

    // this is not a appropriate place to be here..
    return result;
}

force_inline LaneU32
operator+(LaneU32 a, LaneU32 b)
{
    LaneU32 result = {};
    result.v = vaddq_u32(a.v, b.v);

    return result;
}

force_inline LaneU32
operator-(LaneU32 a, LaneU32 b)
{
    LaneU32 result = {};
    result.v = vsubq_u32(a.v, b.v);

    return result;
}

force_inline LaneU32
operator*(LaneU32 a, LaneU32 b)
{
    LaneU32 result = {};
    result.v = vmulq_u32(a.v, b.v);

    return result;
}

/*
   NOTE(joon): Types of shifting
    1. Rounding : When underflowed, the bits are rounded instead of being discarded
    ex) 1111 >> 1 = 1000(which is 111 + 1)
    2. Saturating : When overflowed, all bits are set to 1 
    ex) 0xff000000 << 1 = 0xffffffff

    3. Normal : When left shifted, LSB are set to 0
                When right shifted, MSB are set to 1
*/

// TODO(joon): There's no simple instruction in ARM for left & right shift just like in intel,
// as the simple shift requires compile time constant shift values.
// For now, I am using normal shift for shifting left, and saturating shift for shifting right
// Which is error prone when if the sign of the input shift values are negative
force_inline LaneU32
operator<<(LaneU32 a, u32 shift_amount)
{
    LaneU32 result = {};

    int32x4_t shift = vdupq_n_u32(shift_amount);

    // NOTE(joon): using normal shift for shifting left
    result.v = vshlq_s32(a.v, shift);

    return result;
}

// NOTE(joon): logical shift, does not keep the sign bit
force_inline LaneU32
operator>>(LaneU32 a, u32 shift_amount)
{
    LaneU32 result = {};

    int32x4_t shift = vdupq_n_u32(-shift_amount);
    // NOTE(joon): using saturating shift for shifting right, as saturation does not have any effect when underflowing
    result.v = vqshlq_u32(a.v, shift);

    return result;
}

force_inline LaneU32
operator|(LaneU32 a, LaneU32 b)
{
    LaneU32 result = {};
    result.v = vorrq_u32(a.v, b.v);

    return result;
}

force_inline LaneU32
operator&(LaneU32 a, LaneU32 b)
{
    LaneU32 result = {};
    result.v = vandq_u32(a.v, b.v);

    return result;
}



force_inline LaneU32
compare_equal(LaneU32 a, LaneU32 b)
{
    LaneU32 result = {};

    result.v = vceqq_u32(a.v, b.v);

    return result;
}

force_inline LaneU32
compare_greater_equal(LaneU32 a, LaneU32 b)
{
    LaneU32 result = {};

    result.v = vcgeq_u32(a.v, b.v);

    return result;
}

force_inline LaneU32
compare_greater(LaneU32 a, LaneU32 b)
{
    LaneU32 result = {};

    result.v = vcgtq_u32(a.v, b.v);

    return result;
}

force_inline LaneU32
compare_less_equal(LaneU32 a, LaneU32 b)
{
    LaneU32 result = {};

    result.v = vcleq_u32(a.v, b.v);

    return result;
}

force_inline LaneU32
compare_less(LaneU32 a, LaneU32 b)
{
    LaneU32 result = {};

    result.v = vcltq_u32(a.v, b.v);

    return result;
}

// NOTE(joon): clear the values that will be overwritten by the new value, or them to overwrite
force_inline LaneU32
overwrite(LaneU32 dest, LaneU32 mask, LaneU32 source)
{
    LaneU32 result = {};

    // TODO(joonl: might be better to just write the reinterpret function...
    result = (dest & (~mask)) | (source & mask);

    return result;
}

// NOTE(joon): Fills the mask lane with 1 if the lane has non-zero value
// if not, fills with 0
force_inline LaneU32
is_lane_non_zero(LaneU32 a)
{
    LaneU32 result = {};
    
    LaneU32 LaneU32_0 = LaneU32_(0);
    result = ~compare_equal(a, LaneU32_0);

    return result;
}

force_inline u32
add_all_lanes(LaneU32 a)
{
    u32 result = vaddvq_u32(a.v);

    return result;
}

// TODO(joon): Should be a magic intrinsic function
// that does this for me
// NOTE(joon): This only works if the lane is all 1 or 0
force_inline u32
get_non_zero_lane_count_from_all_set_bit(LaneU32 a)
{
    u32 result = 0;

    LaneU32 LaneU32_1 = LaneU32_(1);

    result = add_all_lanes(a & LaneU32_1);
    //result = get_lane(a, 0) || get_lane(a, 1) || get_lane(a, 2) || get_lane(a, 3);

    return result;
}

// TODO(joon): Should be a magic intrinsic function
// that does this for me
// TODO(joon): This routine is not tested yet
force_inline u32
get_non_zero_lane_count(LaneU32 a)
{
    u32 result = 0;

    LaneU32 LaneU32_0 = LaneU32_(0);
    result = add_all_lanes(~compare_equal(a, LaneU32_0))/8;

    return result;
}

force_inline b32
all_lanes_zero(LaneU32 a)
{
    b32 result = !(add_all_lanes(a));

    return result;
}

//////////////////// LaneF32 //////////////////// 

force_inline LaneF32
LaneF32_(r32 dup)
{
    LaneF32 result = {};
    result.v = vdupq_n_f32(dup);

    return result;
}

force_inline LaneF32
LaneF32_(r32 value0, r32 value1, r32 value2, r32 value3)
{
    // TODO(joon): I assume that instead of using this path, it's much better to structure some kind of SOA
    LaneF32 result = {};
    result.v = {value0, value1, value2, value3};

    return result;
}

force_inline LaneF32
LaneF32_load(r32 *ptr)
{
    LaneF32 result = {};
    result.v = vld1q_f32(ptr);

    return result;
}

force_inline r32
get_lane(LaneF32 a, u32 lane)
{
    return (((r32 *)&a)[lane]);
}

force_inline LaneF32
operator+(LaneF32 a, LaneF32 b)
{
    LaneF32 result = {};
    result.v = vaddq_f32(a.v, b.v);

    return result;
}

force_inline LaneF32
operator-(LaneF32 a, LaneF32 b)
{
    LaneF32 result = {};
    result.v = vsubq_f32(a.v, b.v);

    return result;
}

force_inline LaneF32
operator-(LaneF32 a)
{
    // TODO(joon): this means we create this constant value every time...
    // any way to do this more efficiently?
    float32x4_t minus_1 = vdupq_n_f32(-1.0f);

    a.v = vmulq_f32(a.v, minus_1);

    return a;
}

force_inline LaneF32
operator*(LaneF32 a, LaneF32 b)
{
    LaneF32 result = {};
    result.v = vmulq_f32(a.v, b.v);

    return result;
}
force_inline LaneF32
operator/(LaneF32 a, LaneF32 b)
{
    LaneF32 result = {};

    result.v = vdivq_f32(a.v, b.v);

    return result;
}

force_inline LaneF32 &
operator+=(LaneF32 &a, LaneF32 b)
{
    a.v = vaddq_f32(a.v, b.v);

    return a;
}

force_inline LaneF32 &
operator-=(LaneF32 &a, LaneF32 b)
{
    a.v = vsubq_f32(a.v, b.v);

    return a;
}

force_inline LaneF32 &
operator/=(LaneF32 &a, LaneF32 b)
{
    a.v = vdivq_f32(a.v, b.v);

    return a;
}

force_inline LaneU32
convert_u32_from_f32(LaneF32 value)
{
    LaneU32 result = {};
    result.v = vcvtq_u32_f32(value.v);

    return result;
}

force_inline LaneF32
convert_f32_from_u32(LaneU32 value)
{
    LaneF32 result = {};
    result.v = vcvtq_f32_u32(value.v);

    return result;
}

force_inline LaneF32
operator|(LaneF32 a, LaneU32 mask)
{
    LaneF32 result = {};

    result.v = vreinterpretq_f32_u32(vorrq_u32(vreinterpretq_u32_f32(a.v), mask.v));

    return result;
}

// NOTE(joon):using this function might add up ambiguity
#if 1
force_inline LaneF32
operator&(LaneF32 a, LaneU32 mask)
{
    LaneF32 result = {};

    result.v = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(a.v), mask.v));

    return result;
}
#endif

force_inline LaneF32
operator|(LaneF32 a, LaneF32 b)
{
    LaneF32 result = {};
    
    result.v = vreinterpretq_f32_u32(vreinterpretq_u32_f32(a.v) | vreinterpretq_u32_f32(b.v));

    return result;
}

force_inline LaneF32
sqrt(LaneF32 a)
{
    LaneF32 result = {};

    result.v = vsqrtq_f32(a.v);

    return result;
}

// NOTE(joon): clear the values that will be overwritten by the new value, or them to overwrite
force_inline LaneF32
overwrite(LaneF32 dest, LaneU32 mask, LaneF32 source)
{
    LaneF32 result = {};

    result.v = vreinterpretq_f32_u32(((vreinterpretq_u32_f32(dest.v) & ~mask.v) | (vreinterpretq_u32_f32(source.v) & mask.v)));

    return result;
}

// TODO(joon): Not sure how accurate this floating point comparison is
force_inline LaneU32
compare_equal(LaneF32 a, LaneF32 b)
{
    LaneU32 result = {};

    result.v = vceqq_f32(a.v, b.v);

    return result;
}

force_inline LaneU32
compare_greater_equal(LaneF32 a, LaneF32 b)
{
    LaneU32 result = {};

    result.v = vcgeq_f32(a.v, b.v);

    return result;
}

force_inline LaneU32
compare_greater(LaneF32 a, LaneF32 b)
{
    LaneU32 result = {};

    result.v = vcgtq_f32(a.v, b.v);

    return result;
}

force_inline LaneU32
compare_less_equal(LaneF32 a, LaneF32 b)
{
    LaneU32 result = {};

    result.v = vcleq_f32(a.v, b.v);

    return result;
}

force_inline LaneU32
compare_less(LaneF32 a, LaneF32 b)
{
    LaneU32 result = {};

    result.v = vcltq_f32(a.v, b.v);

    return result;
}

force_inline LaneU32
compare_not_equal(LaneF32 a, LaneF32 b)
{
    LaneU32 result = {};

    result = ~compare_equal(a, b);

    return result;
}

force_inline LaneF32
lerp(LaneF32 min, LaneF32 t, LaneF32 max)
{
    return min + t*(max-min);
}

force_inline LaneF32
min(LaneF32 a, LaneF32 b)
{
    LaneF32 result = {};

    result.v = vminq_f32(a.v, b.v);

    return result;
}

force_inline f32
min_component(LaneF32 a)
{
    f32 result = vminvq_f32(a.v);

    return result;
}

force_inline f32
max_component(LaneF32 a)
{
    f32 result = vmaxvq_f32(a.v);

    return result;
}

force_inline LaneF32
max(LaneF32 a, LaneF32 b)
{
    LaneF32 result = {};

    result.v = vmaxq_f32(a.v, b.v);

    return result;
}



#if 0
// NOTE(joon): Fills the mask lane with 1 if the lane has non-zero value
// if not, fills with 0
force_inline LaneU32
is_lane_non_zero(LaneF32 a)
{
    LaneU32 result = {};
    
    LaneU32 LaneU32_0 = LaneU32_(0);

    result = ~compare_equal(vreinterpretq_u32_f32(a.v), LaneU32_0);

    return result;
}
#endif

force_inline f32
add_all_lanes(LaneF32 a)
{
    f32 result = vaddvq_f32(a.v);

    return result;
}

force_inline b32
all_lanes_zero(LaneF32 a)
{
    b32 result = !(add_all_lanes(a));

    return result;
}

//////////////////// LaneV3 ////////////////////

// TODO(joon): we can also make this from operator=
force_inline LaneV3
LaneV3_(v3 dup)
{
    LaneV3 result = {};

#if 1
    result.x = vdupq_n_f32(dup.x);
    result.y = vdupq_n_f32(dup.y);
    result.z = vdupq_n_f32(dup.z);
#else

    result.x = {dup.x, dup.x, dup.x, dup.x};
    result.y = {dup.y, dup.y, dup.y, dup.y};
    result.z = {dup.z, dup.z, dup.z, dup.z};
#endif


    return result;
}

// NOTE(joon): This path is most likely for the initial steps of opitimization, consider using the SOA version instead!
force_inline LaneV3
LaneV3_(v3 value0, v3 value1, v3 value2, v3 value3)
{
    LaneV3 result = {};

    result.x = {value0.x, value1.x, value2.x, value3.x};
    result.y = {value0.y, value1.y, value2.y, value3.y};
    result.z = {value0.z, value1.z, value2.z, value3.z};

    return result;
}

force_inline LaneV3
LaneV3_(LaneF32 value0, LaneF32 value1, LaneF32 value2)
{
    LaneV3 result = {};

    result.x = value0.v;    
    result.y = value1.v;
    result.z = value2.v;

    return result;
}

// NOTE(joon): This is a much more preferred version of loading in the data
force_inline LaneV3
LaneV3_load(r32 *array_of_x, r32 *array_of_y, r32 *array_of_z)
{
    LaneV3 result = {};

    result.x = vld1q_f32(array_of_x);
    result.y = vld1q_f32(array_of_y);
    result.z = vld1q_f32(array_of_z);

    return result;
}

force_inline LaneV3
operator+(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};

    result.x = vaddq_f32(a.x, b.x);
    result.y = vaddq_f32(a.y, b.y);
    result.z = vaddq_f32(a.z, b.z);

    return result;
}

force_inline LaneV3
operator-(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};
    result.x = vsubq_f32(a.x, b.x);
    result.y = vsubq_f32(a.y, b.y);
    result.z = vsubq_f32(a.z, b.z);

    return result;
}

force_inline LaneV3
operator*(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};

    result.x = vmulq_f32(a.x, b.x);
    result.y = vmulq_f32(a.y, b.y);
    result.z = vmulq_f32(a.z, b.z);

    return result;
}

force_inline LaneV3
operator*(LaneF32 value, LaneV3 b)
{
    LaneV3 result = {};

    result.x = vmulq_f32(value.v, b.x);
    result.y = vmulq_f32(value.v, b.y);
    result.z = vmulq_f32(value.v, b.z);

    return result;
}

force_inline LaneV3 &
operator*=(LaneV3 &v, LaneF32 value)
{
    v.x = vmulq_f32(v.x, value.v);
    v.y = vmulq_f32(v.y, value.v);
    v.z = vmulq_f32(v.z, value.v);

    return v;
}

force_inline LaneV3 &
operator/=(LaneV3 &v, LaneF32 value)
{
    v.x = vdivq_f32(v.x, value.v);
    v.y = vdivq_f32(v.y, value.v);
    v.z = vdivq_f32(v.z, value.v);

    return v;
}

force_inline LaneV3 &
operator+=(LaneV3 &v, LaneF32 value)
{
    v.x = vaddq_f32(v.x, value.v);
    v.y = vaddq_f32(v.y, value.v);
    v.z = vaddq_f32(v.z, value.v);

    return v;
}


force_inline LaneV3
operator/(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};

    result.x = vdivq_f32(a.x, b.x);
    result.y = vdivq_f32(a.y, b.y);
    result.z = vdivq_f32(a.z, b.z);

    return result;
}

force_inline LaneV3
operator/(LaneV3 a, LaneF32 b)
{
    LaneV3 result = {};

    result.x = vdivq_f32(a.x, b.v);
    result.y = vdivq_f32(a.y, b.v);
    result.z = vdivq_f32(a.z, b.v);

    return result;
}

force_inline LaneV3 
operator-(LaneV3 a)
{
    // TODO(joon): this means we create this constant value every time...
    // any way to do this more efficiently?
    f32 minus_1_[4] = {-1.0f, -1.0f, -1.0f, -1.0f};
    float32x4_t minus_1 = vld1q_f32(minus_1_);

    a.x = vmulq_f32(a.x, minus_1);
    a.y = vmulq_f32(a.y, minus_1);
    a.z = vmulq_f32(a.z, minus_1);

    return a;
}

force_inline LaneV3
overwrite(LaneV3 dest, LaneU32 mask, LaneV3 source)
{
    LaneV3 result = {};

    result.x = vreinterpretq_f32_u32(((vreinterpretq_u32_f32(dest.x) & ~mask.v) | (vreinterpretq_u32_f32(source.x) & mask.v)));
    result.y = vreinterpretq_f32_u32(((vreinterpretq_u32_f32(dest.y) & ~mask.v) | (vreinterpretq_u32_f32(source.y) & mask.v)));
    result.z = vreinterpretq_f32_u32(((vreinterpretq_u32_f32(dest.z) & ~mask.v) | (vreinterpretq_u32_f32(source.z) & mask.v)));

    return result;
}

force_inline LaneV3
lerp(LaneV3 min, LaneF32 t, LaneV3 max)
{
    LaneV3 result = {};

    result.x = vaddq_f32(min.x, vmulq_f32(t.v, vsubq_f32(max.x, min.x)));
    result.y = vaddq_f32(min.y, vmulq_f32(t.v, vsubq_f32(max.y, min.y)));
    result.z = vaddq_f32(min.z, vmulq_f32(t.v, vsubq_f32(max.z, min.z)));

    return result;
}

force_inline LaneF32
dot(LaneV3 a, LaneV3 b)
{
    LaneV3 hadamard = a*b;

    LaneF32 result = {};

    result.v = hadamard.x + hadamard.y + hadamard.z;

    return result;
}

force_inline LaneF32
length_square(LaneV3 a)
{
    return dot(a, a);
}

force_inline LaneF32
length(LaneV3 a)
{
    LaneF32 result = sqrt(length_square(a));

    return result;
}

force_inline LaneV3
normalize(LaneV3 a)
{
    LaneV3 result = a / length(a);

    return result;
}

force_inline LaneV3
cross(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};

    result.x = a.y*b.z - b.y*a.z;
    result.y = b.x*a.z - a.x*b.z;
    result.z = a.x*b.y - b.x*a.y;

    return result;
}

force_inline LaneV3
min(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};

    result.x = vminq_f32(a.x, b.x);
    result.y = vminq_f32(a.y, b.y);
    result.z = vminq_f32(a.z, b.z);

    return result;
}

force_inline LaneV3
max(LaneV3 a, LaneV3 b)
{
    LaneV3 result = {};

    result.x = vmaxq_f32(a.x, b.x);
    result.y = vmaxq_f32(a.y, b.y);
    result.z = vmaxq_f32(a.z, b.z);

    return result;
}

force_inline LaneF32
min_component(LaneV3 a)
{
    LaneF32 result = {};
    result.v = vminq_f32(vminq_f32(a.x, a.y), a.z);

    return result;
}

force_inline LaneF32
max_component(LaneV3 a)
{
    LaneF32 result = {};
    result.v = vmaxq_f32(vmaxq_f32(a.x, a.y), a.z);

    return result;
}


force_inline v3
add_all_lanes(LaneV3 a)
{
    v3 result = {};

    result.x = vaddvq_f32(a.x);
    result.y = vaddvq_f32(a.y);
    result.z = vaddvq_f32(a.z);

    return result;
}

// TODO(joon): Not sure how accurate this floating point comparison is
force_inline LaneU32
compare_equal(LaneV3 a, LaneV3 b)
{
    LaneU32 result = {};

    result.v = vceqq_f32(a.x, b.x) & vceqq_f32(a.y, b.y) & vceqq_f32(a.z, b.z);

    return result;
}

force_inline LaneU32
compare_greater_equal(LaneV3 a, LaneV3 b)
{
    LaneU32 result = {};

    result.v = vcgeq_f32(a.x, b.x) & vcgeq_f32(a.y, b.y) & vcgeq_f32(a.z, b.z);

    return result;
}

force_inline LaneU32
compare_greater(LaneV3 a, LaneV3 b)
{
    LaneU32 result = {};

    result.v = vcgtq_f32(a.x, b.x) & vcgtq_f32(a.y, b.y) & vcgtq_f32(a.z, b.z);

    return result;
}

force_inline LaneU32
compare_less_equal(LaneV3 a, LaneV3 b)
{
    LaneU32 result = {};

    result.v = vcleq_f32(a.x, b.x) & vcleq_f32(a.y, b.y) & vcleq_f32(a.z, b.z);

    return result;
}

force_inline LaneU32
compare_less(LaneV3 a, LaneV3 b)
{
    LaneU32 result = {};

    result.v = vcltq_f32(a.x, b.x) & vcltq_f32(a.y, b.y) & vcltq_f32(a.z, b.z);

    return result;
}

force_inline LaneU32
compare_not_equal(LaneV3 a, LaneV3 b)
{
    LaneU32 result = {};

    result = ~compare_equal(a, b);

    return result;
}

// TODO(joon): reinterpret?

//////////////////// random_series //////////////////// 

struct simd_random_series
{
    LaneU32 next_random;
};

force_inline void
xor_shift_32(LaneU32 *next_random)
{
    LaneU32 x = *next_random;

    x.v = veorq_u32(x.v, vshlq_n_u32(x.v, 13));
    x.v = veorq_u32(x.v, vshrq_n_u32(x.v, 17));
    x.v = veorq_u32(x.v, vshrq_n_u32(x.v, 5));

    *next_random = x;
}

force_inline simd_random_series
start_random_series(u32 seed0, u32 seed1, u32 seed2, u32 seed3)
{
    simd_random_series result = {};
    result.next_random = LaneU32_(seed0, seed1, seed2, seed3);

    return result;
}

force_inline LaneF32
random_between_0_1(simd_random_series *series)
{
    xor_shift_32(&series->next_random);

    LaneF32 max = LaneF32_((r32)U32_Max);

    LaneF32 result = convert_f32_from_u32(series->next_random)/max;

    return result;
}

force_inline LaneF32
random_between(simd_random_series *series, r32 min, r32 max)
{
    LaneF32 simd_min = LaneF32_(min);
    LaneF32 simd_t = LaneF32_(max-min);

    return simd_min + simd_t*random_between_0_1(series);
}

force_inline LaneF32
random_between_minus_1_1(simd_random_series *series)
{
    LaneF32 simd_2 = LaneF32_(2.0f);
    LaneF32 simd_1 = LaneF32_(1.0f);

    LaneF32 result = (simd_2*random_between_0_1(series)) - simd_1;
    return result;
}

#endif // #if ARM

#endif
