# ifndef MEKA_TYPES_H
#define MEKA_TYPES_H

#include <float.h>

typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;
typedef int32_t b32;

typedef uint8_t u8; 
typedef uint16_t u16; 
typedef uint32_t u32;
typedef uint64_t u64;

typedef uintptr_t uintptr;

typedef float r32;
typedef float f32;
#ifndef MEKA_METAL_SHADER
typedef double f64;
#endif

#define Flt_Min FLT_MIN
#define Flt_Max FLT_MAX

#define U8_Max UINT8_MAX
#define U16_Max UINT16_MAX
#define U32_Max UINT32_MAX

#define I32_Min INT32_MIN
#define I32_Max INT32_MAX
#define I16_Min INT16_MIN
#define I16_Max INT16_MAX
#define I8_Min INT8_MIN
#define I8_Max INT8_MAX

#define u8_max UINT8_MAX
#define u16_max UINT16_MAX
#define u32_max UINT32_MAX

#define i32_min INT32_MIN
#define i32_max INT32_MAX
#define i16_min INT16_MIN
#define i16_max INT16_MAX
#define i8_min INT8_MIN
#define i8_max INT8_MAX

#define maximum(a, b) ((a>b)? a:b) 
#define minimum(a, b) ((a<b)? a:b) 
#define sign(a) ((a >= 0.0f) ? 1.0f : -1.0f)

struct v2
{
    r32 x;
    r32 y;
};


struct v3
{
    union
    {
        struct 
        {
            r32 x;
            r32 y;
            r32 z;
        };
        struct 
        {
            r32 r;
            r32 g;
            r32 b;
        };
        struct 
        {
            v2 xy;
            r32 ignored;
        };

        r32 e[3];
    };
};

struct v3u
{
    u32 x;
    u32 y;
    u32 z;
};

struct v4
{
    union
    {
        struct 
        {
            r32 x, y, z, w;
        };

        struct 
        {
            r32 r, g, b, a;
        };
        struct 
        {
            v3 xyz; 
            r32 ignored;
        };
        struct 
        {
            v3 rgb; 
            r32 ignored1;
        };

        r32 e[4];
    };
};

struct m4
{
    union
    {
        struct 
        {
            v4 column[4];
        };
        
        r32 e[16];
    };
};

// NOTE(joon) row major
struct m3
{
    union
    {
        struct 
        {
            /*
               row[0]
               row[1]
               row[2]
            */
            v3 rows[3];
        };

        f32 e[3][3];
    };
};



#endif
