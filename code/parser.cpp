#include "parser.h"

// NOTE(joon) 'peek' function takes a Tokenizer as a value, therefore 
// not advancing the tokenizer
// 'eat' function takes a tokenizer pointer, and will advance the tokenizer
inline void
eat(Tokenizer *tokenizer, i32 a)
{
    tokenizer->at += a;
}

// NOTE(joon) this function has no bound checking
internal b32
string_compare(char *a, char *b)
{
    b32 result = true;

    while((*a != '\0') && 
        ((u8)*b != '\0'))
    {
        if(*a != (u8)*b)
        {
            result = false;
            break;
        }

        a++;
        b++;
    }

    return result;
}

// NOTE(joon) this function has no bound check
internal void
unsafe_string_append(char *dest, const char *source, u32 source_size)
{
    while(*dest != '\0')
    {
        dest++;
    }

    while(source_size-- > 0)
    {
        *dest++ = *source++;
    }
}

// NOTE(joon) this function has no bound check
internal void
unsafe_string_append(char *dest, const char *source)
{
    while(*dest != '\0')
    {
        dest++;
    }
    
    while(*source != '\0')
    {
        *dest++ = *source++;
    }
}

internal void
unsafe_string_append(char *dest, u32 number)
{
    b32 most_significant_number_appeared = false;
    // TODO(joon) complete hack!!!
    for(u32 denom = 100;
            denom > 0;
            denom /= 10)

    {
        u32 nom = (number / denom) % 10;
        if(nom != 0)
        {
            most_significant_number_appeared = true;
        }

        if(most_significant_number_appeared)
        {
            while(*dest != '\0')
            {
                dest++;
            }
            *dest = nom + '0';
        }
    }
}

// NOTE(joon) THIS DOES NOT WORK IF THERE WAS A DIRECTORY WITH .
internal void
get_extension(char *dest, char *source)
{
    assert(*dest == 0);

    while(*source != '.')
    {
        source++;
    }

    source++;

    while(*source != '\0')
    {
        *dest++ = *source++;
    }
}

internal void
eat_until_whitespace(Tokenizer *tokenizer)
{
    while(tokenizer->at < tokenizer->one_past_end)
    {
        if(*tokenizer->at == '\n' ||
            *tokenizer->at == '\r' ||
           *tokenizer->at == ' ')
        {
            break;
        }

        tokenizer->at++;
    }
}

internal void
eat_until_newline(Tokenizer *tokenizer)
{
    while(tokenizer->at != tokenizer->one_past_end)
    {
        if(*tokenizer->at == '\n' ||
            *tokenizer->at == '\r')
        {
            break;
        }

        tokenizer->at++;
    }
}


internal void
eat_all_whitespaces(Tokenizer *tokenizer)
{
    while(tokenizer->at != tokenizer->one_past_end)
    {
        if(!(*tokenizer->at == '\n' ||
            *tokenizer->at == '\r' ||
           *tokenizer->at == ' '))
        {
            break;
        }

        tokenizer->at++;
    }
}

internal ParseNumericResult
eat_numeric(Tokenizer *tokenizer)
{
#if 0
    assert(*start == '0' || 
            *start == '1' || 
            *start == '2' || 
            *start == '3' || 
            *start == '4' || 
            *start == '5' || 
            *start == '6' || 
            *start == '7' || 
            *start == '8' || 
            *start == '9');
#endif

    ParseNumericResult result = {};

    b32 scientific_notation = false;
    f64 decimal_point_adjustment = 10.0f;
    f64 number = 0;
    while(tokenizer->at < tokenizer->one_past_end)
    {
        if(*tokenizer->at >= '0' && *tokenizer->at <= '9' )
        {
            number *= 10;
            number += *tokenizer->at-'0';
        }
        else if(*tokenizer->at == '.')
        {
            result.is_float = true;
        }
        else if(*tokenizer->at == 'e')
        {
            // -5.4335527188698052e-09
            scientific_notation = true;
            break;
        }
        else
        {
            break;
        }

        if(result.is_float)
        {
            decimal_point_adjustment *= 0.1f;
        }

        tokenizer->at++;
    }

    if(result.is_float)
    {
        result.value_f32 = (f32)(number*decimal_point_adjustment);
    }
    else
    {
        result.value_i32 = (i32)number;
    }

    if(scientific_notation)
    {
        assert(*tokenizer->at == 'e');
        eat(tokenizer, 1);

        b32 scientific_plus = false;
        if(*tokenizer->at == '+')
        {
            scientific_plus = true;
        }
        eat(tokenizer, 1);

        f32 scientific_value = 0;
        while(tokenizer->at < tokenizer->one_past_end)
        {
            if(*tokenizer->at >= '0' && *tokenizer->at <= '9' )
            {
                scientific_value *= 10;
                scientific_value += (*tokenizer->at-'0');
            }
            else
            {
                break;
            }

            eat(tokenizer, 1);
        }

        result.value_f32 *= powf(10.0f, (scientific_plus ? 1 : -1) * scientific_value);
    }
    
    return result;
}

internal void
numeric_obj_token_to_f32(ObjToken token, f32 *dest)
{
    if(token.type == obj_token_type_f32)
    {
        *dest = token.value_f32;
    }
    else if(token.type == obj_token_type_i32)
    {
        *dest = (f32)token.value_i32;
    }
    else
    {
        invalid_code_path;
    }
}

internal PlyToken
eat_ply_token(Tokenizer *tokenizer)
{
    PlyToken result = {};

    // NOTE(joon) eat all the unnecessary charcters(whilespace, newline...)
    eat_all_whitespaces(tokenizer);

    if(tokenizer->at < tokenizer->one_past_end)
    {
        if(string_compare((char *)tokenizer->at, "element"))
        {
            result.type = ply_token_type_element;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "vertex"))
        {
            result.type = ply_token_type_vertex;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "face"))
        {
            result.type = ply_token_type_face;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "end_header"))
        {
            result.type = ply_token_type_end_header;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "property"))
        {
            result.type = ply_token_type_property;
            // TODO(joon) might wanna get the rest of the elements after property, too...
            eat_until_newline(tokenizer);
        }
  
        else
        {
            b32 hyphen_appeared = false;
            switch(*tokenizer->at)
            {
                case '-':
                {
                    hyphen_appeared = true;
                    tokenizer->at++;
                };
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                {
                    ParseNumericResult parse_numeric_result = eat_numeric(tokenizer);

                    if(parse_numeric_result.is_float)
                    {
                        result.type = ply_token_type_f32;
                        result.value_f32 = parse_numeric_result.value_f32;
                        result.is_float = true;

                        if(hyphen_appeared)
                        {
                            result.value_f32 *= -1.0f;
                            hyphen_appeared = false;
                        }
                    }
                    else
                    {
                        result.type = ply_token_type_i32;
                        result.value_i32 = parse_numeric_result.value_i32;
                        result.is_float = false;
                        if(hyphen_appeared)
                        {
                            result.value_i32 *= -1;
                            hyphen_appeared = false;
                        }
                    }
                }break;

                default:
                {
                    eat_until_newline(tokenizer);
                }break;
            }

        }

    }

    return result;
}

internal PlyToken
eat_and_check_ply_token(Tokenizer *tokenizer, PlyTokenType expected_type)
{
    PlyToken result = eat_ply_token(tokenizer);
    assert(result.type == expected_type);

    return result;
}

internal PlyToken
peek_ply_token(Tokenizer tokenizer)
{
    PlyToken result = eat_ply_token(&tokenizer);

    return result;
}

internal ParsePlyHeaderResult
parse_ply_header(u8 *memory, u32 file_size)
{
    ParsePlyHeaderResult result = {};

    Tokenizer tokenizer = {};
    tokenizer.at = memory;
    tokenizer.one_past_end = memory + file_size;

    b32 end_header_appeared = false;
    while(tokenizer.at < tokenizer.one_past_end && !end_header_appeared)
    {
        PlyToken token = eat_ply_token(&tokenizer);

        switch(token.type)
        {
            case ply_token_type_element:
            {
                PlyToken element_name = eat_ply_token(&tokenizer);
                if(element_name.type == ply_token_type_vertex)
                {
                    PlyToken vertex_count = eat_and_check_ply_token(&tokenizer, ply_token_type_i32);
                    result.vertex_count = vertex_count.value_i32;
                }
                else if(element_name.type == ply_token_type_face)
                {
                    // It seems like this information is useless?
                }
                else
                {
                    invalid_code_path;
                }


            }break;

            case ply_token_type_property:
            {
                result.vertex_property_count++;
            }break;

            case ply_token_type_end_header:
            {
                end_header_appeared = true;
            }break;
        }
    }

    // TODO(joon) hack, probably not a bad idea (but what if .ply file only has two elements?)
    result.vertex_property_count--;

    // NOTE(joon) ply files do not specify how many indices are there, so we need to get them ourselvese
    for(u32 vertex_index = 0;
            vertex_index < result.vertex_count;
            ++vertex_index)
    {
        for(u32 property_index = 0;
                property_index < result.vertex_property_count;
                ++property_index)
        {
            PlyToken token = eat_ply_token(&tokenizer);

            assert(token.type == ply_token_type_f32 ||
                   token.type == ply_token_type_i32);

        }
    }

    while(tokenizer.at < tokenizer.one_past_end && 
            peek_ply_token(tokenizer).type != 0) // For eof
    {
        PlyToken index_count = eat_ply_token(&tokenizer);
        assert(index_count.type == ply_token_type_i32 &&
                index_count.value_i32 >= 3);

        result.index_count += (u32)(3 * (index_count.value_i32 - 2));

        eat_until_newline(&tokenizer);
    }
    // NOTE(joon) now get how many indices do we need

    return result;
}

// NOTE(joon) minimal ply parser, that only parses vertices for now. 
internal void
parse_ply(ParsePlyHeaderResult *header, u8 *memory, u32 file_size, v3 *vertices, u32 *indices)
{
    Tokenizer tokenizer = {};
    tokenizer.at = memory;
    tokenizer.one_past_end = memory + file_size;

    while(tokenizer.at < tokenizer.one_past_end)
    {
        PlyToken token = eat_ply_token(&tokenizer);

        if(token.type == ply_token_type_end_header)
        {
            break;
        }
    }

    u32 vertex_index = 0;
    for(u32 i = 0;
            i < header->vertex_count;
            ++i)
    {
        // TODO(joon) This does not work, if the position does not come first
        PlyToken p_x = eat_ply_token(&tokenizer);
        PlyToken p_y = eat_ply_token(&tokenizer);
        PlyToken p_z = eat_ply_token(&tokenizer);
        assert(p_x.type == ply_token_type_i32 ||
                p_x.type == ply_token_type_f32);
        assert(p_y.type == ply_token_type_i32 ||
                p_y.type == ply_token_type_f32);
        assert(p_z.type == ply_token_type_i32 ||
                p_z.type == ply_token_type_f32);

        if(p_x.type == ply_token_type_i32)
        {
            vertices[i].x = p_x.value_i32;
        }
        else if(p_x.type == ply_token_type_f32)
        {
            vertices[i].x = p_x.value_f32;
        }

        if(p_y.type == ply_token_type_i32)
        {
            vertices[i].y = p_y.value_i32;
        }
        else if(p_y.type == ply_token_type_f32)
        {
            vertices[i].y = p_y.value_f32;
        }

        if(p_z.type == ply_token_type_i32)
        {
            vertices[i].z = p_z.value_i32;
        }
        else if(p_z.type == ply_token_type_f32)
        {
            vertices[i].z = p_z.value_f32;
        }
        
        vertex_index++;
        eat_until_newline(&tokenizer);
    }

    assert(header->vertex_count == vertex_index);

    u32 index_index = 0;
    // NOTE(joon) this assumes that the indices will always appear at the last
    while(tokenizer.at < tokenizer.one_past_end && 
            peek_ply_token(tokenizer).type != 0) // for eof
    {
        PlyToken index_count = eat_ply_token(&tokenizer);
        assert(index_count.type == ply_token_type_i32 &&
                index_count.value_i32 >= 3);

        PlyToken index_0 = eat_ply_token(&tokenizer); // will be used as a first index to the strip for this line
        PlyToken index_1 = eat_ply_token(&tokenizer);
        PlyToken index_2 = eat_ply_token(&tokenizer);

        assert(index_0.type == ply_token_type_i32 && 
               index_1.type == ply_token_type_i32 &&
               index_2.type == ply_token_type_i32);

        indices[index_index++] = index_0.value_i32;
        indices[index_index++] = index_1.value_i32;
        indices[index_index++] = index_2.value_i32;

        // starting from 1, as we already parsed the first strip
        for(u32 strip_index = 1;
                strip_index < index_count.value_i32 - 2;
                ++strip_index)
        {
            u32 second_index = indices[index_index-1];
            
            indices[index_index++] = index_0.value_i32;
            indices[index_index++] = second_index;
            indices[index_index++] = eat_ply_token(&tokenizer).value_i32;
        }
    }

    assert(index_index == header->index_count);
}

// NOTE/Joon: This function is more like a general purpose token getter, with minimum erro checking.
// the error checking itself will happen inside the parsing loop, not here.
internal ObjToken
eat_obj_token(Tokenizer *tokenizer)
{
    ObjToken result = {};

    eat_all_whitespaces(tokenizer);
    
    if(tokenizer->at < tokenizer->one_past_end)
    {
        b32 hyphen_appeared = false;
        if(*tokenizer->at ==  '-')
        {
            // NOTE(joon) rare case where the case should go on, instead of breaking out
            hyphen_appeared = true;
            tokenizer->at++;
        }

        if(string_compare((char *)tokenizer->at, "v "))
        {
            result.type = obj_token_type_v;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "vt "))
        {
            result.type = obj_token_type_vt;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "vn "))
        {
            result.type = obj_token_type_vn;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "f "))
        {
            result.type = obj_token_type_f;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "mtllib "))
        {
            eat_until_newline(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "o "))
        {
            eat_until_newline(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "usemtl "))
        {
            eat_until_newline(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "g "))
        {
            eat_until_newline(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "body"))
        {
            eat_until_newline(tokenizer);
        }

        else if(*tokenizer->at == '/')
        {
            result.type = obj_token_type_slash;
            eat(tokenizer, 1);
        }
        else if(*tokenizer->at == '#')
        {
            eat_until_newline(tokenizer);
        }
        else if((*tokenizer->at >= '0') && 
            (*tokenizer->at <= '9'))
        {
            ParseNumericResult parse_result = eat_numeric(tokenizer);

            if(parse_result.is_float)
            {
                result.type = obj_token_type_f32;
                result.value_f32 = parse_result.value_f32;

                if(hyphen_appeared)
                {
                    result.value_f32 *= -1.0f;
                    hyphen_appeared = false;
                }
            }

            else
            {
                result.type = obj_token_type_i32;
                result.value_i32 = parse_result.value_i32;

                if(hyphen_appeared)
                {
                    result.value_i32 *= -1;
                    hyphen_appeared = false;
                }
            }

        }
    }

    return result;
}

internal ObjToken
peek_obj_token(Tokenizer tokenizer)
{
    ObjToken result = eat_obj_token(&tokenizer);

    return result;
}

// TODO(joon): parsing positions and vertex normals work just fine,
// but havent yet tested with the texture coords. Will do when I have a png loader :)
// pre_parse returns how many vertices / normals / indices the user needs to allocate.
internal PreParseObjResult
pre_parse_obj(u8 *file, size_t file_size)
{
    assert(file && file_size > 0);

    PreParseObjResult result = {};

    Tokenizer tokenizer = {};
    tokenizer.at = file;
    tokenizer.one_past_end = file + file_size;

    b32 v_appeared = false;
    b32 vn_appeared = false;
    b32 vt_appeared = false;

    while(tokenizer.at < tokenizer.one_past_end)
    {
        ObjToken token = eat_obj_token(&tokenizer);
        switch(token.type)
        {
            case obj_token_type_v:
            {
                result.position_count++;
                v_appeared = true;
            }break;
            case obj_token_type_vn:
            {
                result.normal_count++;
                vn_appeared = true;
            }break;
            case obj_token_type_vt:
            {
                result.texcoord_count++;
                vt_appeared = true;
            }break;

            case obj_token_type_f:
            {
                u32 number_count = 0;
                while(1)
                {
                    ObjToken t = peek_obj_token(tokenizer);

                    if(t.type == obj_token_type_i32)
                    {
                        number_count++;
                        eat_obj_token(&tokenizer); // actually advance
                    }
                    else if(t.type == obj_token_type_slash)
                    {
                        eat_obj_token(&tokenizer); // actually advance
                    }
                    else
                    {
                        break;
                    }
                }

                assert(v_appeared + vt_appeared + vn_appeared != 0);

                u32 strip_count = (number_count / (v_appeared + vt_appeared + vn_appeared));
                result.index_count += 3 * (strip_count - 2);
            }break;
        }
    }

    if(v_appeared)
    {
        if(!vt_appeared && !vn_appeared)
        {
            result.vertex_type = obj_vertex_type_v;
        }
        else if(!vt_appeared && vn_appeared)
        {
            result.vertex_type = obj_vertex_type_v_vn;
        }
        else if(vt_appeared && !vn_appeared)
        {
            result.vertex_type = obj_vertex_type_v_vt;
        }
        else if(vt_appeared && vn_appeared)
        {
            result.vertex_type = obj_vertex_type_v_vt_vn;
        }
    }
    else
    {
        invalid_code_path;
    }

    return result;
}

internal void
parse_obj(PreParseObjResult *pre_parse, u8 *file, u32 file_size, 
        v3 *positions, v3 *normals, v2 *texcoords, u32 *indices)
{
    assert(file && file_size > 0);

    Tokenizer tokenizer = {};
    tokenizer.at = file;
    tokenizer.one_past_end = file + file_size;

    u32 position_index = 0;
    u32 normal_index = 0;
    u32 index_index = 0;
    while(tokenizer.at < tokenizer.one_past_end)
    {
        ObjToken token = eat_obj_token(&tokenizer);
        switch(token.type)
        {
            case obj_token_type_v:
            {
                if(positions)
                {
                    ObjToken p0 = eat_obj_token(&tokenizer);
                    ObjToken p1 = eat_obj_token(&tokenizer);
                    ObjToken p2 = eat_obj_token(&tokenizer);

                    v3 *position = positions + position_index++;

                    numeric_obj_token_to_f32(p0, &position->x);
                    numeric_obj_token_to_f32(p1, &position->y);
                    numeric_obj_token_to_f32(p2, &position->z);
                }
            }break;

            case obj_token_type_vn:
            {
                if(normals)
                {
                    ObjToken n0 = eat_obj_token(&tokenizer);
                    ObjToken n1 = eat_obj_token(&tokenizer);
                    ObjToken n2 = eat_obj_token(&tokenizer);

                    v3 *normal = normals + normal_index++;

                    numeric_obj_token_to_f32(n0, &normal->x);
                    numeric_obj_token_to_f32(n1, &normal->y);
                    numeric_obj_token_to_f32(n2, &normal->z);
                }
            }break;

            case obj_token_type_vt:
            {
                // TODO(joon): add support for vt
            }break;

            case obj_token_type_f:
            {
                if(indices)
                {
                    switch(pre_parse->vertex_type)
                    {
                        case obj_vertex_type_v:
                        {
                            ObjToken i0 = eat_obj_token(&tokenizer);
                            ObjToken i1 = eat_obj_token(&tokenizer);
                            ObjToken i2 = eat_obj_token(&tokenizer);
                            indices[index_index++] = i0.value_i32 - 1;
                            indices[index_index++] = i1.value_i32 - 1;
                            indices[index_index++] = i2.value_i32 - 1;

                            while(1)
                            {
                                ObjToken t = peek_obj_token(tokenizer);
                                if(t.type == obj_token_type_i32)
                                {
                                    indices[index_index++] = i0.value_i32 - 1;
                                    indices[index_index] = indices[index_index-2];
                                    index_index++;
                                    indices[index_index++] = t.value_i32 - 1;

                                    eat_obj_token(&tokenizer);
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }break;
                        case obj_vertex_type_v_vn:
                        {
                            ObjToken i0 = eat_obj_token(&tokenizer);
                            ObjToken s0 = eat_obj_token(&tokenizer);
                            ObjToken s1 = eat_obj_token(&tokenizer);
                            ObjToken vni0 = eat_obj_token(&tokenizer);

                            ObjToken i1 = eat_obj_token(&tokenizer);
                            ObjToken s2 = eat_obj_token(&tokenizer);
                            ObjToken s3 = eat_obj_token(&tokenizer);
                            ObjToken vni1 = eat_obj_token(&tokenizer);

                            ObjToken i2 = eat_obj_token(&tokenizer);
                            ObjToken s4 = eat_obj_token(&tokenizer);
                            ObjToken s5 = eat_obj_token(&tokenizer);
                            ObjToken vni2 = eat_obj_token(&tokenizer);

                            indices[index_index++] = i0.value_i32 - 1;
                            indices[index_index++] = i1.value_i32 - 1;
                            indices[index_index++] = i2.value_i32 - 1;

                            //f 304808//304808 304802//304802 304823//304823 

                            u32 number_count = 0;
                            while(1)
                            {
                                ObjToken t = peek_obj_token(tokenizer);
                                if(t.type == obj_token_type_i32)
                                {
                                    if((number_count % 2) == 0)
                                    {
                                        indices[index_index++] = i0.value_i32 - 1;
                                        indices[index_index] = indices[index_index-2];
                                        index_index++;
                                        indices[index_index++] = t.value_i32 - 1;
                                    }

                                    number_count++;

                                    eat_obj_token(&tokenizer);
                                }
                                else if(t.type == obj_token_type_slash)
                                {
                                    eat_obj_token(&tokenizer);
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }break;

                            // TODO(joon) implement these!
                        case obj_vertex_type_v_vt:
                        {
                        }break;

                        case obj_vertex_type_v_vt_vn:
                        {
                            ObjToken i0 = eat_obj_token(&tokenizer);
                            ObjToken s0 = eat_obj_token(&tokenizer);
                            ObjToken vti0 = eat_obj_token(&tokenizer);
                            ObjToken s1 = eat_obj_token(&tokenizer);
                            ObjToken vni0 = eat_obj_token(&tokenizer);

                            ObjToken i1 = eat_obj_token(&tokenizer);
                            ObjToken s2 = eat_obj_token(&tokenizer);
                            ObjToken vti1 = eat_obj_token(&tokenizer);
                            ObjToken s3 = eat_obj_token(&tokenizer);
                            ObjToken vni1 = eat_obj_token(&tokenizer);

                            ObjToken i2 = eat_obj_token(&tokenizer);
                            ObjToken s4 = eat_obj_token(&tokenizer);
                            ObjToken vti2 = eat_obj_token(&tokenizer);
                            ObjToken s5 = eat_obj_token(&tokenizer);
                            ObjToken vni2 = eat_obj_token(&tokenizer);

                            indices[index_index++] = i0.value_i32 - 1;
                            indices[index_index++] = i1.value_i32 - 1;
                            indices[index_index++] = i2.value_i32 - 1;

                            u32 number_count = 0;
                            while(1)
                            {
                                ObjToken t = peek_obj_token(tokenizer);
                                if(t.type == obj_token_type_i32)
                                {
                                    if((number_count % 3) == 0)
                                    {
                                        indices[index_index++] = i0.value_i32 - 1;
                                        indices[index_index] = indices[index_index-2];
                                        index_index++;
                                        indices[index_index++] = t.value_i32 - 1;
                                    }

                                    number_count++;

                                    eat_obj_token(&tokenizer);
                                }
                                else if(t.type == obj_token_type_slash)
                                {
                                    eat_obj_token(&tokenizer);
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }break;
                    }
                }
            }break;
        }
    }
}

#if 1
internal SceneToken
eat_scn_token(Tokenizer *tokenizer)
{
    SceneToken result = {};

    // NOTE(joon) eat all the unnecessary charcters(whilespace, newline...)
    eat_all_whitespaces(tokenizer);

    if(tokenizer->at < tokenizer->one_past_end)
    {
        if(string_compare((char *)tokenizer->at, "screen"))
        {
            result.type = scene_token_type_screen;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "camera"))
        {
            result.type = scene_token_type_camera;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "ambient"))
        {
            result.type = scene_token_type_ambient;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "light"))
        {
            result.type = scene_token_type_light;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "sphere"))
        {
            result.type = scene_token_type_sphere;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "brdf"))
        {
            result.type = scene_token_type_brdf;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "box"))
        {
            result.type = scene_token_type_box;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "cylinder"))
        {
            result.type = scene_token_type_cylinder;
            eat_until_whitespace(tokenizer);
        }
        else if(string_compare((char *)tokenizer->at, "mesh"))
        {
            result.type = scene_token_type_mesh;
            eat_until_whitespace(tokenizer);
        }

        else if(*tokenizer->at == 'b' && *(tokenizer->at + 1) == ' ')
        {
            // TODO(joon) more solid way to differentiate this with the string that has this letter as a starting character? (i.e box?)
            result.type = scene_token_type_b;
            eat_until_whitespace(tokenizer);
        }
        else if(*tokenizer->at == 'q' && *(tokenizer->at + 1) == ' ')
        {
            result.type = scene_token_type_q;
            eat_until_whitespace(tokenizer);
        }
        else if(*tokenizer->at == 'z' && *(tokenizer->at + 1) == ' ')
        {
            result.type = scene_token_type_z;
            eat_until_whitespace(tokenizer);
        }

        // NOTE(joon) This check should happen after the other string checks
        else if((*tokenizer->at >= 'a' && *tokenizer->at < 'z') || 
                (*tokenizer->at >= 'A' && *tokenizer->at < 'Z'))
        {
            // string start
            result.type = scene_token_type_string;
            result.start = tokenizer->at;
            eat_until_whitespace(tokenizer);
            result.value_i32 = tokenizer->at - result.start;
        }

        else
        {
            b32 hyphen_appeared = false;
            switch(*tokenizer->at)
            {
                case '-':
                {
                    hyphen_appeared = true;
                    tokenizer->at++;
                };
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                {
                    ParseNumericResult parse_numeric_result = eat_numeric(tokenizer);

                    if(parse_numeric_result.is_float)
                    {
                        result.type = scene_token_type_f32;
                        result.value_f32 = parse_numeric_result.value_f32;

                        if(hyphen_appeared)
                        {
                            result.value_f32 *= -1.0f;
                            hyphen_appeared = false;
                        }
                    }
                    else
                    {
                        result.type = scene_token_type_i32;
                        result.value_i32 = parse_numeric_result.value_i32;
                        if(hyphen_appeared)
                        {
                            result.value_i32 *= -1;
                            hyphen_appeared = false;
                        }
                    }
                }break;

                default:
                {
                }break;
            }
        }

        eat_until_whitespace(tokenizer);
    }

    return result;
}

internal SceneToken
peek_scn_token(Tokenizer tokenizer)
{
    SceneToken result = eat_scn_token(&tokenizer);

    return result;
}

internal SceneToken
eat_and_check_scn_token(Tokenizer *tokenizer, SceneTokenType expected_type)
{
    SceneToken result = eat_scn_token(tokenizer);
    assert(result.type == expected_type);

    return result;
}

internal void
push_light(ParseSceneResult *scene, void *shape, ShapeType type)
{
    // NOTE(joon) We are saving header(type of the shape) + pointer to the shape
    // save the type
    ShapeType *t = (ShapeType *)push_size(&scene->light_push_buffer, sizeof(ShapeType));
    *t = type;

    switch(type)
    {
        case Shape_Type_Sphere:
        {
            // save the pointer
            Sphere **sphere_ptr = (Sphere **)push_size(&scene->light_push_buffer, sizeof(Sphere *));
            *sphere_ptr = (Sphere *)shape;
        }break;

        case Shape_Type_Mesh:
        {
            invalid_code_path;
            // save the pointer
            Mesh **mesh_ptr = (Mesh **)push_size(&scene->light_push_buffer, sizeof(Mesh *));
            *mesh_ptr = (Mesh *)shape;
        }break;

        case Shape_Type_Cylinder:
        {
            Cylinder **cylinder_ptr = (Cylinder **)push_size(&scene->light_push_buffer, sizeof(Cylinder *));
            *cylinder_ptr = (Cylinder *)shape;
        }break;

        default:
        {
            invalid_code_path;
        }
    }

    scene->light_count++;
}

internal void
parse_scene(ParseSceneResult *scene, u8 *file, u32 file_size, char *base_file_path)
{
    scene->mat_count = 1; // 0th mat is a default mat 

    Tokenizer tokenizer = {};
    tokenizer.at = file;
    tokenizer.one_past_end = file + file_size;
    
    while(tokenizer.at != tokenizer.one_past_end)
    {
        SceneToken token = eat_scn_token(&tokenizer);

        switch(token.type)
        {
            case scene_token_type_screen:
            {
                // syntax: screen width height
                SceneToken output_width = eat_and_check_scn_token(&tokenizer, scene_token_type_i32);
                SceneToken output_height = eat_and_check_scn_token(&tokenizer, scene_token_type_i32);

                scene->output_width = output_width.value_i32;
                scene->output_height = output_height.value_i32;
            }break;
            case scene_token_type_camera:
            {
                // syntax: camera x y z   ry   <orientation spec>
                SceneToken p_x = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken p_y = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken p_z = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken b = eat_and_check_scn_token(&tokenizer, scene_token_type_b);
                SceneToken height_ratio = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken q = eat_and_check_scn_token(&tokenizer, scene_token_type_q);

                SceneToken orientation_w = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken orientation_x = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken orientation_y = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken orientation_z = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);

                scene->camera_p = v3_(p_x.value_f32, p_y.value_f32, p_z.value_f32);
                scene->camera_height_ratio = height_ratio.value_f32;
                scene->camera_quaternion = v4_(orientation_x.value_f32, orientation_y.value_f32, orientation_z.value_f32, orientation_w.value_f32);
            }break;
            case scene_token_type_ambient:
            {
                // syntax: ambient r g b
                SceneToken r = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken g = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken b = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);

                scene->ambient = v3_(r.value_f32, g.value_f32, b.value_f32);
            }break;
            case scene_token_type_light:
            {
                // syntax: light  r g b   
                SceneToken r = eat_and_check_scn_token(&tokenizer, scene_token_type_i32);
                SceneToken g = eat_and_check_scn_token(&tokenizer, scene_token_type_i32);
                SceneToken b = eat_and_check_scn_token(&tokenizer, scene_token_type_i32);
                
                Material *mat = scene->materials + scene->mat_count++;
                zero(mat);
                mat->is_light = true;
                mat->emit_color = v3_(r.value_i32, g.value_i32, b.value_i32);
                assert(scene->mat_count <= array_count(scene->materials));
            }break;
            case scene_token_type_sphere:
            {
                // syntax: sphere x y z   r
                SceneToken x = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken y = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken z = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken r = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);

                Sphere *sphere = scene->spheres + scene->sphere_count++;
                sphere->center = v3_(x.value_f32, y.value_f32, z.value_f32);
                sphere->r = r.value_f32;
                sphere->mat_index = scene->mat_count - 1;

                Material *mat = scene->materials + sphere->mat_index;
                if(mat->is_light)
                {
                    push_light(scene, sphere, Shape_Type_Sphere);
                }

                assert(scene->sphere_count <= array_count(scene->spheres));
            }break;
            case scene_token_type_brdf:
            {
                Material *mat = scene->materials + scene->mat_count++;
                zero(mat);
                
                if(scene->mat_count == 15)
                {
                    int break_here = 0;
                }

                // syntax: brdf  r g b   r g b  alpha
                // later:  brdf  (diffuse )r g b   (spercular)r g b  alpha  refractionr g b ior
                // NOTE(joon) Diffuse Specular Refraction
                SceneToken diffuse_r = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken diffuse_g = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken diffuse_b = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);

                SceneToken specular_r = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken specular_g = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken specular_b = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken specular_a = eat_and_check_scn_token(&tokenizer, scene_token_type_i32);

                SceneToken check = peek_scn_token(tokenizer);
                if(check.type == scene_token_type_f32 || check.type == scene_token_type_i32)
                {
                    // TODO(joon) This will be used when we start to use the proper brdf equation!
                    SceneToken transmission_r = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                    SceneToken transmission_g = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                    SceneToken transmission_b = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                    SceneToken ior = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                    mat->transmission = v3_(transmission_r.value_f32, transmission_g.value_f32, transmission_b.value_f32);
                    mat->ior = ior.value_f32;

                }

                mat->diffuse = v3_(diffuse_r.value_f32, diffuse_g.value_f32, diffuse_b.value_f32);
                mat->specular = v4_(specular_r.value_f32, specular_g.value_f32, specular_b.value_f32, (f32)specular_a.value_i32);

                assert(scene->mat_count <= array_count(scene->materials));
            }break;
            case scene_token_type_box:
            
            {
                // syntax: box bx by bz   dx dy dz
                SceneToken min_x = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken min_y = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken min_z = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken max_x = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken max_y = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken max_z = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);

                AAB *box = scene->boxes + scene->box_count++;
                box->min = v3_(min_x.value_f32, min_y.value_f32, min_z.value_f32);
                box->max = box->min + v3_(max_x.value_f32, max_y.value_f32, max_z.value_f32);
                box->mat_index = scene->mat_count - 1;

                assert(scene->box_count <= array_count(scene->boxes));
            }break;
            case scene_token_type_cylinder:
            {
                // syntax: cylinder bx by bz   ax ay az  r
                SceneToken base_x = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken base_y = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken base_z = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken axis_x = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken axis_y = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken axis_z = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken r = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);

                Cylinder *cylinder = scene->cylinders + scene->cylinder_count++;
                cylinder->base = v3_(base_x.value_f32, base_y.value_f32, base_z.value_f32);
                cylinder->axis = v3_(axis_x.value_f32, axis_y.value_f32, axis_z.value_f32);
                cylinder->r = r.value_f32;
                cylinder->mat_index = scene->mat_count - 1;

                if(cylinder->mat_index)
                {
                    push_light(scene, cylinder, Shape_Type_Cylinder);
                }

                assert(scene->cylinder_count <= array_count(scene->cylinders));
            }break;
            case scene_token_type_mesh:
            {
                MeshInfo *mesh_info = scene->mesh_infos + scene->mesh_count++;

                // TODO(joon) we also need to parse position & 
                // syntax: mesh   filename   tx ty tz   s   <orientation>
                SceneToken file_name = eat_and_check_scn_token(&tokenizer, scene_token_type_string);

                SceneToken translate_x = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken translate_y = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                SceneToken translate_z = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                mesh_info->translate = v3_(translate_x.value_f32, translate_y.value_f32, translate_z.value_f32);

                SceneToken scale = eat_and_check_scn_token(&tokenizer, scene_token_type_f32);
                mesh_info->scale = scale.value_f32;

                SceneToken det = eat_scn_token(&tokenizer);
                if(det.type == scene_token_type_z)
                {
                    SceneToken axis = eat_scn_token(&tokenizer);

                    if(axis.type == scene_token_type_f32)
                    {
                        mesh_info->axis.degree = axis.value_f32;
                    }
                    else if(axis.type == scene_token_type_i32)
                    {
                        mesh_info->axis.degree = axis.value_i32;
                    }

                    SceneToken q = eat_and_check_scn_token(&tokenizer, scene_token_type_q);
                }
                else if(det.type == scene_token_type_q)
                {
                }
                else
                {
                    assert(0);
                }

                SceneToken q_w = eat_scn_token(&tokenizer);
                SceneToken q_x = eat_scn_token(&tokenizer);
                SceneToken q_y = eat_scn_token(&tokenizer);
                SceneToken q_z = eat_scn_token(&tokenizer);

                if(q_x.type == scene_token_type_f32)
                {
                    mesh_info->quaternion.x = q_x.value_f32;
                }
                else if(q_x.type == scene_token_type_i32)
                {
                    mesh_info->quaternion.x = q_x.value_i32;
                }

                if(q_y.type == scene_token_type_f32)
                {
                    mesh_info->quaternion.y = q_y.value_f32;
                }
                else if(q_y.type == scene_token_type_i32)
                {
                    mesh_info->quaternion.y = q_y.value_i32;
                }                

                if(q_z.type == scene_token_type_f32)
                {
                    mesh_info->quaternion.z = q_z.value_f32;
                }
                else if(q_z.type == scene_token_type_i32)
                {
                    mesh_info->quaternion.z = q_z.value_i32;
                }                

                if(q_w.type == scene_token_type_f32)
                {
                    mesh_info->quaternion.w = q_w.value_f32;
                }
                else if(q_w.type == scene_token_type_i32)
                {
                    mesh_info->quaternion.w = q_w.value_i32;
                }

                // TODO(joon) We can just not do it this way, and just store the file start & end
                u32 file_path_size = 256;
                mesh_info->file_path =(char *)malloc(sizeof(char) * file_path_size);
                zero_memory((void *)mesh_info->file_path, file_path_size);
                unsafe_string_append(mesh_info->file_path, base_file_path);
                unsafe_string_append(mesh_info->file_path, (char *)file_name.start, file_name.value_i32);

                mesh_info->mat_index = scene->mat_count - 1;
                
                assert(scene->mesh_count < array_count(scene->mesh_infos));
            }break;
        }
    }
}
#endif

struct PushedString
{
    u8 *start;
    u8 *one_past_end;
};

internal void
parse_obj_path(MemoryArena *arena, u8 *file, u32 file_size)
{
    Tokenizer tokenizer = {};
    tokenizer.at = file;
    tokenizer.one_past_end = file + file_size;

    while(tokenizer.at < tokenizer.one_past_end)
    {
        if((*tokenizer.at >= 'a' && *tokenizer.at <= 'z') || 
           (*tokenizer.at >= 'A' && *tokenizer.at <= 'Z'))
        {
            PushedString *string = push_struct(arena, PushedString);
            string->start = tokenizer.at;
            eat_until_newline(&tokenizer);
            string->one_past_end = tokenizer.at;
        }

        tokenizer.at++;
    }
}


































