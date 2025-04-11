#ifndef OPS_H
#define OPS_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "parser.h"

using namespace parser;
using namespace std;

inline Vec3f addition(const Vec3f &v1, const Vec3f &v2) {
    return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

inline Vec3f subtraction(const Vec3f &v1, const Vec3f &v2) {
    return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

inline float dot_product(const Vec3f &v1, const Vec3f &v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline Vec3f cross_product(const Vec3f &v1, const Vec3f &v2) {
    return {
        v1.y * v2.z - v1.z * v2.y,
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x
    };
}

inline float vector_length(const Vec3f &v) {
    return sqrt(dot_product(v, v)); 
}

inline float vector_distance(const Vec3f &v1, const Vec3f &v2) {
    Vec3f diff = subtraction(v1, v2);
    return vector_length(diff);
}

inline Vec3f normalization(const Vec3f &v) {
    float length = vector_length(v);
    return {v.x / length, v.y / length, v.z / length};
}

inline float find_determinant(const Vec3f &v1, const Vec3f &v2, const Vec3f &v3) {
    return v1.x * (v2.y * v3.z - v3.y * v2.z)
         - v1.y * (v2.x * v3.z - v3.x * v2.z)
         + v1.z * (v2.x * v3.y - v3.x * v2.y);
}

#endif