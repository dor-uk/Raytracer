#ifndef RAY_H
#define RAY_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "parser.h"
#include "ops.h"

using namespace parser;
using namespace std;

typedef struct Ray {
    Vec3f origin;
    Vec3f direction;
    bool isShadowRay;
    int depth;
} Ray;

inline Vec3f rayEquation(const Ray &r, float t) {
    return {
        r.origin.x + t * r.direction.x,
        r.origin.y + t * r.direction.y,
        r.origin.z + t * r.direction.z
    };
}

inline void set_plane(const Camera &camera, float &left, float &right, float &bottom, float &top) {
    left = camera.near_plane.x;
    right = camera.near_plane.y;
    bottom = camera.near_plane.z;
    top = camera.near_plane.w;
}

inline void calc_center(const Camera &camera, Vec3f &center, const Vec3f &gaze) {
    float near_dist = camera.near_distance;
    center = {camera.position.x + gaze.x * near_dist, 
              camera.position.y + gaze.y * near_dist, 
              camera.position.z + gaze.z * near_dist};
}

inline void calc_top_left(const Camera &camera, Vec3f &top_left, const Vec3f &center, const Vec3f &u, const Vec3f &v, float left, float top) {
    top_left = {center.x + u.x * left + v.x * top, 
                center.y + u.y * left + v.y * top, 
                center.z + u.z * left + v.z * top};
}

inline void calc_s(Vec3f &s, const Vec3f &q, const Vec3f &u, const Vec3f &v, float su, float sv) {
    s = {q.x + u.x * su - v.x * sv,
         q.y + u.y * su - v.y * sv,
         q.z + u.z * su - v.z * sv};
}

Ray produce_ray(const Camera &camera, int i, int j) {
    float left, right, bottom, top, su, sv;
    Ray ray;
    set_plane(camera, left, right, bottom, top);

    Vec3f gaze = normalization(camera.gaze);
    su = (right - left) / camera.image_width * (j + 0.5);
    sv = (top - bottom) / camera.image_height * (i + 0.5);

    Vec3f m;
    calc_center(camera, m, gaze);

    Vec3f u = normalization(cross_product(gaze, camera.up)); 
    Vec3f v = cross_product(u, gaze); 

    Vec3f q;
    calc_top_left(camera, q, m, u, v, left, top);

    Vec3f s;
    calc_s(s, q, u, v, su, sv);

    ray.origin = camera.position;
    ray.direction = normalization(subtraction(s, camera.position));
    ray.depth = 0;
    ray.isShadowRay = false;

    return ray;
}

#endif