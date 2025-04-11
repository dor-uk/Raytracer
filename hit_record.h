#ifndef HIT_RECORD_H
#define HIT_RECORD_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "parser.h"
#include "ops.h"
#include "ray.h"

using namespace parser;
using namespace std;


enum ObjType {
    SPHERE, TRIANGLE, MESH
};

typedef struct HitRecord {
    ObjType objectType;
    Vec3f normal;
    Vec3f intersectPoint;
    int obj_ID;
    int material_ID;
    float t;
    bool is_hit;
} HitRecord;

inline void sphere_normal_normalization(HitRecord &hit, const Vec3f &center, float r) {
    Vec3f normal = subtraction(hit.intersectPoint, center);
    hit.normal = {normal.x / r, normal.y / r, normal.z / r};
}

inline void assignRecord(HitRecord &hit, ObjType type, int matID, int objID, float t) {
    hit.material_ID = matID;
    hit.obj_ID = objID;
    hit.objectType = type;
    hit.t = t;
    hit.is_hit = true;
}

HitRecord sphere_intersection(const Ray &ray, const Vec3f &center, float r, int matID, int objID) {
    HitRecord hit;
    hit.is_hit = false;

    const float A = dot_product(ray.direction, ray.direction);
    Vec3f sphere_normal = subtraction(ray.origin, center);
    const float B = 2 * dot_product(ray.direction, sphere_normal);
    const float C = dot_product(sphere_normal, sphere_normal) - r * r;
    const float disk = B * B - 4 * A * C;

    if (disk >= 0) {
        float t1, t2, t;
        const float sqrt_disc = sqrtf(disk);
        const float inv_2A = 0.5f / A;
        t1 = (-B - sqrt_disc) * inv_2A;
        t2 = (-B + sqrt_disc) * inv_2A;
        t = (t1 > 0) ? t1 : ((t2 > 0) ? t2 : -1);

        if (t > 0) {
            assignRecord(hit, SPHERE, matID, objID, t);
            hit.intersectPoint = rayEquation(ray, t);
            sphere_normal_normalization(hit, center, r);
        }
    }

    return hit;
}

HitRecord tri_intersection(const Ray &ray, const Vec3f &a, const Vec3f &b, const Vec3f &c, int matID, int objID)
{
    Vec3f o, d, a_b, a_c, a_o, normal;
    float dot_n_ray, det_A, t , gamma, beta;
	HitRecord hit_rec;
	hit_rec.is_hit = false;

	d = ray.direction;
    o = ray.origin;

    a_c = subtraction(a, c);
	a_b = subtraction(a, b);
	a_o = subtraction(a, o);

	normal = normalization(cross_product(a_b, a_c));

    dot_n_ray = dot_product(normal, ray.direction);
	if (dot_n_ray > 0) {return hit_rec;}

	det_A = find_determinant(a_b, a_c, d);
	if(det_A == 0.0){return hit_rec;}

	t = (find_determinant(a_b, a_c, a_o)) / det_A;
	if(t <= 0.0) {return hit_rec;}

	gamma = (find_determinant(a_b, a_o, d)) / det_A;  
    if(gamma > 1 || gamma < 0) {return hit_rec;}

	beta = (find_determinant(a_o, a_c, d)) / det_A;
    if(beta > (1 - gamma) || beta < 0) {return hit_rec;}

	assignRecord(hit_rec, TRIANGLE, matID, objID, t);
    hit_rec.intersectPoint = rayEquation(ray, t);
	hit_rec.normal = normalization(cross_product(subtraction(b, a), subtraction(c, a)));

	return hit_rec;
}

HitRecord meshIntersection(const Ray &ray, const Mesh &mesh, const Scene &scene, int matID, int objID) {
    HitRecord closest_hit;
    Vec3f v0, v1, v2, edge1, edge2, normal;
    int face_num, mesh_size;
    closest_hit.is_hit = false; 
    float dot_n_ray;

    mesh_size = mesh.faces.size();

    for (face_num = 0; face_num < mesh_size; face_num++) {
        const Face &face = mesh.faces[face_num]; 

        v0 = scene.vertex_data[face.v0_id - 1];
        v1 = scene.vertex_data[face.v1_id - 1];
        v2 = scene.vertex_data[face.v2_id - 1];

        edge1 = subtraction(v0, v1);
        edge2 = subtraction(v0, v2);
        

        normal = normalization(cross_product(edge1, edge2));

        dot_n_ray = dot_product(normal, ray.direction);
        if (dot_n_ray > 0) {continue;}

        HitRecord hit = tri_intersection(ray, v0, v1, v2, matID, objID);

        if (hit.is_hit && hit.t >= 0) {
            if (!closest_hit.is_hit || hit.t < closest_hit.t) {
                closest_hit = hit; 
            }
        }
    }

    return closest_hit; 
}


#endif
