#ifndef SHADE_H
#define SHADE_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "parser.h"
#include "ops.h"
#include "ray.h"
#include "hit_record.h"

using namespace parser;
using namespace std;

inline Vec3f calc_irradiance(const PointLight &point_light, const Vec3f &intersect_point) {
    Vec3f d = subtraction(point_light.position, intersect_point);
    float d_square = dot_product(d, d);

    if (d_square != 0.0) {
        return {point_light.intensity.x / d_square,
                point_light.intensity.y / d_square,
                point_light.intensity.z / d_square};
    }
    return {0, 0, 0}; 
}

const Vec3f calc_diffuse(const PointLight &curr_light, const Scene &scene, int matID, const Vec3f &normal, const Vec3f &intersect_point) {
    Vec3f irradiance, wi;
    float dot_product_result;
    const Material &material = scene.materials[matID - 1];

    irradiance = calc_irradiance(curr_light, intersect_point);
    wi = normalization(subtraction(curr_light.position, intersect_point));
    dot_product_result = fmax(dot_product(wi, normal), 0.0f); 

    return {material.diffuse.x * dot_product_result * irradiance.x,
            material.diffuse.y * dot_product_result * irradiance.y,
            material.diffuse.z * dot_product_result * irradiance.z};
}

Vec3f calc_specular(const PointLight &curr_light, const Scene &scene, const Ray &ray, int matID, const Vec3f &normal, const Vec3f &intersect_point) {
    Vec3f irradiance, wi, wo, half_vector;
    float dot_product_result;
    const Material &material = scene.materials[matID - 1];

    irradiance = calc_irradiance(curr_light, intersect_point);
    wi = normalization(subtraction(curr_light.position, intersect_point));
    wo = {-ray.direction.x, -ray.direction.y, -ray.direction.z}; 

    half_vector = normalization(addition(wi, wo)); 
    dot_product_result = fmax(dot_product(normal, half_vector), 0.0f);

    return {material.specular.x * pow(dot_product_result, material.phong_exponent) * irradiance.x,
            material.specular.y * pow(dot_product_result, material.phong_exponent) * irradiance.y,
            material.specular.z * pow(dot_product_result, material.phong_exponent) * irradiance.z};
}

inline bool is_mirror(const Scene &scene, int matID) {
    const Vec3f &mirror = scene.materials[matID - 1].mirror; 
    return (mirror.x > 0 || mirror.y > 0 || mirror.z > 0);
}

#endif