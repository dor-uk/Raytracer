#ifndef BVH_H
#define BVH_H

#include <iostream>
#include <algorithm>
#include "parser.h"
#include "ppm.h"
#include "ops.h"
#include "ray.h"
#include "hit_record.h"
#include "shading.h"


const int MAX_BVH_DEPTH = 16;

struct BoundingBox {
    Vec3f min; 
    Vec3f max; 

    BoundingBox() : min{INFINITY, INFINITY, INFINITY}, max{-INFINITY, -INFINITY, -INFINITY} {}
};

inline void expandBoundingBox(BoundingBox &bbox, const Vec3f &point) {
    bbox.min.x = min(bbox.min.x, point.x);
    bbox.min.y = min(bbox.min.y, point.y);
    bbox.min.z = min(bbox.min.z, point.z);

    bbox.max.x = max(bbox.max.x, point.x);
    bbox.max.y = max(bbox.max.y, point.y);
    bbox.max.z = max(bbox.max.z, point.z);
}

inline BoundingBox mergeBoundingBoxes(const BoundingBox &bbox1, const BoundingBox &bbox2) {
    BoundingBox result;
    result.min = {min(bbox1.min.x, bbox2.min.x), min(bbox1.min.y, bbox2.min.y), min(bbox1.min.z, bbox2.min.z)};
    result.max = {max(bbox1.max.x, bbox2.max.x), max(bbox1.max.y, bbox2.max.y), max(bbox1.max.z, bbox2.max.z)};
    return result;
}


struct BVHNode {
    BoundingBox bbox;             
    BVHNode* left = nullptr;       
    BVHNode* right = nullptr;      
    vector<int> objectIndices;     
};

inline BoundingBox calculateBoundingBoxForSphere(const Vec3f &center, float radius) {
    BoundingBox bbox;
    bbox.min = subtraction(center, {radius, radius, radius});
    bbox.max = addition(center, {radius, radius, radius});
    return bbox;
}

inline BoundingBox calculateBoundingBoxForTriangle(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2) {
    BoundingBox bbox;
    
    bbox.min.x = min(v0.x, min(v1.x, v2.x));
    bbox.min.y = min(v0.y, min(v1.y, v2.y));
    bbox.min.z = min(v0.z, min(v1.z, v2.z));

    bbox.max.x = max(v0.x, max(v1.x, v2.x));
    bbox.max.y = max(v0.y, max(v1.y, v2.y));
    bbox.max.z = max(v0.z, max(v1.z, v2.z));

    return bbox;
}

inline BoundingBox calculateBoundingBox(const Vec3f &center, float radius) {
    return calculateBoundingBoxForSphere(center, radius);
}
inline BoundingBox calculateBoundingBox(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2) {
    return calculateBoundingBoxForTriangle(v0, v1, v2);
}

inline BoundingBox calculateBoundingBoxForMesh(const Mesh &mesh, const Scene &scene) {
    BoundingBox bbox;
    for (const auto &face : mesh.faces) {
        Vec3f v0 = scene.vertex_data[face.v0_id - 1];
        Vec3f v1 = scene.vertex_data[face.v1_id - 1];
        Vec3f v2 = scene.vertex_data[face.v2_id - 1];

        BoundingBox triangleBbox = calculateBoundingBox(v0, v1, v2);
        bbox = mergeBoundingBoxes(bbox, triangleBbox);
    }
    return bbox;
}

void deleteBVHTree(BVHNode* node) {
    if (node) {
        deleteBVHTree(node->left);
        deleteBVHTree(node->right);
        delete node;
    }
}

inline bool intersectBoundingBox(const Ray &ray, const BoundingBox &bbox, float &tmin, float &tmax) {
    tmin = (bbox.min.x - ray.origin.x) / ray.direction.x;
    tmax = (bbox.max.x - ray.origin.x) / ray.direction.x;

    if (tmin > tmax) swap(tmin, tmax); 

    float tymin = (bbox.min.y - ray.origin.y) / ray.direction.y;
    float tymax = (bbox.max.y - ray.origin.y) / ray.direction.y;

    if (tymin > tymax) swap(tymin, tymax);

    if ((tmin > tymax) || (tymin > tmax))
        return false;

    if (tymin > tmin)
        tmin = tymin;

    if (tymax < tmax)
        tmax = tymax;

    float tzmin = (bbox.min.z - ray.origin.z) / ray.direction.z;
    float tzmax = (bbox.max.z - ray.origin.z) / ray.direction.z;

    if (tzmin > tzmax) swap(tzmin, tzmax);

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    if (tzmin > tmin)
        tmin = tzmin;

    if (tzmax < tmax)
        tmax = tzmax;

    return true;
}

BVHNode* buildBVH(const vector<int>& objectIndices, int depth, const Scene& scene) {
    if (objectIndices.size() <= 1 || depth > MAX_BVH_DEPTH) {
        BVHNode* leaf = new BVHNode();
        leaf->objectIndices = objectIndices;

        for (int index : objectIndices) {
            BoundingBox objectBbox;

            if (index < scene.spheres.size()) {
                Vec3f center = scene.vertex_data[scene.spheres[index].center_vertex_id - 1];
                objectBbox = calculateBoundingBox(center, scene.spheres[index].radius);
            } 
            else if (index < scene.spheres.size() + scene.triangles.size()) {
                int triangleIndex = index - scene.spheres.size();
                const Triangle &triangle = scene.triangles[triangleIndex];
                Vec3f v0 = scene.vertex_data[triangle.indices.v0_id - 1];
                Vec3f v1 = scene.vertex_data[triangle.indices.v1_id - 1];
                Vec3f v2 = scene.vertex_data[triangle.indices.v2_id - 1];
                objectBbox = calculateBoundingBox(v0, v1, v2);
            } 
            else {
                int meshIndex = index - scene.spheres.size() - scene.triangles.size();
                const Mesh &mesh = scene.meshes[meshIndex];
                objectBbox = calculateBoundingBoxForMesh(mesh, scene);
            }

            leaf->bbox = mergeBoundingBoxes(leaf->bbox, objectBbox);
        }
        return leaf;
    }

    int axis = depth % 3; 
    vector<int> sortedIndices = objectIndices;

    std::sort(sortedIndices.begin(), sortedIndices.end(),
              [&scene, axis](int a, int b) {
                  Vec3f centerA, centerB;

                  if (a < scene.spheres.size()) {
                      centerA = scene.vertex_data[scene.spheres[a].center_vertex_id - 1];
                  } 
                  else if (a < scene.spheres.size() + scene.triangles.size()) {
                      int triangleIndex = a - scene.spheres.size();
                      const Triangle &triangle = scene.triangles[triangleIndex];
                      Vec3f v0 = scene.vertex_data[triangle.indices.v0_id - 1];
                      Vec3f v1 = scene.vertex_data[triangle.indices.v1_id - 1];
                      Vec3f v2 = scene.vertex_data[triangle.indices.v2_id - 1];
                      centerA = {(v0.x + v1.x + v2.x) / 3, (v0.y + v1.y + v2.y) / 3, (v0.z + v1.z + v2.z) / 3};
                  } 
                  else {
                      int meshIndex = a - scene.spheres.size() - scene.triangles.size();
                      const Mesh &mesh = scene.meshes[meshIndex];
                      BoundingBox meshBbox = calculateBoundingBoxForMesh(mesh, scene);
                      centerA = {(meshBbox.min.x + meshBbox.max.x) / 2, 
                                 (meshBbox.min.y + meshBbox.max.y) / 2, 
                                 (meshBbox.min.z + meshBbox.max.z) / 2};
                  }

                  if (b < scene.spheres.size()) {
                      centerB = scene.vertex_data[scene.spheres[b].center_vertex_id - 1];
                  } 
                  else if (b < scene.spheres.size() + scene.triangles.size()) {
                      int triangleIndex = b - scene.spheres.size();
                      const Triangle &triangle = scene.triangles[triangleIndex];
                      Vec3f v0 = scene.vertex_data[triangle.indices.v0_id - 1];
                      Vec3f v1 = scene.vertex_data[triangle.indices.v1_id - 1];
                      Vec3f v2 = scene.vertex_data[triangle.indices.v2_id - 1];
                      centerB = {(v0.x + v1.x + v2.x) / 3, (v0.y + v1.y + v2.y) / 3, (v0.z + v1.z + v2.z) / 3};
                  } 
                  else {
                      int meshIndex = b - scene.spheres.size() - scene.triangles.size();
                      const Mesh &mesh = scene.meshes[meshIndex];
                      BoundingBox meshBbox = calculateBoundingBoxForMesh(mesh, scene);
                      centerB = {(meshBbox.min.x + meshBbox.max.x) / 2, 
                                 (meshBbox.min.y + meshBbox.max.y) / 2, 
                                 (meshBbox.min.z + meshBbox.max.z) / 2};
                  }

                  if (axis == 0) {
                      return centerA.x < centerB.x;
                  } else if (axis == 1) {
                      return centerA.y < centerB.y;
                  } else {
                      return centerA.z < centerB.z;
                  }
              });

    size_t mid = sortedIndices.size() / 2;
    vector<int> leftObjects(sortedIndices.begin(), sortedIndices.begin() + mid);
    vector<int> rightObjects(sortedIndices.begin() + mid, sortedIndices.end());

    BVHNode* node = new BVHNode();
    node->left = buildBVH(leftObjects, depth + 1, scene);
    node->right = buildBVH(rightObjects, depth + 1, scene);

    node->bbox = mergeBoundingBoxes(node->left->bbox, node->right->bbox);

    return node;
}


HitRecord closestHit2(const Scene &scene, const Ray &ray, BVHNode* node) {
    HitRecord closest_hit;
    closest_hit.is_hit = false;

    float tmin, tmax;
    if (!node || !intersectBoundingBox(ray, node->bbox, tmin, tmax)) {
        return closest_hit;
    }

    if (!node->left && !node->right) {
        for (int index : node->objectIndices) {
            HitRecord hit;

            if (index < scene.spheres.size()) {
                const Sphere &sphere = scene.spheres[index];
                Vec3f center = scene.vertex_data[sphere.center_vertex_id - 1];
                hit = sphere_intersection(ray, center, sphere.radius, sphere.material_id, index);
            } 
            else if (index < scene.spheres.size() + scene.triangles.size()) {
                int triangleIndex = index - scene.spheres.size();
                const Triangle &triangle = scene.triangles[triangleIndex];
                Vec3f v0 = scene.vertex_data[triangle.indices.v0_id - 1];
                Vec3f v1 = scene.vertex_data[triangle.indices.v1_id - 1];
                Vec3f v2 = scene.vertex_data[triangle.indices.v2_id - 1];
                hit = tri_intersection(ray, v0, v1, v2, triangle.material_id, triangleIndex);
            } 
            else {
                int meshIndex = index - scene.spheres.size() - scene.triangles.size();
                const Mesh &mesh = scene.meshes[meshIndex];
                hit = meshIntersection(ray, mesh, scene, mesh.material_id, meshIndex);
            }

            if (hit.is_hit && (!closest_hit.is_hit || hit.t < closest_hit.t)) {
                closest_hit = hit;
            }
        }
    } else {
        HitRecord left_hit = closestHit2(scene, ray, node->left);
        HitRecord right_hit = closestHit2(scene, ray, node->right);
        closest_hit = (left_hit.is_hit && (!right_hit.is_hit || left_hit.t < right_hit.t)) ? left_hit : right_hit;
    }

    return closest_hit;
}

#endif