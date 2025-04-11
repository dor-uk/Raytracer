#include <iostream>
#include <thread> 
#include <vector> 
#include "parser.h"
#include "ppm.h"
#include "ops.h"
#include "ray.h"
#include "hit_record.h"
#include "shading.h"
#include "bvh.h"

using namespace parser;
using namespace std;

typedef unsigned char RGB[3];


HitRecord closestHit(const Scene &scene, const Ray &ray) {
    int nos, nott, nom;
    int triNo, meshNo, sphereNo;
   
    nott = scene.triangles.size();
    nom = scene.meshes.size();
    nos = scene.spheres.size();
    
    HitRecord closest_hit;
    closest_hit.is_hit = false;  

    for (triNo = 0; triNo < nott; ++triNo) {
        const Triangle &curr_tri = scene.triangles[triNo];
        Vec3f v0, v1, v2;  
        v0 = scene.vertex_data[curr_tri.indices.v0_id - 1];
        v1 = scene.vertex_data[curr_tri.indices.v1_id - 1];
        v2 = scene.vertex_data[curr_tri.indices.v2_id - 1];

        HitRecord hit_rec = tri_intersection(ray, v0, v1, v2, curr_tri.material_id, triNo);

        if (hit_rec.is_hit && hit_rec.t >= 0 && (!closest_hit.is_hit || hit_rec.t < closest_hit.t)) {
            closest_hit = hit_rec;
        }
    }

    for (meshNo = 0; meshNo < nom; ++meshNo) {
        const Mesh &curr_mesh = scene.meshes[meshNo]; 

        HitRecord hit_rec = meshIntersection(ray, curr_mesh, scene, curr_mesh.material_id, meshNo);

        if (hit_rec.is_hit && hit_rec.t >= 0 && (!closest_hit.is_hit || hit_rec.t < closest_hit.t)) {
            closest_hit = hit_rec;
        }
    }

    for (sphereNo = 0; sphereNo < nos; ++sphereNo) {
        const Sphere &curr_sphere = scene.spheres[sphereNo];  
        Vec3f center = scene.vertex_data[curr_sphere.center_vertex_id - 1];
        float radius = curr_sphere.radius;

        HitRecord hit_rec = sphere_intersection(ray, center, radius, curr_sphere.material_id, sphereNo);

        if (hit_rec.is_hit && hit_rec.t >= 0 && (!closest_hit.is_hit || hit_rec.t < closest_hit.t)) {
            closest_hit = hit_rec;
        }
    }

    return closest_hit; 
}




bool is_shadowed(const Scene& scene, const PointLight& curr_light, const HitRecord &hit_record, int nos, int nott, int nom) {
    Vec3f wi,wi_n, n_epsilon, diff, shadowRayOrigin, normal;
    int tri_num, mesh_num, sphere_num;
    //wi = subtraction(curr_light.position, hit_record.intersectPoint);
    wi_n = normalization(subtraction(curr_light.position, hit_record.intersectPoint));

    normal = normalization(hit_record.normal);
    n_epsilon = {normal.x * scene.shadow_ray_epsilon, normal.y * scene.shadow_ray_epsilon, normal.z * scene.shadow_ray_epsilon};
    shadowRayOrigin = addition(hit_record.intersectPoint, n_epsilon);
    //wi = normalization(subtraction(curr_light.position, shadowRayOrigin));

    Ray shadowRay;
    shadowRay.origin = shadowRayOrigin;
    shadowRay.direction = wi_n;
    shadowRay.isShadowRay = true;

    diff = subtraction(curr_light.position, shadowRay.origin);
    /*float tLight = diff.x / shadowRay.direction.x;
    if (fabs(shadowRay.direction.y) > fabs(shadowRay.direction.x) && fabs(shadowRay.direction.y) > fabs(shadowRay.direction.z)) {
        tLight = diff.y / shadowRay.direction.y;
    } else if (fabs(shadowRay.direction.z) > fabs(shadowRay.direction.x)) {
        tLight = diff.z / shadowRay.direction.z;
    }*/
    float tLight = vector_length(diff);

    for (sphere_num = 0; sphere_num < nos; ++sphere_num) {
        const Sphere &curr_sphere = scene.spheres[sphere_num];
        Vec3f center = scene.vertex_data[curr_sphere.center_vertex_id - 1];
        float radius = curr_sphere.radius;

        HitRecord shadow_hit = sphere_intersection(shadowRay, center, radius, curr_sphere.material_id, sphere_num);

        if (shadow_hit.is_hit && shadow_hit.t >= 0 && shadow_hit.t < tLight) {
            return true; 
        }
    }

    for (tri_num = 0; tri_num < nott; ++tri_num) {
        Vec3f v0, v1, v2;
        const Triangle &curr_tri = scene.triangles[tri_num];
        v0 = scene.vertex_data[curr_tri.indices.v0_id - 1];
        v1 = scene.vertex_data[curr_tri.indices.v1_id - 1];
        v2 = scene.vertex_data[curr_tri.indices.v2_id - 1];

        HitRecord shadow_hit = tri_intersection(shadowRay, v0, v1, v2, curr_tri.material_id, tri_num);

        if (shadow_hit.is_hit && shadow_hit.t >= 0 && shadow_hit.t < tLight) {
            return true; 
        }
    }

    for (mesh_num = 0; mesh_num < nom; ++mesh_num) {
        const Mesh &curr_mesh = scene.meshes[mesh_num];

        HitRecord shadow_hit = meshIntersection(shadowRay, curr_mesh, scene, curr_mesh.material_id, mesh_num);

        if (shadow_hit.is_hit && shadow_hit.t >= 0 && shadow_hit.t < tLight) {
            return true; 
        }
    }

    

    return false; 
}


Vec3f computeColor(const Scene& scene, const HitRecord &hit_record, int max_depth, const Ray &ray, const Camera& camera);

Vec3f apply_shading(const Scene& scene, const HitRecord &hit_record, int max_depth, const Ray &ray, const Camera& camera) {
    Vec3f pixel_color, ambient_light, w_o, n, w_r, reflection_origin;
    float pix_1, pix_2, pix_3;

    int matID = hit_record.material_ID;
    const Material &material = scene.materials[matID - 1];
    ambient_light = scene.ambient_light;

    pix_1 = material.ambient.x * ambient_light.x;
    pix_2 = material.ambient.y * ambient_light.y;
    pix_3 = material.ambient.z * ambient_light.z;

    if (is_mirror(scene, matID)) {
        w_o = {-ray.direction.x, -ray.direction.y, -ray.direction.z};
        n = hit_record.normal;
        float dot_n_w_o = dot_product(n, w_o);

        w_r = { -w_o.x + 2 * n.x * dot_n_w_o,
                      -w_o.y + 2 * n.y * dot_n_w_o,
                      -w_o.z + 2 * n.z * dot_n_w_o };
        w_r = normalization(w_r);

        reflection_origin = addition(hit_record.intersectPoint, {
            w_r.x * scene.shadow_ray_epsilon,
            w_r.y * scene.shadow_ray_epsilon,
            w_r.z * scene.shadow_ray_epsilon
        });

        Ray reflectionRay = {reflection_origin, w_r, false, ray.depth + 1};
        HitRecord reflection_hit = closestHit(scene, reflectionRay);

        if (!(reflection_hit.obj_ID == hit_record.obj_ID && reflection_hit.objectType == hit_record.objectType)) {
            Vec3f reflection = computeColor(scene, reflection_hit, max_depth, reflectionRay, camera);
            pix_1 += reflection.x * material.mirror.x;
            pix_2 += reflection.y * material.mirror.y;
            pix_3 += reflection.z * material.mirror.z;
        }
    }

    for (const PointLight &curr_light : scene.point_lights) {
        if (!is_shadowed(scene, curr_light, hit_record, scene.spheres.size(), scene.triangles.size(), scene.meshes.size())) {
            Vec3f diffusion, specular;
            
            specular = calc_specular(curr_light, scene, ray, matID, hit_record.normal, hit_record.intersectPoint);
            diffusion = calc_diffuse(curr_light, scene, matID, hit_record.normal, hit_record.intersectPoint);

            pix_1 += diffusion.x + specular.x;
            pix_2 += diffusion.y + specular.y;
            pix_3 += diffusion.z + specular.z;
        }
    }

    pixel_color = {pix_1, pix_2, pix_3};
    return pixel_color;
}



Vec3f computeColor(const Scene& scene, const HitRecord &hit_record, int max_depth, const Ray &ray, const Camera& cam){
    Vec3f color;
    Vec3i bg_c;
    if(ray.depth > max_depth){
        return {0, 0, 0};
    }

    if(hit_record.is_hit){
        color = apply_shading(scene, hit_record, max_depth, ray, cam);
        return color;
    }
    else if(ray.depth == 0){
        bg_c = scene.background_color;
        color = {float(bg_c.x) , float(bg_c.y), float(bg_c.z)};
        
        return color;
    }
    else{
        return {0, 0, 0};
    }
}

unsigned char clamp_color(float color_comp) {
    if (color_comp > 255.0f) {
        color_comp = 255.0f;
    } else if (color_comp < 0.0f) {
        color_comp = 0.0f;
    }
    return (unsigned char)(round(color_comp));
}

void render_chunk(const Scene& scene, const Camera& curr_cam, unsigned char* image, int start_row, int end_row, int width, BVHNode* root) {
    int pixel_index, i, j;
    pixel_index = start_row * width * 3;

    for (i = start_row; i < end_row; ++i) {
        for (j = 0; j < width; ++j) {
            Ray ray = produce_ray(curr_cam, i, j);
            HitRecord hitResult = closestHit2(scene, ray, root); 
            Vec3f color = computeColor(scene, hitResult, scene.max_recursion_depth, ray, curr_cam);

            image[pixel_index] = clamp_color(color.x);
            image[pixel_index + 1] = clamp_color(color.y);
            image[pixel_index + 2] = clamp_color(color.z);

            pixel_index += 3;
        }
    }
}

int main(int argc, char* argv[]) {
    parser::Scene scene;
    scene.loadFromXml(argv[1]);

    vector<int> obj_indices;

    int i, triangleOffset, meshOffset, num_cams, cam_no;
    int sphere_size, triangle_size, mesh_size;

    sphere_size = scene.spheres.size();
    triangle_size = scene.triangles.size();
    mesh_size = scene.meshes.size();

    for (i = 0; i < sphere_size; ++i) obj_indices.push_back(i);
    triangleOffset = sphere_size;
    for (i = 0; i < triangle_size; ++i) obj_indices.push_back(triangleOffset + i);
    meshOffset = sphere_size + triangle_size;
    for (i = 0; i < mesh_size; ++i) obj_indices.push_back(meshOffset + i);
    
    BVHNode* root = buildBVH(obj_indices, 0, scene);

    num_cams = scene.cameras.size();

    for (cam_no = 0; cam_no < num_cams; cam_no++) {
        int width, height, num_threads, rows_per_thread;
        const Camera &curr_cam = scene.cameras[cam_no]; 

        height = curr_cam.image_height;
        width = curr_cam.image_width;

        unsigned char* image = new unsigned char[width * height * 3]; 

        num_threads = thread::hardware_concurrency(); 
        rows_per_thread = height / num_threads;

        vector<thread> threads;

        int t;
        for (t = 0; t < num_threads; t++) {
            int start_row, end_row;
            start_row = t * rows_per_thread;
            end_row = (t == num_threads - 1) ? height : (t + 1) * rows_per_thread;

            threads.emplace_back(render_chunk, std::ref(scene), std::ref(curr_cam), image, start_row, end_row, width, root);
        }

        for (auto& th : threads) {
            th.join();
        }

        write_ppm(curr_cam.image_name.c_str(), image, width, height);
        delete[] image;
    }

    deleteBVHTree(root);

    return 0;
}

