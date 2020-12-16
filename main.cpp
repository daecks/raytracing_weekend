#include <iostream>
#include "sphere.h"
#include "hitable_list.h"
#include "float.h"
#include "camera.h"

vec3 random_in_unit_sphere() {
    vec3 p;
    do {
        p = 2.0*vec3(drand48(), drand48(), drand48()) - vec3(1,1,1);
    } while (p.squared_length() >= 1.0);
    return p;
}

vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2*dot(v,n)*n;
}

class material {
    public:
        virtual bool scatter(const ray& r_in, 
                const hit_record& rec, 
                vec3& attenuation, 
                ray& scattered) const = 0;
};

bool refract(const vec3& v, const vec3& n, float ni_over_nt, vec3& refracted) {
    vec3 uv = unit_vector(v);
    float dt = dot(uv, n);
    float discriminant = 1.0 - ni_over_nt*ni_over_nt*(1-dt*dt);
    if (discriminant > 0) {
        refracted = ni_over_nt*(uv - n*dt) - n*sqrt(discriminant);
        return true;
    }else{
        return false;
    }
}

vec3 color(const ray& r, hitable *world, int depth) {
    hit_record rec;
    if (world->hit(r, 0.001, MAXFLOAT, rec)) {
        ray scattered;
        vec3 attenuation;
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return attenuation*color(scattered, world, depth+1);
        } else {
            return vec3(0,0,0);
        }
    }
    else {
        vec3 unit_direction = unit_vector(r.direction());
        float t = 0.5*(unit_direction.y() + 1.0);
        return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
    }
}

class lambertian : public material {
    public:
        lambertian(const vec3& a) : albedo(a) {}
        virtual bool scatter(const ray& r_in, 
                const hit_record& rec, 
                vec3& attenuation, 
                ray& scattered) const {
            vec3 target = rec.p + rec.normal + random_in_unit_sphere();
            scattered = ray(rec.p, target - rec.p);
            attenuation = albedo;
            return true;
        }

        vec3 albedo;
};

class metal : public material {
    public:
        metal(const vec3& a, float f) : albedo(a) { if (f < 1) fuzz = f; else fuzz = 1; }
        virtual bool scatter(const ray& r_in, 
                const hit_record& rec, 
                vec3& attenuation,
                ray& scattered) const {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
            scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere());
            attenuation = albedo;
            return (dot(scattered.direction(), rec.normal) > 0);
        }
        vec3 albedo;
        float fuzz;
};

float schlick(float cosine, float ref_idx){
    float r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0 * r0;
    return r0 + (1-r0)*pow((1 - cosine), 5);
}

class dielectric : public material {
    public:
        dielectric(float ri) : ref_idx(ri) {}
        virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
            vec3 outward_normal;
            vec3 reflected = reflect(r_in.direction(), rec.normal);
            float ni_over_nt;
            attenuation = vec3(1.0, 1.0, 1.0);
            vec3 refracted;
            float reflect_prob;
            float cosine;
            if (dot(r_in.direction(), rec.normal) > 0) {
                outward_normal = -rec.normal;
                ni_over_nt = ref_idx;
                cosine = ref_idx * dot(r_in.direction(), rec.normal) / r_in.direction().length();
            } else {
                outward_normal = rec.normal;
                ni_over_nt = 1.0 / ref_idx;
                cosine = -dot(r_in.direction(), rec.normal) / r_in.direction().length();
            }

            if ( refract(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
                reflect_prob =schlick(cosine, ref_idx);
            } else {
                scattered = ray(rec.p, reflected);
                reflect_prob = 1.0;
            }

            if(drand48() < reflect_prob) {
                scattered = ray(rec.p, reflected);
            }else{
                scattered = ray(rec.p, refracted);
            }
            return true;

        }
        float ref_idx;
};

hitable *random_scene(){
    static const int num_spheres = 500;
    hitable **list = new hitable*[num_spheres+1]; // +1 for the big sphere (ground), +3 for the example spheres
    static const int something = 11;
    static const float small_sphere_radius = 0.2;

    // Create the sphere upon which everything sits.
    list[0] = new sphere(vec3(0, -10000, 0), 10000, new lambertian(vec3(0.5, 0.5, 0.5)));

    int i = 1;
    for (int a = -something; a < something; a++) {
        for (int b = -something; b < something; b++) {
            float material_choice = drand48();
            vec3 center(a+0.9*drand48(), 0.2, b+0.9*drand48()); // random location for each sphere
            //vec3 center(a, small_sphere_radius, b-2); // fixed location for each sphere
            if ((center - vec3(4, small_sphere_radius, 0)).length() > 0.9) {
                if (material_choice < 0.8) { //diffuse
                    list[i++] = new sphere(center, 
                            small_sphere_radius, 
                            new lambertian(vec3( drand48() * drand48(), drand48() * drand48(), drand48() * drand48())));
                } else if (material_choice < 0.95) { // metal
                    list[i++] = new sphere(center, 
                            small_sphere_radius, 
                            new metal(vec3( 0.5*(1 + drand48()), 0.5*(1 + drand48()), 0.5*(1 + drand48())), 0.5 * drand48()));
                } else { //
                    list[i++] = new sphere(center, small_sphere_radius, new dielectric(1.5));
                }
            }
        }
    }

    list[i++] = new sphere(vec3(0,1,0), 1.0, new dielectric(1.5));
    list[i++] = new sphere(vec3(-4,1,0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
    list[i++] = new sphere(vec3(4,1,0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

    return new hitable_list(list, i);
}


int main(){
    static const int nx = 640;
    static const int ny = 480;
    static const int ns = 1; //TODO: some magic number I need to remind myself about. Higher is better quality.

    static const int vfov = 20; // degrees
    static const float aspect_ratio = float(nx)/float(ny);
    static const float aperture = 0.05;

    hitable *world = random_scene();

    const vec3 lookfrom(13,6,6);
    const vec3 lookat(0,0,0);
    const vec3 vup(0,1,0);
    const float dist_to_focus = (lookfrom - lookat).length();
    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus);
    
    std::cout << "P3\n" << nx << " " << ny << "\n255\n";

    for (int j = ny-1; j >=0; j--) {
        for (int i = 0; i < nx; i++) {
            vec3 col(0, 0, 0);
            for (int s=0; s < ns; s++) {
                float u = float(i + drand48()) / float(nx);
                float v = float(j + drand48()) / float(ny);
                ray r = cam.get_ray(u, v);
                vec3 p = r.point_at_parameter(2.0);
                col += color(r, world, 0);
            }

            col /= float(ns);
            col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
            int ir = int(255.99*col[0]);
            int ig = int(255.99*col[1]);
            int ib = int(255.99*col[2]);

            std::cout << ir << " " << ig << " " << ib << "\n";
        }
    }
}
