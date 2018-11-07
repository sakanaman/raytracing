#ifndef FIGURE_H
#define FIGURE_H

#include "vec3.h"
#include "texture.h"
#include "ray.h"
#include "hit.h"
#include <cmath>
#include <memory>
#include <utility>

class Figure
{
  public:
    int material;
    std::shared_ptr<Texture> texture;
    Figure(){};
    Figure(std::shared_ptr<Texture> texture, int material) : texture(texture), material(material){};
    virtual bool intersect(const Ray &ray, Hit &hit) const = 0;
};

class Sphere : public Figure
{
  public:
    Vec3 center;
    double radius;
    Sphere(const Vec3 &center, double radius, std::shared_ptr<Texture> texture, int material)
        : center(center), radius(radius), Figure(texture, material){};

    virtual bool intersect(const Ray &ray, Hit &hit) const
    {
        double d_norm = ray.direction.length();
        double oc_norm = (ray.origin - center).length();

        double a = d_norm * d_norm;
        double b = 2 * dot(ray.direction, ray.origin - center);
        double c = oc_norm * oc_norm - radius * radius;
        double d = b * b - 4 * a * c;
        if (d < 0)
            return false;

        double t1 = (-b - sqrt(d)) / (2 * a);
        double t2 = (-b + sqrt(d)) / (2 * a);

        double t = t1;
        if (t < 0)
        {
            t = t2;
            if (t < 0)
                return false;
        }

        hit.t = t;
        hit.hitPos = ray.origin + t * ray.direction;
        hit.hitNormal = normalize(hit.hitPos - center);
        hit.hitShape = this;
        double phi = std::atan2(hit.hitNormal.z, hit.hitNormal.x) - M_PI / 2.0;
        if (phi < 0)
            phi += 2 * M_PI;
        if (phi > 2 * M_PI)
            phi -= 2 * M_PI;
        double theta = std::acos(hit.hitNormal.y);
        hit.u = phi / (2 * M_PI);
        hit.v = theta / M_PI;
        return true;
    };
};

void renritsu_3(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &origin, const Vec3 &direction, double &m, double &n, double &t)
{
    Vec3 E1 = b - a;
    Vec3 E2 = c - a;
    Vec3 T = origin - a;
    Vec3 R = direction;

    Vec3 P = cross(R, E2);
    Vec3 Q = cross(T, E1);

    t = dot(Q, E2) / dot(P, E1);
    m = dot(P, T) / dot(P, E1);
    n = dot(Q, R) / dot(P, E1);
}

class Tri : public Figure
{
  public:
    Vec3 a;
    Vec3 b;
    Vec3 c;

    Tri(const Vec3 &a, const Vec3 &b, const Vec3 &c, std::shared_ptr<Texture> texture, int material) : a(a), b(b), c(c), Figure(texture, material){};

    virtual bool intersect(const Ray &ray, Hit &hit) const
    {
        double m, n, t;
        renritsu_3(a, b, c, ray.origin, ray.direction, m, n, t);
        if (m + n <= 1 && 0 <= m && m <= 1 && n <= 1 && n >= 0 && t >= 0)
        {
            hit.t = t;
            hit.hitPos = ray.origin + t * ray.direction;
            hit.hitNormal = dot(cross(b - a, c - a), ray.direction) < 0 ? normalize(cross(b - a, c - a)) : -1 * normalize(cross(b - a, c - a));

            hit.hitShape = this;
            hit.u = n / (m + n);
            hit.v = m + n;

            return true;
        }
        return false;
    }
};

class Cylinder : public Figure
{
  public:
    double radius, phiMax;
    Vec3 zMin, zMax;
    Vec3 xv, n, zv;
    Vec3 x, y, z;

    Cylinder(const double radius, const double phimax, const Vec3 &zMax, const Vec3 &zMin,
             std::shared_ptr<Texture> texture, int material)
        : radius(radius), phiMax(phiMax), zMax(zMax), zMin(zMin), Figure(texture, material)
    {
        n = normalize(zMax - zMin);
        coordinate_system(n, xv, zv);

        x = xv;
        y = n;
        z = zv;

        Transpose3x3(xv, n, zv);
    };

    virtual bool intersect(const Ray &ray, Hit &hit) const
    {

        Vec3 new_origin = ((ray.origin - zMin).x) * xv + ((ray.origin - zMin).y) * n + ((ray.origin - zMin).z) * zv;
        Vec3 new_direction = (ray.direction.x) * xv + (ray.direction.y) * n + (ray.direction.z) * zv;
        double z_max = (zMax - zMin).length();
        double z_min = 0;
        double a = (new_direction.x) * (new_direction.x) + (new_direction.z) * (new_direction.z);
        double b = 2 * ((new_direction.x) * (new_origin.x) + (new_direction.z) * (new_origin.z));
        double c = (new_origin.x) * (new_origin.x) + (new_origin.z) * (new_origin.z) - radius * radius;
        double d = b * b - 4 * a * c;

        if (d < 0)
            return false;

        double q;
        if (b < 0)
            q = (-b + sqrt(d)) / 2.0;
        else
            q = (-b - sqrt(d)) / 2.0;

        double t0 = q / a;
        double t1 = c / q;
        if (t0 > t1)
            std::swap(t0, t1);

        double t = t0;
        if (t0 < 0)
        {
            t = t1;
            if (t < 0)
                return false;
        }

        Vec3 pHit = new_origin + t * new_direction;

        if (pHit.y < z_min || pHit.y > z_max)
        {
            if (t == t1)
                return false;
            t = t1;

            pHit = new_origin + t * new_direction;

            if (pHit.y < z_min || pHit.y > z_max)
                return false;
        }

        hit.t = t;
        hit.hitPos = pHit.x * x + pHit.y * y + pHit.z * z + zMin;
        Vec3 normal = normalize(Vec3(pHit.x, 0, pHit.z));
        Vec3 normal2 = normalize(normal.x * x + normal.y * y + normal.z * z);
        hit.hitNormal = dot(new_direction, normal) < 0 ? normal2 : -1 * normal2;
        hit.hitShape = this;
        hit.u = 1; //要実装
        hit.v = 1;

        return true;
    }
};

class Disk : public Figure
{
  public:
    Vec3 x, y, z;             //正規直交基底
    Vec3 x_inv, y_inv, z_inv; //逆行列
    Vec3 center;
    double height, radius, innerRadius;

    Disk(const Vec3 &z, const Vec3 &center, double height, double radius, double innerRadius, std::shared_ptr<Texture> texture, int material) : z(z), center(center), height(height), radius(radius), innerRadius(innerRadius), Figure(texture, material)
    {
        coordinate_system(z, x, y);
        x_inv = x;
        y_inv = y;
        z_inv = z;
        Transpose3x3(x_inv, y_inv, z_inv);
    };
    virtual bool intersect(const Ray &ray, Hit &hit) const
    {
        Vec3 new_direction = (ray.direction.x) * x_inv + (ray.direction.y) * y_inv + (ray.direction.z) * z_inv;
        Vec3 new_origin = ((ray.origin - center).x) * x_inv + ((ray.origin - center).y) * y_inv + ((ray.origin - center).z) * z_inv;

        if (ray.direction.z == 0)
        {
            return false;
        }

        double t = (height - new_origin.z) / new_direction.z;
        if (t <= 0)
        {
            return false;
        }

        Vec3 pHit = new_origin + t * new_direction;
        double dist2 = pHit.x * pHit.x + pHit.y * pHit.y;
        if (dist2 > radius * radius || dist2 < innerRadius * innerRadius)
        {
            return false;
        }
        hit.t = t;
        hit.hitPos = (pHit.x) * x + (pHit.y) * y + (pHit.z) * z + center;
        hit.hitNormal = dot(ray.direction, z) < 0 ? z : -1 * (z);
        hit.hitShape = this;
        hit.u = 1;
        hit.v = 1;

        return true;
    }
};

#endif
