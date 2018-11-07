#ifndef MYRANDOM_H
#define MYRANDOM_H

#include <random>
#include <cmath>
#include "vec3.h"

std::random_device dev_rnd;
std::mt19937 mt(dev_rnd());
std::uniform_real_distribution<> dist1(0, 1);

inline double rnd()
{
    return dist1(mt);
}

Vec3 random_unit_circle(const Vec3 &n)
{
    double u = rnd();
    double r = rnd();

    double y = 0;
    double x = r * std::cos(2 * M_PI * u);
    double z = r * std::sin(2 * M_PI * u);

    Vec3 xv, zv;
    orthonormalBasis(n, xv, zv);

    return x * xv + y * n + z * zv;
}

Vec3 randomHemisphere(const Vec3 &n)
{
    double u = rnd();
    double v = rnd();

    double y = u;
    double x = std::sqrt(1 - u * u) * std::cos(2 * M_PI * v);
    double z = std::sqrt(1 - u * u) * std::sin(2 * M_PI * v);

    Vec3 xv, zv;
    orthonormalBasis(n, xv, zv);

    return x * xv + y * n + z * zv; //{xv,n,v}を基底とした時の[x,y,z]
}

Vec3 randomCosineHemisphere(double &pdf, const Vec3 &n)
{
    double u = rnd();
    double v = rnd();

    double theta = 0.5 * std::acos(1 - 2 * u);
    double phi = 2 * M_PI * v;
    pdf = 1 / M_PI * std::cos(theta);

    double x = std::cos(phi) * std::sin(theta);
    double y = std::cos(theta);
    double z = std::sin(phi) * std::sin(theta);
    Vec3 xv, zv;
    orthonormalBasis(n, xv, zv);
    return x * xv + y * n + z * zv;
}

#endif