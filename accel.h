#ifndef ACCEL_H
#define ACCEL_H
#include <vector>
#include <memory>
#include "figure.h"
#include "hit.h"
#include "ray.h"

class Accel
{
  public:
    std::vector<std::shared_ptr<Figure>> shapes;
    Accel(){};

    void add(std::shared_ptr<Figure> p)
    {
        shapes.push_back(p);
    };
    bool intersect(const Ray &ray, Hit &hit) const
    {
        bool isHit = false;

        hit.t = 1000000;

        Hit hit_each;
        for (auto p : shapes)
        {
            if (p->intersect(ray, hit_each))
            {
                isHit = true;

                if (hit_each.t < hit.t)
                {
                    hit = hit_each;
                }
            }
        }
        return isHit;
    };
};
#endif
