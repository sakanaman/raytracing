#include "vec3.h"
#include "image.h"
#include "camera.h"
#include "ray.h"
#include "figure.h"
#include "accel.h"
#include <cmath>
#include <iostream>
#include <omp.h>
#include "myrandom.h"
#include "texture.h"
#include "ibl.h"
#include "perlin.h"

Accel accel;

Vec3 getColor(const Ray &ray, IBL &ibl, double roulette = 1.0, int depth = 0)
{
    //if(depth > 100) return Vec3(0,0,0);//深度で打ち切り

    if (rnd() > roulette)
    { //ロシアンルーレットで打ち切り
        return Vec3(0, 0, 0);
    }
    if (depth != 0)
        roulette *= 0.96;
    Hit hit;
    if (accel.intersect(ray, hit))
    {

        if (hit.hitShape->material == 0)
        { //difuse
            double pdf;
            Ray nextRay(hit.hitPos + 0.001 * hit.hitNormal, randomCosineHemisphere(pdf, hit.hitNormal));
            double cos_term = std::max(dot(nextRay.direction, hit.hitNormal), 0.0);
            return 1 / roulette * 1 / pdf * hit.hitShape->texture->getColor1(hit.u, hit.v, hit.hitPos) * getColor(nextRay, ibl, roulette, depth + 1) / M_PI * cos_term;
        }
        else if (hit.hitShape->material == 1)
        { //metal
            Ray nextRay(hit.hitPos + 0.001 * hit.hitNormal, reflect(ray.direction, hit.hitNormal));
            if (dot(nextRay.direction, hit.hitNormal) > 0)
            {
                return 1 / roulette * hit.hitShape->texture->getColor1(hit.u, hit.v, hit.hitPos) * getColor(nextRay, ibl, roulette, depth + 1);
            }
            return Vec3(0, 0, 0);
        }
        else if (hit.hitShape->material == 2)
        { //dielectric
            const Vec3 orienting_normal = dot(hit.hitNormal, ray.direction) < 0.0 ? hit.hitNormal : (-1.0 * hit.hitNormal);
            const Ray reflection_ray = Ray(hit.hitPos + 0.001 * orienting_normal, ray.direction - hit.hitNormal * 2.0 * dot(hit.hitNormal, ray.direction));
            const bool into = dot(hit.hitNormal, orienting_normal) > 0.0;

            const double nc = 1.0;  // 真空の屈折率
            const double nt = 1.48; // オブジェクトの屈折率
            const double nnt = into ? nc / nt : nt / nc;
            const double ddn = dot(ray.direction, orienting_normal);
            const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

            if (cos2t < 0.0)
            { // 全反射
                return getColor(reflection_ray, ibl, roulette, depth + 1) * 1 / roulette;
            }
            const Ray refraction_ray = Ray(hit.hitPos - orienting_normal * 0.001,
                                           normalize(ray.direction * nnt - hit.hitNormal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t))));

            const double a = nt - nc, b = nt + nc;
            const double R0 = (a * a) / (b * b);

            const double c = 1.0 - (into ? -ddn : dot(refraction_ray.direction, -1.0 * orienting_normal));
            const double Re = R0 + (1.0 - R0) * pow(c, 5.0);        // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
            const double nnt2 = pow(into ? nc / nt : nt / nc, 2.0); // レイの運ぶ放射輝度は屈折率の異なる物体間を移動するとき、屈折率の比の二乗の分だけ変化する。
            const double Tr = (1.0 - Re) * nnt2;                    // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

            const double prob = 0.25 + 0.5 * Re;
            if (depth > 2)
            {
                if (rnd() < prob)
                {
                    return getColor(reflection_ray, ibl, roulette, depth + 1) * Re * (1 / prob) * (1 / roulette);
                }
                else
                {
                    return getColor(refraction_ray, ibl, roulette, depth + 1) * Tr * (1 / (1 - prob)) * (1 / roulette);
                }
            }
            else
            {
                return getColor(reflection_ray, ibl, roulette, depth + 1) * Re * 1 / roulette +
                       getColor(refraction_ray, ibl, roulette, depth + 1) * Tr * 1 / roulette;
            }
        }
        else
        {
            return Vec3(0, 0, 0);
        }
    }
    else
    {

        return ibl.getColor2(ray);
    }
}

int main()
{

    IBL ibl("texture_ibl/PaperMill_E_3k.hdr");

    Image img(512, 512); //横・縦幅をもつImageのオブジェクトの生成
    Vec3 lookfrom(0, 6 * std::sin(M_PI / 12), -6 * std::cos(M_PI / 12));
    Thin_Lens_Camera cam(lookfrom, normalize(-1 * lookfrom), 1.5, Vec3(0, 0, 0), 6.5);

    accel.add(std::make_shared<Sphere>(Vec3(0, 0, 0), 1.0, std::make_shared<ImageTexture>("texture_ibl/earthmap.jpg"), 1));
    accel.add(std::make_shared<Tri>(Vec3(-5, -1, 5), Vec3(5, -1, -5), Vec3(-5, -1, -5),
                                    std::make_shared<CheckerTexture>(Vec3(0.99, 0.99, 0.99), Vec3(0.5, 0.5, 0.5)), 0)); //下の地面1
    accel.add(std::make_shared<Tri>(Vec3(-5, -1, 5), Vec3(5, -1, -5), Vec3(5, -1, 5),
                                    std::make_shared<CheckerTexture>(Vec3(0.99, 0.99, 0.99), Vec3(0.5, 0.5, 0.5)), 0)); //下の地面2

    int samples = 1000;
#pragma omp parallel for
    for (int k = 0; k < 1000; k++)
    {
        for (int i = 0; i < img.width; i++)
        {
            for (int j = 0; j < img.height; j++)
            {
                double u = (2.0 * (i + rnd()) - img.width) / img.width;
                double v = (2.0 * (j + rnd()) - img.height) / img.height;
                Ray ray = cam.getRay(u, v);
                Vec3 color = getColor(ray, ibl);
                img.setPixel(i, j, img.getPixel(i, j) + 1.0 / static_cast<double>(samples) * color);
            }
        }
    }
    img.gamma_correction_clamp();
    img.ppm_output();
}