#ifndef IMAGE_H
#define IMAGE_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "vec3.h"

class Image
{
  public:
    int width;
    int height;
    Vec3 *data;

    Image(int _width, int _height)
    {
        width = _width;
        height = _height;
        data = new Vec3[width * height];
    };

    ~Image()
    {
        delete[] data;
    };

    Vec3 getPixel(int i, int j) const
    {
        return data[width * i + j];
    };

    void setPixel(int i, int j, const Vec3 &color_imfo) const
    {
        data[width * i + j] = color_imfo;
    };

    void ppm_output() const
    {
        std::ofstream file("output.ppm");
        file << "P3" << std::endl;
        file << width << " " << height << std::endl;
        file << "255" << std::endl;
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                Vec3 color_screen = 255 * this->getPixel(j, i);
                int r = (int)color_screen.x;
                int g = (int)color_screen.y;
                int b = (int)color_screen.z;
                file << r << " " << g << " " << b << std::endl;
            }
        }
        file.close();
    };

    void gamma_correction_clamp()
    {
        for (int i = 0; i < width; i++)
        {
            for (int j = 0; j < height; j++)
            {
                Vec3 col = this->getPixel(i, j);
                col.x = col.x <= 0 ? 0 : std::min(std::pow(col.x, 1.0 / 2.2), 1.0);
                col.y = col.y <= 0 ? 0 : std::min(std::pow(col.y, 1.0 / 2.2), 1.0);
                col.z = col.z <= 0 ? 0 : std::min(std::pow(col.z, 1.0 / 2.2), 1.0);
                this->setPixel(i, j, col);
            }
        }
    };
};

#endif
