#include "../headers/Particle.h"

#pragma once

class Point
{
    public:
        Point();
        Point(long double x, long double y);
        void resetGVector();
        long double getX();
        long double getY();
        long double getXGV();
        long double getYGV();
        long double getGVM();
        void addGravity(Particle p);
    private:
        long double x_pos;
        long double y_pos;
        long double GXVector = 0;
        long double GYVector = 0;
        long double GVectorMagnitude = 0;
};