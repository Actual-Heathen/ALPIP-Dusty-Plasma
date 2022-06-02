#include "../headers/Particle.h"

class Point
{
    public:
        Point();
        Point(float x, float y);
        void resetGVector();
        float getX();
        float getY();
        float getXGV();
        float getYGV();
        float getGVM();
        void addGravity(Particle p);
    private:
        float x_pos;
        float y_pos;
        float GXVector = 0;
        float GYVector = 0;
        float GVectorMagnitude = 0;
};