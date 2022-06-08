#pragma once
#include <vector>

class Particle
{
    public:
        Particle();
        Particle(long double x, long double y, long double m, double vx, double vy);
        long double getX();
        long double getY();
        long double getMass();
        long double getSpeed();
        void addAcceleration(double spacing, std::vector<std::vector<double>> rho);
        void move(double s, double size);
        void resetAcc();
        void resetVel();
    private:
        void setX(long double x);
        void setY(long double y);
        void setMass(long double m);

        long double x_pos;
        long double y_pos;
        long double mass;
        long double speed = 0;
        long double velX = 0;
        long double velY = 0;
        long double accel = 0;
        long double accelX = 0;
        long double accelY = 0;
};