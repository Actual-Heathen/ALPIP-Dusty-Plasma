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
        void addAcceleration(double spacing, std::vector<std::vector<double>> dpsix, std::vector<std::vector<double>> dpsiy, double E[3], double B[3], double time);
        void move(double s, double size);
        void resetAcc();
        void resetVel();
        void setX(long double x);
        void setY(long double y);
        void setMass(long double m);
    private:

        double x_pos;
        double y_pos;
        double mass;
        double speed = 0;
        double velX = 0;
        double velY = 0;
        double accel = 0;
        double accelX = 0;
        double accelY = 0;
        double charge;
};
