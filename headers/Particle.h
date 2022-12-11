#pragma once
#include <vector>

class Particle
{
    public:
        Particle();
        Particle( double x,  double y,  double m, double vx, double vy);
         double getX();
         double getY();
         double getMass();
         double getSpeed();
        void addAcceleration(double spacing, std::vector<std::vector<double>> dpsix, std::vector<std::vector<double>> dpsiy,std::vector<std::vector<double>> dphix,std::vector<std::vector<double>> dphiy, double E[3], double B[3], double time);
        void move(double s, double size);
        void resetAcc();
        void resetVel();
        void setX( double x);
        void setY( double y);
        void setMass( double m);
    private:

        double x_pos;
        double y_pos;
        double mass;
        double speed = 0;
        double velX = 0;
        double velY = 0;
        double velZ = 0;
        double accel = 0;
        double accelX = 0;
        double accelY = 0;
        double accelZ = 0;
        double charge;
};
