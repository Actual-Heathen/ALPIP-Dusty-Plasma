#pragma once

class Particle
{
    struct Vector
    {
        long double xVec = 0;
        long double yVec = 0;
        long double mag = 0;
    };

    public:
        Particle();
        Particle(long double x, long double y, long double m);
        long double getX();
        long double getY();
        long double getMass();
        long double getSpeed();
        void addAcceleration(double rho[][3], double spacing);
        void move(long double s);
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