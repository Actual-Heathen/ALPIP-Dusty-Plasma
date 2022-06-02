#pragma once

class Particle
{
    struct Vector
    {
        float xVec = 0;
        float yVec = 0;
        float mag = 0;
    };

    public:
        Particle();
        Particle(float x, float y, float m);
        float getX();
        float getY();
        float getMass();
        float getSpeed();
        void addAcceleration(float x, float y);
        void move(float s);
        void resetAcc();
        void resetVel();
    private:
        void setX(float x);
        void setY(float y);
        void setMass(float m);

        float x_pos;
        float y_pos;
        float mass;
        float speed = 0;
        float velX = 0;
        float velY = 0;
        float accel = 0;
        float accelX = 0;
        float accelY = 0;
};