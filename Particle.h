class Particle
{
    public:
        Particle();
        Particle(float x, float y, float m);
        float getX();
        float getY();
        float getMass();
    private:
        void setX(float x);
        void setY(float y);
        void setMass(float m);

        float x_pos;
        float y_pos;
        float mass;
};