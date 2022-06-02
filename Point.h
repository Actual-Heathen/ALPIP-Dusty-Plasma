class Point
{
    public:
        Point(float x, float y);
        void resetGVector();
        float getX();
        float getY();
        float getXGV();
        float getYGV();
        float getGVM();
    private:
        float x_pos;
        float y_pos;
        float GXVector;
        float GYVector;
        float GVectorMagnitude;
};