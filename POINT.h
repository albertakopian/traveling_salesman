class point {
public:
    point(){
        x = 0;
        y = 0;
    }
    point(double x_, double y_){
        x = x_;
        y = y_;
    }
    double x;
    double y;
};
double distance(const point &p1, const point& p2)
{
    return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}