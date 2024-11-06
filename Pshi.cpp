#include "Psi.h"

Point psi1(double x, double x1, double x2)
{
    Point result(0, 0);
    double hx = x2 - x1;
    result.y = (x2 - x) / hx;
    return result;
}


Point psi2(double x, double x1, double x2)
{
    Point result(0, 0);
    double hx = x2 - x1;
    result.y = (x - x1) / hx;
    return result;
}

Point psi3(double y, double y1, double y2)
{
    Point result(0, 0);
    double hy = y2 - y1;
    result.x = (y2 - y) / hy;
    return result;
}

Point psi4(double y, double y1, double y2)
{
    Point result(0, 0);
    double hy = y2 - y1;
    result.x = (y - y1) / hy;
    return result;
}