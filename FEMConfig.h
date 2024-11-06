#pragma once
#include "Point.h"

class FEMConfig
{
public:
    double mu(int type);
    double gamma(int type);
    Point F(double x, double y, int type);
    Point S1(double x, double y, int type);
};