#pragma once
#include <vector>
#include "Edge.h"

class BoundaryCondition
{
public:
    int formulaType; // для кейса в функции
    std::vector<int> edges; // рёбра, в которых задано это условие
};