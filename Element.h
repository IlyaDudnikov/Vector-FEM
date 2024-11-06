#pragma once
#include <vector>

class Element
{
public:
    int node1, node2, node3, node4;
    std::vector<int> globalNumberOfUnknowns;
    int areaType;
};