#pragma once
#include <vector>
#include "element.h"
#include "point.h"

using std::vector;

class Grid
{
private:
    vector<int> igEdge, jgEdge;
    vector<Element> elements;
    vector<Point> points;


public:
    Grid(vector<int> igEdge, vector<int> jgEdge, vector<Element> elements, vector<Point> points)
        : igEdge(igEdge), jgEdge(jgEdge), elements(elements), points(points) {}
    Grid() {}
    //Grid(std::string areaFile = "area.txt", std::string gridFile = "grid.txt", std::string boundaryConditionsFile = "boundaryConditions.txt");
};