#pragma once
#include <vector>
#include <string>
#include "Element.h"
#include "Point.h"
#include "Edge.h"
#include "BoundaryCondition.h"
#include "Subarea.h"

using std::vector;

class Generator
{
private:
    vector<int> igEdge, jgEdge;
    vector<int> ig, jg;
    vector<Element> elements;
    vector<Point> points;
    vector<Edge> edges;
    vector<BoundaryCondition> S1;
    vector<double> XW, YW, X, Y;
    vector<int> XI, YI;
    vector<Subarea> subareas;

    void makeSubareas(std::string areaFile = "area.txt");
    void subdivideGrid(std::string gridFile = "grid.txt");
    int getAreaType(int xNum, int yNum);
    int globalNumberOfNodeOfElement(int xNum, int yNum, int localNodeNum);
    int edgeNumberByNodes(int node1, int node2);
    void makeElements();
    void makePoints();
    void makeEdges();
    void readBoundaryConditions(std::string boundaryConditionsFile = "boundaryConditions.txt");
    void makePortrait(); // ig, jg


public:
    Generator(std::string areaFile = "area.txt", std::string gridFile = "grid.txt", std::string boundaryConditionsFile = "boundaryConditions.txt");
    void getValues(vector<int>& igEdge, vector<int>& jgEdge, vector<int>& ig, vector<int>& jg,
        vector<Element>& elements, vector<Point>& points, vector<Edge>& edges, vector<BoundaryCondition>& S1, vector<double>& X, vector<double>& Y);
};