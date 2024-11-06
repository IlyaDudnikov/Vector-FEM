#pragma once
#include "Generator.h"
#include "FEMConfig.h"

class VectorFEM
{
private:
    vector<int> igEdge, jgEdge;
    vector<double> X, Y;
    vector<int> ig, jg;
    vector<double> di, gg, b, q;
    vector<Element> elements;
    vector<Point> points;
    vector<Edge> edges;
    vector<BoundaryCondition> S1;
    FEMConfig config;

    vector<vector<double>> getLocalG(Element element);
    vector<vector<double>> getC(Element element);
    vector<vector<double>> getLocalM(Element element);
    void addLocalMatrix(vector<vector<double>> G_loc, vector<vector<double>> M_loc, Element element);
    vector<double> getLocalB(Element element);
    void addLocalB(vector<double> b_loc, Element element);
    void setS1Edge(int edgeNum, int formulaType, bool xAxis);

public:
    void init(Generator generator);
    void makeGlobalA();
    void makeGlobalB();
    void setS1();
    void setS1_2();
    void solve();
    void outputResultsToFIle(std::string fileName = "output.csv");
    void outputAllPointsAndVectors(std::string pointsFileName = "points.txt", std::string vectorsFileName = "vectors.txt");
    Point calculatePoint(double x, double y);
    void info();
};
