#include "VectorFEM.h"
#include <fstream>
#include <iomanip>
#include "operations.h"
#include "LOS.h"
#include "Psi.h"

void VectorFEM::init(Generator generator)
{
    generator.getValues(igEdge, jgEdge, ig, jg, elements, points, edges, S1, X, Y);
}

vector<vector<double>> VectorFEM::getLocalG(Element element)
{
    double hx = points[element.node2].x - points[element.node1].x;
    double hy = points[element.node3].y - points[element.node1].y;

    double yDivX = hy / hx;
    double xDivY = hx / hy;

    vector<vector<double>> G_loc(4);
    G_loc[0] = { yDivX, -yDivX, -1, 1 };
    G_loc[1] = { -yDivX ,yDivX, 1, -1 };
    G_loc[2] = { -1, 1, xDivY, -xDivY };
    G_loc[3] = { 1, -1, -xDivY, xDivY };

    G_loc = G_loc * (1 / config.mu(element.areaType));

    return G_loc;
}

vector<vector<double>> VectorFEM::getC(Element element)
{
    vector<vector<double>> C(4);
    C[0] = { 2, 1, 0, 0 };
    C[1] = { 1, 2, 0, 0 };
    C[2] = { 0, 0, 2, 1 };
    C[3] = { 0, 0, 1, 2 };

    double hx = points[element.node2].x - points[element.node1].x;
    double hy = points[element.node3].y - points[element.node1].y;

    C = C * (hx * hy / 6);

    return C;
}

vector<vector<double>> VectorFEM::getLocalM(Element element)
{
    vector<vector<double>> M_loc = getC(element) * config.gamma(element.areaType);
    return M_loc;
}

void VectorFEM::addLocalMatrix(vector<vector<double>> G_loc, vector<vector<double>> M_loc, Element element)
{
    // заносим диагональ
    for (size_t i = 0; i < 4; i++)
    {
        di[element.globalNumberOfUnknowns[i]] += G_loc[i][i] + M_loc[i][i];
    }

    // заносим нижний треугольник
    int ibeg, iend;
    int ind1, ind2; // глоабльные номера рёбер
    int ind; // переменная для цикла дихотомии
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            ind1 = element.globalNumberOfUnknowns[i];
            ind2 = element.globalNumberOfUnknowns[j];

            if (ind1 < ind2)
            {
                int k = ind1;
                ind1 = ind2;
                ind2 = k;
            }

            ibeg = ig[ind1];
            iend = ig[ind1 + 1] - 1;

            while (jg[ibeg] != ind2)
            {
                ind = (ibeg + iend) / 2;
                if (jg[ind] < ind2)
                    ibeg = ind + 1;
                else if (jg[ind] == ind2)
                    ibeg = ind;
                else
                    iend = ind;
            }

            gg[ibeg] += G_loc[i][j] + M_loc[i][j];
        }
    }
}

void VectorFEM::makeGlobalA()
{
    int edgeCount = ig.size() - 1;
    di.resize(edgeCount);
    gg.resize(jg.size());
    vector<vector<double>> G_loc, M_loc;
    for (const Element& element : elements)
    {
        G_loc = getLocalG(element);
        M_loc = getLocalM(element);

        addLocalMatrix(G_loc, M_loc, element);
    }
}

vector<double> VectorFEM::getLocalB(Element element)
{
    double xBeg = points[element.node1].x;
    double xMiddle = points[element.node2].x - points[element.node1].x;
    double xEnd = points[element.node2].x;

    double yBeg = points[element.node1].y;
    double yMiddle = points[element.node3].y - points[element.node1].y;
    double yEnd = points[element.node3].y;

    vector<double> fValues(4);
    fValues[0] = config.F(xBeg, yMiddle, element.areaType).y;
    fValues[1] = config.F(xEnd, yMiddle, element.areaType).y;
    fValues[2] = config.F(xMiddle, yBeg, element.areaType).x;
    fValues[3] = config.F(xMiddle, yEnd, element.areaType).x;

    vector<double> b_loc = getC(element) * fValues;

    return b_loc;
}

void VectorFEM::addLocalB(vector<double> b_loc, Element element)
{
    for (size_t i = 0; i < 4; i++)
    {
        b[element.globalNumberOfUnknowns[i]] += b_loc[i];
    }
}

void VectorFEM::makeGlobalB()
{
    int edgeCount = ig.size() - 1;
    b.resize(edgeCount);
    vector<double> b_loc;
    for (const Element& element : elements)
    {
        b_loc = getLocalB(element);
        addLocalB(b_loc, element);
    }
}

void VectorFEM::setS1()
{
    bool xAxis; // по оси x или y
    for (size_t condNum = 0; condNum < S1.size(); condNum++)
    {
        if ((edges[S1[condNum].edges[0]].node2 - edges[S1[condNum].edges[0]].node1) == 1)
            xAxis = true; // если по оси x
        else
            xAxis = false; // если по оси y

        for (const int& edge : S1[condNum].edges)
        {
            di[edge] = 1e+30;
            if (xAxis)
                b[edge] = config.S1(1, points[edges[edge].node1].y, S1[condNum].formulaType).x * 1e+30;
            else
                b[edge] = config.S1(points[edges[edge].node1].x, 1, S1[condNum].formulaType).y * 1e+30;
        }
    }
}

void VectorFEM::setS1Edge(int edgeNum, int formulaType, bool xAxis)
{
    di[edgeNum] = 1;
    if (xAxis)
    {
        double xMiddle = (points[edges[edgeNum].node2].x + points[edges[edgeNum].node1].x) / 2;
        b[edgeNum] = config.S1(xMiddle, points[edges[edgeNum].node1].y, formulaType).x;
    }
    else
    {
        double yMiddle = (points[edges[edgeNum].node2].y + points[edges[edgeNum].node1].y) / 2;
        b[edgeNum] = config.S1(points[edges[edgeNum].node1].x, yMiddle, formulaType).y;
    }
    for (size_t i = ig[edgeNum]; i < ig[edgeNum + 1]; i++)
    {
        b[jg[i]] -= gg[i] * b[edgeNum];
        gg[i] = 0;
    }
    for (size_t i = edgeNum + 1; i < edges.size(); i++)
    {
        for (size_t j = ig[i]; j < ig[i + 1]; j++)
        {
            if (jg[j] == edgeNum)
            {
                b[i] -= gg[j] * b[edgeNum];
                gg[j] = 0;
                break;
            }
        }
    }
}

void VectorFEM::setS1_2()
{
    bool xAxis; // по оси x или y
    for (size_t condNum = 0; condNum < S1.size(); condNum++)
    {
        if ((edges[S1[condNum].edges[0]].node2 - edges[S1[condNum].edges[0]].node1) == 1)
            xAxis = true; // если по оси x
        else
            xAxis = false; // если по оси y

        for (const int& edge : S1[condNum].edges)
        {
            setS1Edge(edge, S1[condNum].formulaType, xAxis);
        }
    }
}


void VectorFEM::solve()
{
    q.assign(di.size(), 1);
    LOS_sq(di.size(), 1000, 1e-12, q, ig, jg, di, gg, gg, b);
}

void VectorFEM::outputResultsToFIle(std::string fileName)
{
    std::ofstream resultFile(fileName);
    resultFile << "Edges count:," << edges.size() << '\n';
    resultFile << "Element count:," << elements.size() << '\n';
    resultFile << "№,q,q*,|q - q*|\n";
    double qTrue = 0;
    for (size_t i = 0; i < q.size(); i++)
    {
        if ((edges[i].node2 - edges[i].node1) == 1)
        {
            double xMiddle = (points[edges[i].node2].x + points[edges[i].node1].x) / 2;
            qTrue = config.S1(xMiddle, points[edges[i].node1].y, S1[0].formulaType).x; // если горизонтальное ребро
        }
        else
        {
            double yMiddle = (points[edges[i].node2].y + points[edges[i].node1].y) / 2;
            qTrue = config.S1(points[edges[i].node1].x, yMiddle, S1[0].formulaType).y; // если вертикальное ребро
        }
            
        resultFile << i << ',' << std::scientific << std::setprecision(5) << q[i] << ',' << qTrue << ',' << abs(q[i] - qTrue) << '\n';
    }
    resultFile.close();
}

void VectorFEM::outputAllPointsAndVectors(std::string pointsFileName, std::string vectorsFileName)
{
    std::ofstream pointsFile(pointsFileName), vectorsFile(vectorsFileName);
    double hx = 0, hy = 0;
    double xCurrent, yCurrent;
    int numberXIntervals = 3, numberYIntervals = 3;
    for (size_t yNum = 0; yNum < Y.size() - 1; yNum++)
    {
        hy = (Y[yNum + 1] - Y[yNum]) / numberYIntervals;
        for (size_t xNum = 0; xNum < X.size() - 1; xNum++)
        {
            hx = (X[xNum + 1] - X[xNum]) / numberXIntervals;

            for (size_t yStep = 0; yStep < numberYIntervals + 1; yStep++)
            {
                for (size_t xStep = 0; xStep < numberXIntervals + 1; xStep++)
                {
                    xCurrent = X[xNum] + xStep * hx;
                    yCurrent = Y[yNum] + yStep * hy;

                    Point pointValue = calculatePoint(xCurrent, yCurrent);

                    pointsFile << xCurrent << ' ' << yCurrent << '\n';
                    vectorsFile << pointValue.x / 10 << ' ' << pointValue.y / 10 << '\n';
                }
            }
        }
    }
    pointsFile.close();
    vectorsFile.close();
}

Point VectorFEM::calculatePoint(double x, double y)
{
    int elemNum;
    for (elemNum = 0; elemNum < elements.size(); elemNum++)
    {
        if ((x >= points[elements[elemNum].node1].x) && (x <= points[elements[elemNum].node2].x) && // условие попадания в элемент
            (y >= points[elements[elemNum].node1].y) && (y <= points[elements[elemNum].node3].y))
            break;
    }

    double x1 = points[elements[elemNum].node1].x;
    double x2 = points[elements[elemNum].node2].x;
    double y1 = points[elements[elemNum].node1].y;
    double y2 = points[elements[elemNum].node3].y;

    Point result(0, 0);
    result.y += q[elements[elemNum].globalNumberOfUnknowns[0]] * psi1(x, x1, x2).y;
    result.y += q[elements[elemNum].globalNumberOfUnknowns[1]] * psi2(x, x1, x2).y;
    result.x += q[elements[elemNum].globalNumberOfUnknowns[2]] * psi3(y, y1, y2).x;
    result.x += q[elements[elemNum].globalNumberOfUnknowns[3]] * psi4(y, y1, y2).x;

    return result;
}

void VectorFEM::info()
{
    std::cout << "Edges count: " << edges.size() << '\n';
    std::cout << "Element count: " << elements.size() << '\n';
}