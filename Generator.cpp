#include "Generator.h"
#include <fstream>
#include <cassert>
#include <sstream>

void Generator::makeSubareas(std::string areaFile)
{
    std::ifstream area(areaFile);
    int NxW, NyW, L;

    // заполнение массива XW (границы подобластей по x)
    area >> NxW;
    XW.resize(NxW);
    for (size_t i = 0; i < NxW; i++)
    {
        area >> XW[i];
    }

    // заполнение массива YW (границы подобластей по y)
    area >> NyW;
    YW.resize(NyW);
    for (size_t i = 0; i < NyW; i++)
    {
        area >> YW[i];
    }

    // заполнение подобластей
    area >> L; // количество подобластей
    subareas.resize(L);
    for (size_t i = 0; i < L; i++)
    {
        area >> subareas[i].type;

        area >> subareas[i].xStart;
        area >> subareas[i].xEnd;
        area >> subareas[i].yStart;
        area >> subareas[i].yEnd;
    }
    area.close();
}

void Generator::subdivideGrid(std::string gridFile)
{
    std::ifstream grid(gridFile);

    // формирование массива X (всех иксов сетки)
    vector<int> numIntervals(XW.size() - 1);
    vector<double> dischargeCoefficients(XW.size() - 1);
    int totalNumIntervals = 0;
    for (size_t i = 0; i < XW.size() - 1; i++)
    {
        grid >> numIntervals[i];
        grid >> dischargeCoefficients[i];
        totalNumIntervals += numIntervals[i];
    }

    X.resize(totalNumIntervals + 1);
    XI.resize(XW.size());
    X[0] = XW[0];
    XI[0] = 0;
    for (size_t i = 0, currentXPosition = 0; i < XW.size() - 1; i++)
    {
        currentXPosition++;
        if (dischargeCoefficients[i] == 1)
        {
            double hx = (XW[i + 1] - XW[i]) / numIntervals[i];
            for (size_t j = 1; j < numIntervals[i]; j++)
            {
                X[currentXPosition] = XW[i] + j * hx;
                currentXPosition++;
            }
        }
        else
        {
            double hx = (XW[i + 1] - XW[i]) * (dischargeCoefficients[i] - 1)
                / (pow(dischargeCoefficients[i], numIntervals[i]) - 1); // первый шаг
            for (size_t j = 0; j < numIntervals[i] - 1; j++)
            {
                X[currentXPosition] = X[currentXPosition - 1] + hx * pow(dischargeCoefficients[i], j);
                currentXPosition++;
            }
        }
        X[currentXPosition] = XW[i + 1];
        XI[i + 1] = currentXPosition;
    }

    // формирование массива Y (всех игриков сетки)
    numIntervals.resize(YW.size() - 1);
    dischargeCoefficients.resize(YW.size() - 1);
    totalNumIntervals = 0;
    for (size_t i = 0; i < YW.size() - 1; i++)
    {
        grid >> numIntervals[i];
        grid >> dischargeCoefficients[i];
        totalNumIntervals += numIntervals[i];
    }
    grid.close();

    Y.resize(totalNumIntervals + 1);
    YI.resize(YW.size());
    Y[0] = YW[0];
    YI[0] = 0;
    for (size_t i = 0, currentYPosition = 0; i < YW.size() - 1; i++)
    {
        currentYPosition++;
        if (dischargeCoefficients[i] == 1)
        {
            double hy = (YW[i + 1] - YW[i]) / numIntervals[i];
            for (size_t j = 1; j < numIntervals[i]; j++)
            {
                Y[currentYPosition] = YW[i] + j * hy;
                currentYPosition++;
            }
        }
        else
        {
            double hy = (YW[i + 1] - YW[i]) * (dischargeCoefficients[i] - 1)
                / (pow(dischargeCoefficients[i], numIntervals[i]) - 1); // первый шаг
            for (size_t j = 0; j < numIntervals[i] - 1; j++)
            {
                Y[currentYPosition] = Y[currentYPosition - 1] + hy * pow(dischargeCoefficients[i], j);
                currentYPosition++;
            }
        }
        Y[currentYPosition] = YW[i + 1];
        YI[i + 1] = currentYPosition;
    }
}

int Generator::getAreaType(int xNum, int yNum) // не может определить тип узлов, лежащих на правой и верхней границе расчётной области
{
    for (const Subarea& subarea : subareas)
    {
        if (xNum >= XI[subarea.xStart] && xNum < XI[subarea.xEnd] &&
            yNum >= YI[subarea.yStart] && yNum < YI[subarea.yEnd])
            return subarea.type;
    }
    assert(false); // не вошли ни в одну из областей
}

int Generator::globalNumberOfNodeOfElement(int xNum, int yNum, int localNodeNum) // xSize - число узлов по x
{
    assert(localNodeNum >= 0 && localNodeNum <= 3);

    if ((localNodeNum == 0) || (localNodeNum == 1))
        return X.size() * yNum + xNum + localNodeNum;

    // если loсalNodeNum 2 или 3
    return X.size() * (yNum + 1) + xNum + localNodeNum - 2;
}

int Generator::edgeNumberByNodes(int node1, int node2)
{
    if (node1 > node2)
    {
        int k = node1;
        node1 = node2;
        node2 = k;
    }

    for (size_t i = 0; i < edges.size(); i++)
    {
        if ((edges[i].node1 == node1) && (edges[i].node2 == node2))
            return i;
    }

    assert(false);
    return -1;
}

void Generator::makeElements()
{
    int numElements = (X.size() - 1) * (Y.size() - 1);
    elements.resize(numElements);
    int elementNum = 0;
    for (size_t yNum = 0; yNum < Y.size() - 1; yNum++)
    {
        for (size_t xNum = 0; xNum < X.size() - 1; xNum++)
        {
            elements[elementNum].node1 = globalNumberOfNodeOfElement(xNum, yNum, 0);
            elements[elementNum].node2 = globalNumberOfNodeOfElement(xNum, yNum, 1);
            elements[elementNum].node3 = globalNumberOfNodeOfElement(xNum, yNum, 2);
            elements[elementNum].node4 = globalNumberOfNodeOfElement(xNum, yNum, 3);

            elements[elementNum].areaType = getAreaType(xNum, yNum);

            
            /*elements[elementNum].globalNumberOfUnknowns[0] = elements[elementNum].node1;
            elements[elementNum].globalNumberOfUnknowns[1] = elements[elementNum].node2;
            elements[elementNum].globalNumberOfUnknowns[2] = elements[elementNum].node3;
            elements[elementNum].globalNumberOfUnknowns[3] = elements[elementNum].node4;*/

            elements[elementNum].globalNumberOfUnknowns.resize(4);
            elements[elementNum].globalNumberOfUnknowns[0] = (2 * X.size() - 1) * yNum + X.size() - 1 + xNum;
            elements[elementNum].globalNumberOfUnknowns[1] = elements[elementNum].globalNumberOfUnknowns[0] + 1;
            elements[elementNum].globalNumberOfUnknowns[2] = (2 * X.size() - 1) * yNum + xNum;
            elements[elementNum].globalNumberOfUnknowns[3] = elements[elementNum].globalNumberOfUnknowns[2] + (2 * X.size() - 1);

            elementNum++;
        }
    }
}

void Generator::makePoints()
{
    points.resize(X.size() * Y.size());
    int pointNum = 0;
    for (size_t yNum = 0; yNum < Y.size(); yNum++)
    {
        for (size_t xNum = 0; xNum < X.size(); xNum++)
        {
            points[pointNum].x = X[xNum];
            points[pointNum].y = Y[yNum];

            pointNum++;
        }
    }
}

void Generator::makeEdges()
{
    int edgeCount = elements.size() * 2 + X.size() - 1 + Y.size() - 1;
    edges.resize(edgeCount);
    int edgeNum = 0; // номер текущего ребра
    int elementNum; // номер текущего элемента
    for (size_t yNum = 0; yNum < Y.size() - 1; yNum++)
    {
        // рёбра по горизонтали
        elementNum = yNum * (X.size() - 1);
        for (size_t xNum = 0; xNum < X.size() - 1; xNum++)
        {
            edges[edgeNum].node1 = elements[elementNum].node1;
            edges[edgeNum].node2 = elements[elementNum].node2;

            edgeNum++;
            elementNum++;
        }
        
        // рёбра по вертикали
        elementNum = yNum * (X.size() - 1);
        for (size_t xNum = 0; xNum < X.size() - 1; xNum++)
        {
            edges[edgeNum].node1 = elements[elementNum].node1;
            edges[edgeNum].node2 = elements[elementNum].node3;

            edgeNum++;
            elementNum++;
        }

        // последнее ребро по вертикали
        elementNum--;
        edges[edgeNum].node1 = elements[elementNum].node2;
        edges[edgeNum].node2 = elements[elementNum].node4;
        edgeNum++;
    }

    // рёбра по горизонтали самой верхней границы
    elementNum = (Y.size() - 2) * (X.size() - 1);
    for (size_t xNum = 0; xNum < X.size() - 1; xNum++)
    {
        edges[edgeNum].node1 = elements[elementNum].node3;
        edges[edgeNum].node2 = elements[elementNum].node4;

        edgeNum++;
        elementNum++;
    }
}

void Generator::readBoundaryConditions(std::string boundaryConditionsFile)
{
    std::ifstream conditions(boundaryConditionsFile);
    std::string line;
    int conditionType, formula;
    int xBeg, xEnd, yBeg, yEnd;
    int node1, node2;
    while (std::getline(conditions, line)) {
        std::istringstream iss(line);
        int number;
        iss >> conditionType;
        assert(conditionType >= 1 && conditionType <= 3);
        BoundaryCondition boudaryCondition;
        iss >> boudaryCondition.formulaType;
        
        iss >> xBeg;
        iss >> xEnd;
        iss >> yBeg;
        iss >> yEnd;

        int xStartNum, yStartNum;
        node1 = YI[yBeg] * X.size() + XI[xBeg];
        if (xBeg == xEnd)
        {
            xStartNum = XI[xBeg];
            yStartNum = YI[yBeg] + 1;
        }
        else
        {
            xStartNum = XI[xBeg] + 1;
            yStartNum = YI[yBeg];
        }

        int edgeNum;
        for (size_t yNum = yStartNum; yNum <= YI[yEnd]; yNum++)
        {
            for (size_t xNum = xStartNum; xNum <= XI[xEnd]; xNum++)
            {
                node2 = yNum * X.size() + xNum;
                edgeNum = edgeNumberByNodes(node2, node1);
                boudaryCondition.edges.push_back(edgeNum);
                node1 = node2;
            }
        }

        S1.push_back(boudaryCondition);

        /*boudaryCondition.nodes.resize(4);
        for (size_t i = 0; i < 4; i++)
        {
            iss >> boudaryCondition.nodes[i];
        }*/

        //if (conditionType == 1)
            //S1.push_back(boudaryCondition);
        /*else if (conditionType == 2)
            S2.push_back(boudaryCondition);
        else
            S3.push_back(boudaryCondition);*/
    }
    conditions.close();
}

void Generator::makePortrait()
{
    int edgeCount = elements.size() * 2 + X.size() - 1 + Y.size() - 1;
    vector<vector<int>> list(edgeCount);
    int ind1, ind2;
    for (size_t elementNum = 0; elementNum < elements.size(); elementNum++) // цикл по элементам
    {
        for (size_t i = 0; i < 4; i++)
        {
            for (size_t j = i + 1; j < 4; j++)
            {
                ind1 = elements[elementNum].globalNumberOfUnknowns[i];
                ind2 = elements[elementNum].globalNumberOfUnknowns[j];

                if (ind1 > ind2)
                {
                    int val = ind2;
                    ind2 = ind1;
                    ind1 = val;
                }

                // ind2 > ind1
                assert(ind2 > ind1);

                bool added = false;
                for (size_t k = 0; k < list[ind2].size(); k++)
                {
                    if (ind1 < list[ind2][k])
                    {
                        list[ind2].insert(list[ind2].begin() + k, ind1);
                        added = true;
                        break;
                    }
                    if (ind1 == list[ind2][k])
                    {
                        added = true;
                        break;
                    }
                }
                if (!added)
                    list[ind2].push_back(ind1);
            }
        }
    }

    // создание портрета по списку
    ig.resize(edgeCount + 1);
    ig[0] = 0;
    for (size_t i = 0; i < edgeCount; i++)
    {
        ig[i + 1] = ig[i];
        for (size_t j = 0; j < list[i].size(); j++)
        {
            jg.push_back(list[i][j]);
            ig[i + 1]++;
        }
    }
}

void Generator::getValues(vector<int>& igEdge, vector<int>& jgEdge, vector<int>& ig, vector<int>& jg,
    vector<Element>& elements, vector<Point>& points, vector<Edge>& edges, vector<BoundaryCondition>& S1, vector<double>& X, vector<double>& Y)
{
    igEdge = this->igEdge;
    jgEdge = this->jgEdge;
    ig = this->ig;
    jg = this->jg;
    elements = this->elements;
    points = this->points;
    edges = this->edges;
    S1 = this->S1;
    X = this->X;
    Y = this->Y;
}

Generator::Generator(std::string areaFile, std::string gridFile, std::string boundaryConditionsFile)
{
    makeSubareas(areaFile);
    subdivideGrid(gridFile);
    makeElements();
    makePoints();
    makeEdges();
    makePortrait();
    readBoundaryConditions(boundaryConditionsFile);
}