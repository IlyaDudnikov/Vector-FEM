#include "Main.h"
#include "Generator.h"
#include "Grid.h"
#include "VectorFEM.h"
#include "Point.h"
#include <iostream>
#include <iomanip>

int main()
{
    VectorFEM vectorFEM;
    {
        Generator generator;
        vectorFEM.init(generator);
    }
    
    vectorFEM.makeGlobalA();
    vectorFEM.makeGlobalB();

    vectorFEM.setS1_2();

    vectorFEM.solve();
    vectorFEM.outputResultsToFIle();
    vectorFEM.outputAllPointsAndVectors();

    Point resultPoint = vectorFEM.calculatePoint(1.6, 2.6);
    std::cout << std::scientific << std::setprecision(5) << resultPoint.x << '\t' << abs(6.76 - resultPoint.x) << '\n';
    std::cout << std::scientific << std::setprecision(5) << resultPoint.y << '\t' << abs(2.56 - resultPoint.y) << '\n';

    std::cout << "test runner" << '\n';

    //vectorFEM.info();

    return 0;
}