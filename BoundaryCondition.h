#pragma once
#include <vector>
#include "Edge.h"

class BoundaryCondition
{
public:
    int formulaType; // ��� ����� � �������
    std::vector<int> edges; // ����, � ������� ������ ��� �������
};