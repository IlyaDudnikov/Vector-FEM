#pragma once
class Subarea
{
public:
    int type;
    int xStart, xEnd;
    int yStart, yEnd;

    Subarea() {}
    Subarea(int type, int xStart, int xEnd, int yStart, int yEnd)
        : type(type), xStart(xStart), xEnd(xEnd), yStart(yStart), yEnd(yEnd) {}
};

