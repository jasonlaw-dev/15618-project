#ifndef GRID_ITERATOR_H
#define GRID_ITERATOR_H

#include "partitioninfo.h"
#include <utility>

class GridIterator {
private:
    PartitionInfo partition;
    int borderWidth;
    int i;
    int j;
    bool innerFinished;
    int innerTopPixel;
    int innerBottomPixel;
    int innerLeftPixel;
    int innerRightPixel;
    int PROCID;
    bool communicate;
public:
    GridIterator(PartitionInfo partition, int borderWidth, int PROCID);
    std::pair<int, int> next();
    bool shouldCommunicate();
};

#endif