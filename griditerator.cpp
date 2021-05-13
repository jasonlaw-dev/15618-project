#include <algorithm>
#include "griditerator.h"
#include <iostream>
#include <stdio.h>

// this initializes the iterator with variables inner region. Inner region refers to the region that does not require communication
// size of inner region is determined by borderWidth, which is the width of the border of region that requires communication
// if border width is 0, it iterates top to bottom, left to right.
GridIterator::GridIterator(PartitionInfo partition, int borderWidth, int PROCID): partition(partition), borderWidth(borderWidth), innerFinished(false), PROCID(PROCID), communicate(false) {
    innerTopPixel = std::min(partition.topPixel + borderWidth, partition.bottomPixel);
    if (partition.neigbhors[1] == -1) {
        innerTopPixel = partition.topPixel;
    }
    innerBottomPixel = std::max(partition.bottomPixel - borderWidth, partition.topPixel);
    if (partition.neigbhors[7] == -1) {
        innerBottomPixel = partition.bottomPixel;
    }
    innerLeftPixel = std::min(partition.leftPixel + borderWidth, partition.rightPixel);
    if (partition.neigbhors[3] == -1) {
        innerLeftPixel = partition.leftPixel;
    }
    innerRightPixel = std::max(partition.rightPixel - borderWidth, partition.leftPixel);
    if (partition.neigbhors[5] == -1) {
        innerRightPixel = partition.rightPixel;
    }
    i = innerTopPixel;
    j = innerLeftPixel;
}

std::pair<int, int> GridIterator::next() {
    std::pair<int, int> ret(i, j);
    if (i == -1) {
        return ret;
    }
    if ((i == innerBottomPixel && j == innerRightPixel) || innerFinished) {
        // this is the outer region
        if (innerFinished) {
            j++;
            communicate = false;
        } else {
            // when last pixel of the inner region has reached, we then start from top left of the outer region
            i = partition.topPixel;
            j = partition.leftPixel;
            innerFinished = true;
            communicate = partition.neighborCount > 0;
        }
        bool done = false;
        // this for loop searches for the next valid outer region cell
        for (; i <= partition.bottomPixel; i++) {
            for (; j <= partition.rightPixel; j++) {
                if (i >= innerTopPixel && i <= innerBottomPixel && j >= innerLeftPixel && j <= innerRightPixel) {
                    j = innerRightPixel;
                } else {
                    done = true;
                    break;
                }
            }
            if (done) {
                break;
            }
            j = partition.leftPixel;
        }
        if (!done) {
            i = -1;
        }
    } else {
        // this is the inner region
        j++;
        if (j > innerRightPixel) {
            i++;
            j = innerLeftPixel;
        }
    }
    return ret;
}

bool GridIterator::shouldCommunicate() {
    return communicate;
}
