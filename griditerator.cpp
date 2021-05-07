#include <algorithm>
#include "griditerator.h"
#include <iostream>
#include <stdio.h>

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
    // std::cout << "top: " <<  innerTopPixel << " bottom: "<< innerBottomPixel << " left: " << innerLeftPixel << " right: " << innerRightPixel << std::endl;
    i = innerTopPixel;
    j = innerLeftPixel;
}

std::pair<int, int> GridIterator::next() {
    std::pair<int, int> ret(i, j);
    if (i == -1) {
        return ret;
    }
    if ((i == innerBottomPixel && j == innerRightPixel) || innerFinished) {
        if (innerFinished) {
            j++;
            communicate = false;
        } else {
            i = partition.topPixel;
            j = partition.leftPixel;
            innerFinished = true;
            communicate = partition.neighborCount > 0;
        }
        bool done = false;
        for (; i <= partition.bottomPixel; i++) {
            for (; j <= partition.rightPixel; j++) {
                if (!(i >= innerTopPixel && i <= innerBottomPixel && j >= innerLeftPixel && j <= innerRightPixel)) {
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
