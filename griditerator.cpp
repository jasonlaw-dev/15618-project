#include <algorithm>
#include "griditerator.h"
#include <iostream>

GridIterator::GridIterator(PartitionInfo partition, int borderWidth, int pitch): partition(partition), borderWidth(borderWidth), pitch(pitch), innerFinished(false) {
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
    std::cout << "top: " <<  innerTopPixel << " bottom: "<< innerBottomPixel << " left: " << innerLeftPixel << " right: " << innerRightPixel << std::endl;
    i = innerTopPixel;
    j = innerLeftPixel;
}

bool GridIterator::hasNext() {
    return !(i == -1 && j == -1);
}

std::pair<int, int> GridIterator::next() {
    std::pair<int, int> ret(i, j);
    // std::cout << "i: "  << i << " j: " << j << std::endl;

    if ((i == innerBottomPixel && j == innerRightPixel) || innerFinished) {
        if (innerFinished) {
            j++;
        } else {
            i = partition.topPixel;
            j = partition.leftPixel;
            innerFinished = true;
        }
        bool done = false;
        for (; !done && i <= partition.bottomPixel; i++) {
            for (; !done && j <= partition.rightPixel; j++) {
                if (!(i >= innerTopPixel && i <= innerBottomPixel && j >= innerLeftPixel && j <= innerRightPixel)) {
                    done = true;
                }
            }
            if (!done) {
                j = partition.leftPixel;
            }
        }
        if (!done) {
            i = j = -1;
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
