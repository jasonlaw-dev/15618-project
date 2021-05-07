#ifndef PARTITION_INFO_H
#define PARTITION_INFO_H

struct PartitionInfo {
    // pixel index
    int topPixel;
    int bottomPixel;
    int leftPixel;
    int rightPixel;

    int neigbhors[9];
    int neighborCount;
};
#endif