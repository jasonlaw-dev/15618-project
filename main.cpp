#include <FreeImage.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <locale>
#include <assert.h>
#include <stdio.h>
#include "griditerator.h"
#include "partitioninfo.h"
#include <utility>

#include "mpi.h"

#define SYSEXPECT(expr)       \
    do {                      \
        if (!(expr)) {        \
            perror(__func__); \
            exit(1);          \
        }                     \
    } while (0)
#define error_exit(fmt, ...)                                        \
    do {                                                            \
        fprintf(stderr, "%s error: " fmt, __func__, ##__VA_ARGS__); \
        exit(1);                                                    \
    } while (0);

int BPP, WIDTH, HEIGHT, PITCH;
int BLOCK_ROWS, BLOCK_COLS;
int BLOCK_HEIGHT, BLOCK_WIDTH;
int NPROC = -1;
int PROCID = -1;
int WEAK = 25;
int STRONG = 255;
BYTE *OUTPUT_BUFFER = NULL;
BYTE *INPUT_BUFFER = NULL;
float *OUTPUT_BUFFER_FLOAT = NULL;
float *INPUT_BUFFER_FLOAT = NULL;
FREE_IMAGE_FORMAT IMAGE_FORMAT = FIF_UNKNOWN;

int TAG_OUTPUT_INIT = 500;
int TAG_OUTPUT_STEP0 = 1000;
int TAG_OUTPUT_STEP1 = 1001;
int TAG_OUTPUT_STEP2 = 1002;
int TAG_OUTPUT_STEP2_DIRECTION = 2002;
int TAG_OUTPUT_STEP4 = 1004;
int TAG_OUTPUT_FINAL = 1005;

double startTime;
double endTime;

const char* filepath = "";
int stage = 0;
bool useGaussian = false;
bool useGaussianTwoPass = false;
bool isOptimized = false;

PartitionInfo PARTITION;
PartitionInfo *PARTITIONS = NULL;

void setupBlockSize() {
    BLOCK_ROWS = 1;
    BLOCK_COLS = NPROC;

    int diff = abs(HEIGHT - (WIDTH + BLOCK_COLS - 1) / BLOCK_COLS);

    // Try to find the best combination of rows / cols such that a block is closer to a square.
    for (int r = 2; r <= NPROC; r++) {
        if (NPROC % r != 0) {
            continue;
        }
        int c = NPROC / r;
        int temp = abs((HEIGHT + r - 1) / r - (WIDTH + c - 1) / c);
        if (temp < diff) {
            diff = temp;
            BLOCK_ROWS = r;
            BLOCK_COLS = c;
        }
    }

    BLOCK_HEIGHT = (HEIGHT + BLOCK_ROWS - 1) / BLOCK_ROWS;
    BLOCK_WIDTH = (WIDTH + BLOCK_COLS - 1) / BLOCK_COLS;
}

PartitionInfo *getPartition(int procId) {
    PartitionInfo *partition = new PartitionInfo;

    int row = procId / BLOCK_COLS;
    int col = procId % BLOCK_COLS;

    partition->topPixel = row * BLOCK_HEIGHT;
    partition->bottomPixel = partition->topPixel + BLOCK_HEIGHT - 1;
    if (partition->bottomPixel >= HEIGHT) {
        partition->bottomPixel = HEIGHT - 1;
    }

    partition->leftPixel = col * BLOCK_WIDTH;
    partition->rightPixel = partition->leftPixel + BLOCK_WIDTH - 1;
    if (partition->rightPixel >= WIDTH) {
        partition->rightPixel = WIDTH - 1;
    }

    /*
        3x3 grid, NPROC = 9, procId = 2 (top right cell)
        Neighbors:
        -1      -1      -1
        1       -1      -1
        4       5       -1
    */
    partition->neighborCount = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == 1 && j == 1) {
                partition->neigbhors[i * 3 + j] = -1;
            } else {
                int neighbor = procId + BLOCK_COLS * (i - 1) + (j - 1);
                if (neighbor < 0 || neighbor >= NPROC) {
                    neighbor = -1;
                } else if (row + (i-1) != neighbor / BLOCK_COLS) {
                    neighbor = -1;
                }
                if (neighbor != -1) {
                    partition->neighborCount++;
                }
                partition->neigbhors[i * 3 + j] = neighbor;
            }
        }
    }
    return partition;
}


void getCommunicationInfo(int *info, PartitionInfo partition, int radius, int direction) {
    int width = 0;
    int height = 0;
    switch(direction) {
        case 0:
        case 2:
        case 6:
        case 8:
            width = height = radius;
            break;
        case 1:
        case 7:
            width = partition.rightPixel - partition.leftPixel + 1;
            height = radius;
            break;
        case 3:
        case 5:
            width = radius;
            height = partition.bottomPixel - partition.topPixel + 1;
            break;
    }
    int startI = 0;
    int startJ = 0;
    switch(direction) {
        case 0:
        case 1:
        case 3:
            startI = 0;
            startJ = 0;
            break;
        case 2:
        case 5:
            startI = 0;
            startJ = partition.rightPixel - partition.leftPixel + 1 - radius;
            break;
        case 6:
        case 7:
            startI = partition.bottomPixel - partition.topPixel + 1 - radius;
            startJ = 0;
            break;
        case 8:
            startI = partition.bottomPixel - partition.topPixel + 1 - radius;
            startJ = partition.rightPixel - partition.leftPixel + 1 - radius;
            break;
    }
    info[0] = width;
    info[1] = height;
    info[2] = startI;
    info[3] = startJ;
}

void getCommunicationInfoForInitialization(int *info, PartitionInfo partition, int radius) {
    int startI = std::max(partition.topPixel - radius, 0);
    int endI = std::min(partition.bottomPixel + radius, HEIGHT-1);
    int startJ = std::max(partition.leftPixel - radius, 0);
    int endJ = std::min(partition.rightPixel + radius, WIDTH-1);

    int width = endJ - startJ + 1;
    int height = endI - startI + 1;
    
    info[0] = width;
    info[1] = height;
    info[2] = startI;
    info[3] = startJ;
}

template <typename T>
void sendBlock(T *buffer, MPI_Datatype type, int dest, int tag, int radius, int direction) {
    MPI_Request *request = new MPI_Request;
    if (radius == 0) {
        MPI_Isend(buffer, BLOCK_WIDTH * BLOCK_HEIGHT, type, dest, tag, MPI_COMM_WORLD, request);
    } else {
        int info[4];
        getCommunicationInfo(info, PARTITION, radius, direction);
        int width = info[0];
        int height = info[1];
        int startI = info[2];
        int startJ = info[3];
        int k = 0;
        T *bufferToSend = new T[width * height];
        for (int i = startI; i < startI + height; i++) {
            std::copy(buffer + i * BLOCK_WIDTH + startJ, buffer + i * BLOCK_WIDTH + startJ + width, bufferToSend + k);
            k += width;
        }
        MPI_Isend(bufferToSend, width * height, type, dest, tag, MPI_COMM_WORLD, request);
    }
}

void sendBlockByte(int dest, int tag, int radius = 0, int direction = -1) {
    sendBlock<BYTE>(OUTPUT_BUFFER, MPI_BYTE, dest, tag, radius, direction);
}

void sendBlockFloat(int dest, int tag, int radius = 0, int direction = -1) {
    sendBlock<float>(OUTPUT_BUFFER_FLOAT, MPI_FLOAT, dest, tag, radius, direction);
}

template <typename T>
void copyBlockToBuffer(T *image, T *buffer, int radius) {
    if (radius == 0) {
        for (int i = PARTITION.topPixel; i <= PARTITION.bottomPixel; i++) {
            std::copy(image + i * PITCH + PARTITION.leftPixel, image + i * PITCH + PARTITION.rightPixel + 1, buffer + (i - PARTITION.topPixel) * BLOCK_WIDTH);
        }
    } else {
        for (int i = PARTITION.topPixel; i <= PARTITION.bottomPixel; i++) {
            if (i < PARTITION.topPixel + radius || i > PARTITION.bottomPixel - radius) {
                std::copy(image + i * PITCH + PARTITION.leftPixel, image + i * PITCH + PARTITION.rightPixel + 1, buffer + (i - PARTITION.topPixel) * BLOCK_WIDTH);
            } else {
                std::copy(image + i * PITCH + PARTITION.leftPixel, image + i * PITCH + PARTITION.leftPixel + radius, buffer + (i - PARTITION.topPixel) * BLOCK_WIDTH);
                std::copy(image + i * PITCH + PARTITION.rightPixel - radius + 1, image + i * PITCH + PARTITION.rightPixel + 1, buffer + (i - PARTITION.topPixel) * BLOCK_WIDTH + (PARTITION.rightPixel - PARTITION.leftPixel) - radius + 1);
            }
        }
    }
}

void copyBlockToBufferByte(BYTE *image, int radius = 0) {
    copyBlockToBuffer<BYTE>(image, OUTPUT_BUFFER, radius);
}

void copyBlockToBufferFloat(float *image, int radius = 0) {
    copyBlockToBuffer<float>(image, OUTPUT_BUFFER_FLOAT, radius);
}

template <typename T>
void recvBlocks(T *image, T *buffer, MPI_Datatype type, int recvCount, int tag, int radius) {
    MPI_Status *status = new MPI_Status;
    for (; recvCount > 0; recvCount--) {
        MPI_Recv(buffer, BLOCK_WIDTH * BLOCK_HEIGHT, type, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status);
        int source = status->MPI_SOURCE;
        if (radius == 0) {
            for (int i = PARTITIONS[source].topPixel; i <= PARTITIONS[source].bottomPixel; i++) {
                T *bufferStart = buffer + (i - PARTITIONS[source].topPixel) * BLOCK_WIDTH;
                T *bufferEnd = bufferStart + BLOCK_WIDTH;
                std::copy(bufferStart, bufferEnd, image + i * PITCH + PARTITIONS[source].leftPixel);
            }    
        } else {
            int direction;
            for (direction = 0; direction < 9; direction++) {
                if (PARTITIONS[source].neigbhors[direction] == PROCID) {
                    break;
                }
            }
            int info[4];
            getCommunicationInfo(info, PARTITIONS[source], radius, direction);
            int width = info[0];
            int height = info[1];
            int startI = info[2];
            int startJ = info[3];
            int k = 0;
            for (int i = startI; i < startI + height; i++) {
                std::copy(buffer + k, buffer + k + width, image + (PARTITIONS[source].topPixel + i) * PITCH + PARTITIONS[source].leftPixel + startJ);
                k += width;
            }
        }
    }
}

void recvBlocksByte(BYTE *image, int recvCount, int tag, int radius = 0) {
    recvBlocks<BYTE>(image, INPUT_BUFFER, MPI_BYTE, recvCount, tag, radius);
}

void recvBlocksFloat(float *image, int recvCount, int tag, int radius = 0) {
    recvBlocks<float>(image, INPUT_BUFFER_FLOAT, MPI_FLOAT, recvCount, tag, radius);
}

// https://github.com/aekanman/gaussian-blur/blob/master/helper.cpp
float* getGaussianFilterKernel(int diameter, float sigma) {
    float *filter = new float[diameter * diameter];
    int half = diameter / 2;
    float filterSum = 0.f; //for normalization
    for (int r = -half; r <= half; ++r) {
        for (int c = -half; c <= half; ++c) {
            float filterValue = expf( -(float)(c * c + r * r) / (2.f * sigma * sigma));
            filter[(r + half) * diameter + c + half] = filterValue;
            filterSum += filterValue;
        }
    }

    float normalizationFactor = 1.f / filterSum;
    for (int r = -half; r <= half; ++r) {
        for (int c = -half; c <= half; ++c) {
            filter[(r + half) * diameter + c + half] *= normalizationFactor;
        }
    }
    return filter;
}
float* getGaussianFilterKernel1D(int diameter, float sigma) {
    float *filter = new float[diameter];
    int half = diameter / 2;
    float filterSum = 0.f; //for normalization
    for (int r = -half; r <= half; ++r) {
        int c = -half;
        float filterValue = expf( -(float)(c * c + r * r) / (2.f * sigma * sigma));
        filter[r + half] = filterValue;
        filterSum += filterValue;
    }

    float normalizationFactor = 1.f / filterSum;
    for (int r = -half; r <= half; ++r) {
        filter[r + half] *= normalizationFactor;
    }
    return filter;
}

void applyGaussianFilter(BYTE *src, BYTE *dst, int diameter, float *filter) {
    int half = diameter / 2;
    for (int i = PARTITION.topPixel; i <= PARTITION.bottomPixel; i++) {
        for (int j = PARTITION.leftPixel; j <= PARTITION.rightPixel; j++) {
            float iFiltered = 0;
            int neighbor_i = 0;
            int neighbor_j = 0;

            for (int ii = 0; ii < diameter; ii++) {
                for (int jj = 0; jj < diameter; jj++) {
                    neighbor_i = i - (half - ii);
                    neighbor_j = j - (half - jj);
                    if (neighbor_i < 0) {
                        neighbor_i += diameter + 1;
                    } else if (neighbor_i >= HEIGHT) {
                        neighbor_i -= diameter + 1;
                    }
                    if (neighbor_j < 0) {
                        neighbor_j += diameter + 1;
                    } else if (neighbor_j >= WIDTH) {
                        neighbor_j -= diameter + 1;
                    }

                    iFiltered += src[neighbor_i * PITCH + neighbor_j] * filter[ii * diameter + jj];
                }
            }
            dst[i * PITCH + j] = iFiltered;
        }
    }
}
void applyGaussianFilterTwoPass(BYTE *src, BYTE *dst, BYTE *temp, int diameter, float *filter) {
    int half = diameter / 2;
    for (int i = PARTITION.topPixel; i <= PARTITION.bottomPixel; i++) {
        for (int j = PARTITION.leftPixel; j <= PARTITION.rightPixel; j++) {
            float iFiltered = 0;
            int neighbor_i = 0;
            int neighbor_j = 0;

            for (int ii = 0; ii < diameter; ii++) {
                int jj = half;
                neighbor_i = i - (half - ii);
                neighbor_j = j - (half - jj);
                if (neighbor_i < 0) {
                    neighbor_i += diameter + 1;
                } else if (neighbor_i >= HEIGHT) {
                    neighbor_i -= diameter + 1;
                }
                if (neighbor_j < 0) {
                    neighbor_j += diameter + 1;
                } else if (neighbor_j >= WIDTH) {
                    neighbor_j -= diameter + 1;
                }

                iFiltered += src[neighbor_i * PITCH + neighbor_j] * filter[ii];
            }
            temp[i * PITCH + j] = iFiltered;
        }
    }
    copyBlockToBufferByte(temp, isOptimized ? half : 0);
    for (int i = 0; i < 9; i++) {
        if (PARTITION.neigbhors[i] == -1) {
            continue;
        }
        sendBlockByte(PARTITION.neigbhors[i], TAG_OUTPUT_STEP0, isOptimized ? half : 0, i);
    }
    if (!isOptimized) {
        recvBlocksByte(temp, PARTITION.neighborCount, TAG_OUTPUT_STEP0);
    }

    GridIterator iter(PARTITION, half, PROCID);
    for (std::pair<int,int> ij = iter.next(); ij.first != -1; ij = iter.next()) {
        if (isOptimized && iter.shouldCommunicate()) {
            recvBlocksByte(temp, PARTITION.neighborCount, TAG_OUTPUT_STEP0, half);
        }
        int i = ij.first;
        int j = ij.second;
        float iFiltered = 0;
        int neighbor_i = 0;
        int neighbor_j = 0;

        for (int jj = 0; jj < diameter; jj++) {
            int ii = half;
            neighbor_i = i - (half - ii);
            neighbor_j = j - (half - jj);
            if (neighbor_i < 0) {
                neighbor_i += diameter + 1;
            } else if (neighbor_i >= HEIGHT) {
                neighbor_i -= diameter + 1;
            }
            if (neighbor_j < 0) {
                neighbor_j += diameter + 1;
            } else if (neighbor_j >= WIDTH) {
                neighbor_j -= diameter + 1;
            }

            iFiltered += temp[neighbor_i * PITCH + neighbor_j] * filter[jj];
        }
        dst[i * PITCH + j] = iFiltered;
    }
}

// https://github.com/anlcnydn/bilateral/blob/master/bilateral_filter.cpp
float distance(int x, int y, int i, int j) {
    return sqrtf(powf(x - i, 2) + powf(y - j, 2));
}

float gaussian(float x, double sigma) {
    return expf(-(powf(x, 2)) / (2 * powf(sigma, 2))) / (2 * M_PI * powf(sigma, 2));
}

void applyBilateralFilter(BYTE *src, BYTE *dst, int diameter, float sigmaI, float sigmaS) {
    int half = diameter / 2;
    for (int i = PARTITION.topPixel; i <= PARTITION.bottomPixel; i++) {
        for (int j = PARTITION.leftPixel; j <= PARTITION.rightPixel; j++) {
            float iFiltered = 0;
            float wP = 0;
            int neighbor_i = 0;
            int neighbor_j = 0;

            for (int ii = 0; ii < diameter; ii++) {
                for (int jj = 0; jj < diameter; jj++) {
                    neighbor_i = i - (half - ii);
                    neighbor_j = j - (half - jj);
                    if (neighbor_i < 0) {
                        neighbor_i += diameter + 1;
                    } else if (neighbor_i >= HEIGHT) {
                        neighbor_i -= diameter + 1;
                    }
                    if (neighbor_j < 0) {
                        neighbor_j += diameter + 1;
                    } else if (neighbor_j >= WIDTH) {
                        neighbor_j -= diameter + 1;
                    }

                    float gi = gaussian(src[neighbor_i * PITCH + neighbor_j] - src[i * PITCH + j], sigmaI);
                    float gs = gaussian(distance(i, j, neighbor_i, neighbor_j), sigmaS);
                    float w = gi * gs;
                    iFiltered = iFiltered + src[neighbor_i * PITCH + neighbor_j] * w;
                    wP = wP + w;
                }
            }
            iFiltered = iFiltered / wP;
            dst[i * PITCH + j] = iFiltered;
        }
    }
}

// https://github.com/brunokeymolen/canny/blob/master/canny.cpp
void applySobelFilter(BYTE *src, BYTE *gradient, float *direction) {
    int xFilter[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    int yFilter[3][3] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};
    float maxGradient = 0.f;
    float *gradientTemp = new float[PITCH * HEIGHT];
    GridIterator iter(PARTITION, isOptimized ? 1 : 0, PROCID);
    if (!isOptimized) {
        recvBlocksByte(src, PARTITION.neighborCount, TAG_OUTPUT_STEP1);
    }
    for (std::pair<int,int> ij = iter.next(); ij.first != -1; ij = iter.next()) {
        int i = ij.first;
        int j = ij.second;
        if (isOptimized && iter.shouldCommunicate()) {
            recvBlocksByte(src, PARTITION.neighborCount, TAG_OUTPUT_STEP1, 1);
        }
        if (i == 0 || i == HEIGHT - 1 || j == 0 || j == WIDTH - 1) {
            continue;
        }
        float xGradient = 0;
        float yGradient = 0;
        for (int ii = 0; ii < 3; ii++) {
            for (int jj = 0; jj < 3; jj+=2) {
                int neighbor_i = i + ii - 1;
                int neighbor_j = j + jj - 1;
                xGradient += xFilter[ii][jj] * src[neighbor_i * PITCH + neighbor_j];
            }
        }
        for (int ii = 0; ii < 3; ii+=2) {
            for (int jj = 0; jj < 3; jj++) {
                int neighbor_i = i + ii - 1;
                int neighbor_j = j + jj - 1;
                yGradient += yFilter[ii][jj] * src[neighbor_i * PITCH + neighbor_j];
            }
        }
        float grad = sqrtf(xGradient * xGradient + yGradient * yGradient);

        gradientTemp[i * PITCH + j] = grad;
        maxGradient = std::max(grad, maxGradient);

        direction[i * PITCH + j] = atan2f(yGradient, xGradient) * 180.f / M_PI;
        if (direction[i * PITCH + j] < 0) {
            direction[i * PITCH + j] += 180.f;
        }
    }
    float maxGradientToSend = maxGradient;
    MPI_Allreduce(&maxGradientToSend, &maxGradient, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

    iter = GridIterator(PARTITION, 0, PROCID);
    float factor = 1.f / maxGradient * 255.f;
    for (int i = PARTITION.topPixel; i <= PARTITION.bottomPixel; i++) {
        for (int j = PARTITION.leftPixel; j <= PARTITION.rightPixel; j++) {
            if (i == 0 || i == HEIGHT - 1 || j == 0 || j == WIDTH - 1) {
                continue;
            }
            gradient[i * PITCH + j] = gradientTemp[i * PITCH + j] * factor;
        }
    }
}

void applyNonMaxSuppressionAndThreshold(BYTE *src, BYTE *dst, float *direction, bool applyThreshold = true) {
    GridIterator iter(PARTITION, isOptimized ? 1 : 0, PROCID);
    if (!isOptimized) {
        recvBlocksByte(src, PARTITION.neighborCount, TAG_OUTPUT_STEP2);
        recvBlocksFloat(direction, PARTITION.neighborCount, TAG_OUTPUT_STEP2_DIRECTION);
    }
    const int highThreshold = 22; // 255 * 0.09f
    const int lowThreshold = 1; // highThreshold * 0.05f
    for (std::pair<int,int> ij = iter.next(); ij.first != -1; ij = iter.next()) {
        int i = ij.first;
        int j = ij.second;
        if (isOptimized && iter.shouldCommunicate()) {
            recvBlocksByte(src, PARTITION.neighborCount, TAG_OUTPUT_STEP2, 1);
            recvBlocksFloat(direction, PARTITION.neighborCount, TAG_OUTPUT_STEP2_DIRECTION, 1);
        }
        if (i == 0 || i == HEIGHT - 1 || j == 0 || j == WIDTH - 1) {
            continue;
        }
        int q = 255;
        int r = 255;

        float dir = direction[i * PITCH + j];

        if ((0 <= dir && dir < 22.5) || (157.5 <= dir && dir <= 180)) {
            q = src[i * PITCH + j + 1];
            r = src[i * PITCH + j - 1];
        } else if (22.5 <= dir && dir < 67.5) {
            q = src[(i + 1) * PITCH + j - 1];
            r = src[(i - 1) * PITCH + j + 1];
        } else if (67.5 <= dir && dir < 112.5) {
            q = src[(i + 1) * PITCH + j];
            r = src[(i - 1) * PITCH + j];
        } else if (112.5 <= dir && dir < 157.5) {
            q = src[(i - 1) * PITCH + j - 1];
            r = src[(i + 1) * PITCH + j + 1];
        }
        int val = src[i * PITCH + j];
        if (val >= q && val >= r) {
            if (applyThreshold) {
                if (val >= highThreshold) {
                    dst[i * PITCH + j] = STRONG;
                } else if (val >= lowThreshold) {
                    dst[i * PITCH + j] = WEAK;
                }
            } else {
                dst[i * PITCH + j] = val;
            }
        }
    }
}

void applyNonMaxSuppression(BYTE *src, BYTE *dst, float *direction) {
    applyNonMaxSuppressionAndThreshold(src, dst, direction, false);
}

void applyThreshold(BYTE *image) {
    const int highThreshold = 22; // 255 * 0.09f
    const int lowThreshold = 1; // highThreshold * 0.05f

    for (int i = PARTITION.topPixel; i <= PARTITION.bottomPixel; i++) {
        for (int j = PARTITION.leftPixel; j <= PARTITION.rightPixel; j++) {
            if (i == 0 || i == HEIGHT - 1 || j == 0 || j == WIDTH - 1) {
                continue;
            }
            int center = i * PITCH + j;
            if (image[center] >= highThreshold) {
                image[center] = STRONG;
            } else if (image[center] >= lowThreshold) {
                image[center] = WEAK;
            }
        }
    }
}

void applyHysteresis(BYTE *image) {
    GridIterator iter(PARTITION, isOptimized ? 1 : 0, PROCID);
    if (!isOptimized) {
        recvBlocksByte(image, PARTITION.neighborCount, TAG_OUTPUT_STEP4);
    }
    for (std::pair<int,int> ij = iter.next(); ij.first != -1; ij = iter.next()) {
        int i = ij.first;
        int j = ij.second;
        if (isOptimized && iter.shouldCommunicate()) {
            recvBlocksByte(image, PARTITION.neighborCount, TAG_OUTPUT_STEP4, 1);
        }
        if (i == 0 || i == HEIGHT - 1 || j == 0 || j == WIDTH - 1) {
            continue;
        }
        int bottom = (i + 1) * PITCH + j;
        int center = i * PITCH + j;
        int top = (i - 1) * PITCH + j;
        if (image[center] == WEAK) {
            if ((image[bottom - 1] == STRONG) || (image[bottom] == STRONG) 
            || (image[bottom + 1] == STRONG) || (image[center - 1] == STRONG) 
            || (image[center + 1] == STRONG) || (image[top - 1] == STRONG) 
            || (image[top] == STRONG) || (image[top + 1] == STRONG)) {
                image[center] = STRONG;
            } else {
                image[center] = 0;
            }
        }
    }
}

void printImageInfo(FIBITMAP *image) {
    std::cout << "BPP: " << BPP << " Width: " << WIDTH << " Height: " << HEIGHT << " Pitch: " << PITCH << " Red mask: " << FreeImage_GetRedMask(image) << std::endl;
}

void saveImage(const char *filepath, FIBITMAP *image, const char *stageName){
    std::string outpath;
    outpath = std::string(filepath); 
    outpath.insert(outpath.length() - 4, stageName);
    FreeImage_Save(IMAGE_FORMAT, image, outpath.c_str());
    std::cout << outpath << " complete" << std::endl;
}

void saveImage(const char *filepath, BYTE *dst, const char *stageName){
    FIBITMAP *image = FreeImage_ConvertFromRawBits(dst, WIDTH, HEIGHT, PITCH, BPP, 0, 0, 0, true);
    saveImage(filepath, image, stageName);
}

void aggregateOutputAndSaveImage(const char *filepath, BYTE *dst, const char *stageName) {
    if (PROCID != 0) {
        copyBlockToBufferByte(dst);
        sendBlockByte(0, TAG_OUTPUT_FINAL);
    } else {
        recvBlocksByte(dst, NPROC - 1, TAG_OUTPUT_FINAL);
        saveImage(filepath, dst, stageName);
    }
}

void start(){
    FIBITMAP *image = NULL;
    int info[4];
    if (PROCID == 0) {
        int len = strlen(filepath);
        if (len < 5){
            error_exit("image file not found \"%s\"\n", filepath);
        }
        const char *ext = filepath + len - 4;
        if (strcmp(ext, ".jpg") == 0 || strcmp(ext, ".JPG") == 0) {
            IMAGE_FORMAT = FIF_JPEG;
        } else if (strcmp(ext, ".png") == 0 || strcmp(ext, ".PNG") == 0) {
            IMAGE_FORMAT = FIF_PNG;
        } else if (strcmp(ext, ".bmp") == 0 || strcmp(ext, ".bmp") == 0) {
            IMAGE_FORMAT = FIF_BMP;
        } else {
            error_exit("unsupported image extension \"%s\"\n", ext);
        }
        image = FreeImage_Load(IMAGE_FORMAT, filepath);
        if(image == NULL){
            error_exit("image file not found \"%s\"\n", filepath);
        }
        image = FreeImage_ConvertToGreyscale(image);
        image = FreeImage_ConvertTo8Bits(image);

        info[0] = BPP = FreeImage_GetBPP(image);
        info[1] = WIDTH = FreeImage_GetWidth(image);
        info[2] = HEIGHT = FreeImage_GetHeight(image);
        info[3] = PITCH = FreeImage_GetPitch(image);
    }
    MPI_Bcast(info, 4, MPI_INT, 0, MPI_COMM_WORLD);

    if (PROCID != 0) {
        BPP = info[0];
        WIDTH = info[1];
        HEIGHT = info[2];
        PITCH = info[3];
    }

    BYTE *src = new BYTE[PITCH * HEIGHT];
    BYTE *dst = new BYTE[PITCH * HEIGHT];
    PARTITIONS = new PartitionInfo[NPROC];
    setupBlockSize();

    for (int i = 0; i < NPROC; i++) {
        PARTITIONS[i] = *getPartition(i);
        if (i == PROCID) {
            PARTITION = PARTITIONS[i];
        }
    }
    OUTPUT_BUFFER = new BYTE[BLOCK_WIDTH * BLOCK_HEIGHT];
    INPUT_BUFFER = new BYTE[BLOCK_WIDTH * BLOCK_HEIGHT];
    OUTPUT_BUFFER_FLOAT = new float[BLOCK_WIDTH * BLOCK_HEIGHT];
    INPUT_BUFFER_FLOAT = new float[BLOCK_WIDTH * BLOCK_HEIGHT];

    if (PROCID == 0) {
        FreeImage_ConvertToRawBits(src, image, PITCH, BPP, 0, 0, 0, true);
    }
    
    int blurDiameter = 15;
    
    // MPI_Barrier(MPI_COMM_WORLD);
    startTime = MPI_Wtime();
    if (isOptimized) {
        if (PROCID == 0) {
            MPI_Request *request = new MPI_Request;
            for (int proc = 1; proc < NPROC; proc++) {
                getCommunicationInfoForInitialization(info, PARTITIONS[proc], blurDiameter / 2);
                int width = info[0];
                int height = info[1];
                int startI = info[2];
                int startJ = info[3];
                int k = 0;
                BYTE *buffer = new BYTE[width * height];
                for (int i = startI; i < startI + height; i++) {
                    std::copy(src + i * PITCH + startJ, src + i * PITCH + startJ + width, buffer + k);
                    k += width;
                }
                MPI_Isend(buffer, width * height, MPI_BYTE, proc, TAG_OUTPUT_INIT, MPI_COMM_WORLD, request);
            }
        } else {
            getCommunicationInfoForInitialization(info, PARTITION, blurDiameter / 2);
            int width = info[0];
            int height = info[1];
            int startI = info[2];
            int startJ = info[3];
            BYTE *buffer = new BYTE[width * height];
            MPI_Status status;
            MPI_Recv(buffer, width * height, MPI_BYTE, 0, TAG_OUTPUT_INIT, MPI_COMM_WORLD, &status);
            int k = 0;
            for (int i = startI; i < startI + height; i++) {
                std::copy(buffer + k, buffer + k + width, src + i * PITCH + startJ);
                k += width;
            }
        }
    } else {
        MPI_Bcast(src, PITCH * HEIGHT, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    
    if(useGaussian){
        if (useGaussianTwoPass) {
            BYTE *temp = new BYTE[PITCH * HEIGHT];
            float *gaussianFilter1D = getGaussianFilterKernel1D(blurDiameter, 2.f);
            applyGaussianFilterTwoPass(src, dst, temp, blurDiameter, gaussianFilter1D);
        } else {
            float *gaussianFilter = getGaussianFilterKernel(blurDiameter, 2.f);
            applyGaussianFilter(src, dst, blurDiameter, gaussianFilter);
        }
    } else {
        applyBilateralFilter(src, dst, blurDiameter, 30.f, 20.f);
    }
    copyBlockToBufferByte(dst, isOptimized ? 1 : 0);
    for (int i = 0; i < 9; i++) {
        if (PARTITION.neigbhors[i] == -1) {
            continue;
        }
        sendBlockByte(PARTITION.neigbhors[i], TAG_OUTPUT_STEP1, isOptimized ? 1 : 0, i);
    }

    if(stage == 1){
        printf("stage 1: elapsed time for proc %d: %f\n", PROCID, MPI_Wtime() - startTime);
        if(useGaussian){
            aggregateOutputAndSaveImage(filepath, dst, "-1-gaussian");
        } else {
            aggregateOutputAndSaveImage(filepath, dst, "-1-bilateral");
        }
        return;
    }

    double startTime2 = MPI_Wtime();

    src = dst;
    dst = new BYTE[PITCH * HEIGHT]();
    float *direction = new float[PITCH * HEIGHT];

    applySobelFilter(src, dst, direction);

    copyBlockToBufferByte(dst, isOptimized ? 1 : 0);
    copyBlockToBufferFloat(direction, isOptimized ? 1 : 0);
    for (int i = 0; i < 9; i++) {
        if (PARTITION.neigbhors[i] == -1) {
            continue;
        }
        sendBlockByte(PARTITION.neigbhors[i], TAG_OUTPUT_STEP2, isOptimized ? 1 : 0, i);
        sendBlockFloat(PARTITION.neigbhors[i], TAG_OUTPUT_STEP2_DIRECTION, isOptimized ? 1 : 0, i);   
    }

    if(stage == 2){
        printf("stage 2: elapsed time for proc %d: %f\n", PROCID, MPI_Wtime() - startTime2);
        aggregateOutputAndSaveImage(filepath, dst, "-2-sobel");
        return;
    }

    double startTime3 = MPI_Wtime();
    src = dst;
    dst = new BYTE[PITCH * HEIGHT]();

    if (isOptimized) {
        applyNonMaxSuppressionAndThreshold(src, dst, direction);
    } else {
        applyNonMaxSuppression(src, dst, direction);

        if(stage == 3){
            printf("stage 3: elapsed time for proc %d: %f\n", PROCID, MPI_Wtime() - startTime3);
            aggregateOutputAndSaveImage(filepath, dst, "-3-nonmax");
            return;
        }
        double startTime4 = MPI_Wtime();

        applyThreshold(dst);

        if(stage == 4){
            printf("stage 4: elapsed time for proc %d: %f\n", PROCID, MPI_Wtime() - startTime4);
            aggregateOutputAndSaveImage(filepath, dst, "-4-thres");
            return;
        }
    }
    
    double startTime5 = MPI_Wtime();

    copyBlockToBufferByte(dst, isOptimized ? 1 : 0);
    for (int i = 0; i < 9; i++) {
        if (PARTITION.neigbhors[i] == -1) {
            continue;
        }
        sendBlockByte(PARTITION.neigbhors[i], TAG_OUTPUT_STEP4, isOptimized ? 1 : 0, i);
    }

    applyHysteresis(dst);
    
    if(stage == 5){
        printf("stage 5: elapsed time for proc %d: %f\n", PROCID, MPI_Wtime() - startTime5);
        aggregateOutputAndSaveImage(filepath, dst, "-5-hyster");
        return;
    }
    printf("stage 5: elapsed time for proc %d: %f\n", PROCID, MPI_Wtime() - startTime);
    aggregateOutputAndSaveImage(filepath, dst, "-5-hyster");
}

int main(int argc, char* argv[]) {
    int procID;
    int nproc;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Run computation
    PROCID = procID;
    NPROC = nproc;

    if (argc < 3) {
        std::cout << "Usage: mpirun -n <# of cores> ./main -f <image file>" << std::endl;
        exit(-1);
    }
    int ch = 0;
    while ((ch = getopt(argc, argv, "f:s:got")) != -1)
    {
        switch (ch)
        {
            case 'f':
                filepath = optarg;
                break;
            case 's':
                stage = atoi(optarg);
                break;
            case 'g':
                useGaussian = true;
                break;
            case 'o':
                isOptimized = true;
                break;
            case 't':
                useGaussianTwoPass = true;
        }
    }
    start();
    
    // Compute running time
    MPI_Finalize();
    return 0;
}