#include <FreeImage.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>

#include <assert.h>
#include <stdio.h>

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
int WEAK = 0;
int STRONG = 255;
BYTE *OUTPUT_BUFFER = NULL;

int TAG_OUTPUT = 1001;

struct PartitionInfo {
    // pixel index
    int topPixel;
    int bottomPixel;
    int leftPixel;
    int rightPixel;

    int neigbhors[8];
    int neighborCount;
};

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

    // printf("NPROC: %d\n", NPROC);
    // printf("PROCID: %d\n", procId);
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         printf("%d\t", partition->neigbhors[i*3+j]);
    //     }
    //     printf("\n");
    // }

    return partition;
}

void sendOutputBlock(BYTE *image) {
    for (int i = PARTITION.topPixel; i <= PARTITION.bottomPixel; i++) {
        std::copy(image + i * PITCH + PARTITION.leftPixel, image + i * PITCH + PARTITION.rightPixel + 1, OUTPUT_BUFFER + (i - PARTITION.topPixel) * BLOCK_WIDTH);
    }
    MPI_Request *request = new MPI_Request;
    MPI_Isend(OUTPUT_BUFFER, BLOCK_WIDTH * BLOCK_HEIGHT, MPI_BYTE, 0, TAG_OUTPUT, MPI_COMM_WORLD, request);
}

void recvOutputBlocks(BYTE *image) {
    MPI_Status status;
    for (int x = 0; x < NPROC - 1; x++) {
        MPI_Recv(OUTPUT_BUFFER, BLOCK_WIDTH * BLOCK_HEIGHT, MPI_BYTE, MPI_ANY_SOURCE, TAG_OUTPUT, MPI_COMM_WORLD, &status);
        int source = status.MPI_SOURCE;
        for (int i = PARTITIONS[source].topPixel; i <= PARTITIONS[source].bottomPixel; i++) {
            BYTE *bufferStart = OUTPUT_BUFFER + (i - PARTITIONS[source].topPixel) * BLOCK_WIDTH;
            BYTE *bufferEnd = bufferStart + BLOCK_WIDTH;
            std::copy(bufferStart, bufferEnd, image + i * PITCH + PARTITIONS[source].leftPixel);
        }
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
    for (int i = 1; i < HEIGHT - 1; i++) {
        for (int j = 1; j < WIDTH - 1; j++) {
            float xGradient = 0;
            float yGradient = 0;
            for (int ii = 0; ii < 3; ii++) {
                for (int jj = 0; jj < 3; jj++) {
                    int neighbor_i = i + ii - 1;
                    int neighbor_j = j + jj - 1;
                    xGradient += xFilter[ii][jj] * src[neighbor_i * PITCH + neighbor_j];
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
    }

    for (int i = 1; i < HEIGHT - 1; i++) {
        for (int j = 1; j < WIDTH - 1; j++) {
            gradient[i * PITCH + j] = gradientTemp[i * PITCH + j] / maxGradient * 255;
        }
    }
}

void applyNonMaxSuppression(BYTE *src, BYTE *dst, float *direction) {
    for (int i = 1; i < HEIGHT - 1; i++) {
        for (int j = 1; j < WIDTH - 1; j++) {
            int q = 255;
            int r = 255;

            if ((0 <= direction[i * PITCH + j] && direction[i * PITCH + j] < 22.5) || (157.5 <= direction[i * PITCH + j] && direction[i * PITCH + j] <= 180)) {
                q = src[i * PITCH + j + 1];
                r = src[i * PITCH + j - 1];
            } else if (22.5 <= direction[i * PITCH + j] && direction[i * PITCH + j] < 67.5) {
                q = src[(i + 1) * PITCH + j - 1];
                r = src[(i - 1) * PITCH + j + 1];
            } else if (67.5 <= direction[i * PITCH + j] && direction[i * PITCH + j] < 112.5) {
                q = src[(i + 1) * PITCH + j];
                r = src[(i - 1) * PITCH + j];
            } else if (112.5 <= direction[i * PITCH + j] && direction[i * PITCH + j] < 157.5) {
                q = src[(i - 1) * PITCH + j - 1];
                r = src[(i + 1) * PITCH + j + 1];
            }
            if ((src[i * PITCH + j] >= q) && (src[i * PITCH + j] >= r)) {
                dst[i * PITCH + j] = src[i * PITCH + j];
            } else {
                dst[i * PITCH + j] = 0;
            }
        }
    }
}

void applyThreshold(BYTE *image) {
    float lowThresholdRatio = 0.05f;
    float highThresholdRatio = 0.09f;
    int maxPixel = 0;
    for (int i = 1; i < HEIGHT - 1; i++) {
        for (int j = 1; j < WIDTH - 1; j++) {
            if (image[i * PITCH + j] > maxPixel) {
                maxPixel = image[i * PITCH + j];
            }
        }
    }
    int highThreshold = maxPixel * highThresholdRatio;
    int lowThreshold = highThreshold * lowThresholdRatio;
    for (int i = 1; i < HEIGHT - 1; i++) {
        for (int j = 1; j < WIDTH - 1; j++) {
            if (image[i * PITCH + j] >= highThreshold) {
                image[i * PITCH + j] = STRONG;
            } else if (image[i * PITCH + j] >= lowThreshold) {
                image[i * PITCH + j] = WEAK;
            } else {
                image[i * PITCH + j] = 0;
            }
        }
    }
}

void applyHysteresis(BYTE *image) {
    for (int i = 1; i < HEIGHT - 1; i++) {
        for (int j = 1; j < WIDTH - 1; j++) {
            if (image[i * PITCH+ j] == WEAK) {
                if ((image[(i + 1) * PITCH + j - 1] == STRONG) || (image[(i + 1) * PITCH + j] == STRONG) 
                || (image[(i + 1) * PITCH + j + 1] == STRONG) || (image[i * PITCH + j - 1] == STRONG) 
                || (image[i * PITCH + j + 1] == STRONG) || (image[(i - 1) * PITCH + j - 1] == STRONG) 
                || (image[(i - 1) * PITCH + j] == STRONG) || (image[(i - 1) * PITCH + j + 1] == STRONG)) {
                    image[i* PITCH+ j] = STRONG;
                } else {
                    image[i* PITCH+ j] = 0;
                }
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
    FreeImage_Save(FIF_PNG, image, outpath.c_str());
    std::cout << outpath << " complete" << std::endl;
}

void saveImage(const char *filepath, BYTE *dst, const char *stageName){
    FIBITMAP *image = FreeImage_ConvertFromRawBits(dst, WIDTH, HEIGHT, PITCH, BPP, 0, 0, 0, true);
    saveImage(filepath, image, stageName);
}

void start(int argc, char* argv[]){
    const char* filepath = "";

    if (argc < 3) {
        std::cout << "Usage: mpirun -n <# of cores> ./main -f <image file>" << std::endl;
        exit(-1);
    }
    if (strcmp(argv[1], "-f") == 0) {
        filepath = argv[2];
    } else {
        std::cout << "Usage: mpirun -n <# of cores> ./main -f <image file>" << std::endl;
        exit(-1);
    }

    FIBITMAP *image = NULL;
    int info[4];
    if (PROCID == 0) {
        image = FreeImage_Load(FIF_PNG, filepath);
        if(image == NULL){
            error_exit("image file not found \"%s\"\n", filepath);
        }
        image = FreeImage_ConvertToGreyscale(image);
        image = FreeImage_ConvertTo8Bits(image);

        info[0] = BPP = FreeImage_GetBPP(image);
        info[1] = WIDTH = FreeImage_GetWidth(image);
        info[2] = HEIGHT = FreeImage_GetHeight(image);
        info[3] = PITCH = FreeImage_GetPitch(image);

        printImageInfo(image);
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

    if (PROCID == 0) {
        FreeImage_ConvertToRawBits(src, image, PITCH, BPP, 0, 0, 0, true);
    }
    MPI_Bcast(src, PITCH * HEIGHT, MPI_BYTE, 0, MPI_COMM_WORLD);

    applyBilateralFilter(src, dst, 9, 30.f, 20.f);
    if (PROCID != 0) {
        sendOutputBlock(dst);
    } else {
        recvOutputBlocks(dst);
    }

    if (PROCID == 0) {
        saveImage(filepath, dst, "-1-bilateral");
    }

    // src = dst;
    // dst = new BYTE[PITCH * HEIGHT]();
    // float *direction = new float[PITCH * HEIGHT];

    // applySobelFilter(src, dst, direction);
    // // saveImage(filepath, dst, "-2-sobel");

    // src = dst;
    // dst = new BYTE[PITCH * HEIGHT]();

    // applyNonMaxSuppression(src, dst, direction);
    // // saveImage(filepath, dst, "-3-nonmax");

    // applyThreshold(dst);
    // saveImage(filepath, dst, "-4-thres");

    // src = dst;
    // dst = new BYTE[PITCH * HEIGHT]();

    // applyHysteresis(dst);
    // saveImage(filepath, dst, "-5-hyster");
}

int main(int argc, char* argv[]) {
    int procID;
    int nproc;
    double startTime;
    double endTime;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Run computation
    startTime = MPI_Wtime();
    PROCID = procID;
    NPROC = nproc;
    start(argc, argv);
    endTime = MPI_Wtime();

    // Compute running time
    MPI_Finalize();
    printf("elapsed time for proc %d: %f\n", procID, endTime - startTime);
    return 0;
}
