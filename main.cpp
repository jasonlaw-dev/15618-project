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

int bpp, width, height, pitch;
int NPROC = -1;
int PROCID = -1;
int WEAK = 0;
int STRONG = 255;

// https://github.com/anlcnydn/bilateral/blob/master/bilateral_filter.cpp
float distance(int x, int y, int i, int j) {
    return sqrtf(powf(x - i, 2) + powf(y - j, 2));
}

float gaussian(float x, double sigma) {
    return expf(-(powf(x, 2)) / (2 * powf(sigma, 2))) / (2 * M_PI * powf(sigma, 2));
}

void applyBilateralFilter(BYTE *src, BYTE *dst, int imageWidth, int imageHeight, int imagePitch, int diameter, float sigmaI, float sigmaS) {
    for (int i = 0; i < imageHeight; i++) {
        for (int j = 0; j < imageWidth; j++) {
            float iFiltered = 0;
            float wP = 0;
            int neighbor_i = 0;
            int neighbor_j = 0;
            int half = diameter / 2;

            for (int ii = 0; ii < diameter; ii++) {
                for (int jj = 0; jj < diameter; jj++) {
                    neighbor_i = i - (half - ii);
                    neighbor_j = j - (half - jj);
                    if (neighbor_i < 0) {
                        neighbor_i += diameter + 1;
                    } else if (neighbor_i >= imageHeight) {
                        neighbor_i -= diameter + 1;
                    }
                    if (neighbor_j < 0) {
                        neighbor_j += diameter + 1;
                    } else if (neighbor_j >= imageWidth) {
                        neighbor_j -= diameter + 1;
                    }

                    float gi = gaussian(src[neighbor_i * imagePitch + neighbor_j] - src[i * imagePitch + j], sigmaI);
                    float gs = gaussian(distance(i, j, neighbor_i, neighbor_j), sigmaS);
                    float w = gi * gs;
                    iFiltered = iFiltered + src[neighbor_i * imagePitch + neighbor_j] * w;
                    wP = wP + w;
                }
            }
            iFiltered = iFiltered / wP;
            dst[i * imagePitch + j] = iFiltered;
        }
    }
}

// https://github.com/brunokeymolen/canny/blob/master/canny.cpp
void applySobelFilter(BYTE *src, BYTE *gradient, float *direction, int imageWidth, int imageHeight, int imagePitch) {
    int xFilter[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    int yFilter[3][3] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};
    float maxGradient = 0.f;
    float *gradientTemp = new float[imagePitch * imageHeight];
    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            float xGradient = 0;
            float yGradient = 0;
            for (int ii = 0; ii < 3; ii++) {
                for (int jj = 0; jj < 3; jj++) {
                    int neighbor_i = i + ii - 1;
                    int neighbor_j = j + jj - 1;
                    xGradient += xFilter[ii][jj] * src[neighbor_i * imagePitch + neighbor_j];
                    yGradient += yFilter[ii][jj] * src[neighbor_i * imagePitch + neighbor_j];
                }
            }
            float grad = sqrtf(xGradient * xGradient + yGradient * yGradient);

            gradientTemp[i * imagePitch + j] = grad;
            maxGradient = std::max(grad, maxGradient);

            direction[i * imagePitch + j] = atan2f(yGradient, xGradient) * 180.f / M_PI;
            if (direction[i * imagePitch + j] < 0) {
                direction[i * imagePitch + j] += 180.f;
            }
        }
    }

    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            gradient[i * imagePitch + j] = gradientTemp[i * imagePitch + j] / maxGradient * 255;
        }
    }
}

void applyNonMaxSuppression(BYTE *src, BYTE *dst, float *direction, int imageWidth, int imageHeight, int imagePitch) {
    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            int q = 255;
            int r = 255;

            if ((0 <= direction[i * imagePitch + j] && direction[i * imagePitch + j] < 22.5) || (157.5 <= direction[i * imagePitch + j] && direction[i * imagePitch + j] <= 180)) {
                q = src[i * imagePitch + j + 1];
                r = src[i * imagePitch + j - 1];
            } else if (22.5 <= direction[i * imagePitch + j] && direction[i * imagePitch + j] < 67.5) {
                q = src[(i + 1) * imagePitch + j - 1];
                r = src[(i - 1) * imagePitch + j + 1];
            } else if (67.5 <= direction[i * imagePitch + j] && direction[i * imagePitch + j] < 112.5) {
                q = src[(i + 1) * imagePitch + j];
                r = src[(i - 1) * imagePitch + j];
            } else if (112.5 <= direction[i * imagePitch + j] && direction[i * imagePitch + j] < 157.5) {
                q = src[(i - 1) * imagePitch + j - 1];
                r = src[(i + 1) * imagePitch + j + 1];
            }
            if ((src[i * imagePitch + j] >= q) && (src[i * imagePitch + j] >= r)) {
                dst[i * imagePitch + j] = src[i * imagePitch + j];
            } else {
                dst[i * imagePitch + j] = 0;
            }
        }
    }
}

void applyThreshold(BYTE *image, int imageWidth, int imageHeight, int imagePitch) {
    float lowThresholdRatio = 0.05f;
    float highThresholdRatio = 0.09f;
    int maxPixel = 0;
    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            if (image[i * imagePitch + j] > maxPixel) {
                maxPixel = image[i * imagePitch + j];
            }
        }
    }
    int highThreshold = maxPixel * highThresholdRatio;
    int lowThreshold = highThreshold * lowThresholdRatio;
    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            if (image[i * imagePitch + j] >= highThreshold) {
                image[i * imagePitch + j] = STRONG;
            } else if (image[i * imagePitch + j] >= lowThreshold) {
                image[i * imagePitch + j] = WEAK;
            } else {
                image[i * imagePitch + j] = 0;
            }
        }
    }
}

void applyHysteresis(BYTE *image, int imageWidth, int imageHeight, int imagePitch) {
    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            if (image[i * imagePitch+ j] == WEAK) {
                if ((image[(i + 1) * imagePitch + j - 1] == STRONG) || (image[(i + 1) * imagePitch + j] == STRONG) 
                || (image[(i + 1) * imagePitch + j + 1] == STRONG) || (image[i * imagePitch + j - 1] == STRONG) 
                || (image[i * imagePitch + j + 1] == STRONG) || (image[(i - 1) * imagePitch + j - 1] == STRONG) 
                || (image[(i - 1) * imagePitch + j] == STRONG) || (image[(i - 1) * imagePitch + j + 1] == STRONG)) {
                    image[i* imagePitch+ j] = STRONG;
                } else {
                    image[i* imagePitch+ j] = 0;
                }
            }
        }
    }
}

void printImageInfo(FIBITMAP *image) {
    std::cout << "BPP: " << bpp << " Width: " << width << " Height: " << height << " Pitch: " << pitch << " Red mask: " << FreeImage_GetRedMask(image) << std::endl;
}

void saveImage(const char *filepath, FIBITMAP *image, const char *stageName){
    std::string outpath;
    outpath = std::string(filepath); 
    outpath.insert(outpath.length() - 4, stageName);
    FreeImage_Save(FIF_PNG, image, outpath.c_str());
    std::cout << outpath << " complete" << std::endl;
}

void saveImage(const char *filepath, BYTE *dst, const char *stageName){
    FIBITMAP *image = FreeImage_ConvertFromRawBits(dst, width, height, pitch, bpp, 0, 0, 0, true);
    saveImage(filepath, image, stageName);
}

void start(int procID, int nproc, int argc, char* argv[]){

    PROCID = procID;
    NPROC = nproc;
    //const int root = 0;  // Set the rank 0 process as the root process

    // assign tasks
    const char* filepath = "";

    if(argc < 3){
        std::cout << "Usage: mpirun -n <# of cores> ./main -f <image file>" << std::endl;
        exit(-1);
    }
    if(strcmp(argv[1], "-f") == 0){
        filepath = argv[2];
    }
    else{
        std::cout << "Usage: mpirun -n <# of cores> ./main -f <image file>" << std::endl;
        exit(-1);
    }

    FIBITMAP *image = FreeImage_Load(FIF_PNG, filepath);
    if(image == NULL){
        error_exit("image file not found \"%s\"\n", filepath);
    }
    image = FreeImage_ConvertToGreyscale(image);
    image = FreeImage_ConvertTo8Bits(image);

    bpp = FreeImage_GetBPP(image);
    width = FreeImage_GetWidth(image);
    height = FreeImage_GetHeight(image);
    pitch = FreeImage_GetPitch(image);

    printImageInfo(image);

    BYTE *src = new BYTE[pitch * height];
    BYTE *dst = new BYTE[pitch * height]();

    FreeImage_ConvertToRawBits(src, image, pitch, bpp, 0, 0, 0, true);

    applyBilateralFilter(src, dst, width, height, pitch, 9, 30.f, 20.f);
    saveImage(filepath, dst, "-1-bilateral");

    src = dst;
    dst = new BYTE[pitch * height]();
    float *direction = new float[width * height];

    applySobelFilter(src, dst, direction, width, height, pitch);
    saveImage(filepath, dst, "-2-sobel");

    src = dst;
    dst = new BYTE[pitch * height]();

    applyNonMaxSuppression(src, dst, direction, width, height, pitch);
    saveImage(filepath, dst, "-3-nonmax");

    applyThreshold(dst, width, height, pitch);
    saveImage(filepath, dst, "-4-thres");

    // src = dst;
    // dst = new BYTE[width * height]();

    // applyHysteresis(dst, width,height, pitch);
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
    start(procID, nproc, argc, argv);
    endTime = MPI_Wtime();

    // Compute running time
    MPI_Finalize();
    printf("elapsed time for proc %d: %f\n", procID, endTime - startTime);
    return 0;
}
