#include <FreeImage.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>

int bpp, width, height, pitch;

// https://github.com/anlcnydn/bilateral/blob/master/bilateral_filter.cpp
float distance(int x, int y, int i, int j) {
    return sqrtf(powf(x - i, 2) + powf(y - j, 2));
}

float gaussian(float x, double sigma) {
    return expf(-(powf(x, 2)) / (2 * powf(sigma, 2))) / (2 * M_PI * powf(sigma, 2));
}

void applyBilateralFilterPixel(BYTE *src, BYTE *dst, int imageWidth, int imageHeight, int diameter, float sigmaI, float sigmaS, int i, int j) {
    float iFiltered = 0;
    float wP = 0;
    int neighbor_i = 0;
    int neighbor_j = 0;
    int half = diameter / 2;

    for (int ii = 0; ii < diameter; ii++) {
        for (int jj = 0; jj < diameter; jj++) {
            neighbor_i = i - (half - ii);
            neighbor_j = j - (half - jj);
            float gi = gaussian(src[neighbor_i * imageWidth + neighbor_j] - src[i * imageWidth + j], sigmaI);
            float gs = gaussian(distance(i, j, neighbor_i, neighbor_j), sigmaS);
            float w = gi * gs;
            iFiltered = iFiltered + src[neighbor_i * imageWidth + neighbor_j] * w;
            wP = wP + w;
        }
    }
    iFiltered = iFiltered / wP;
    dst[i * imageWidth + j] = iFiltered;
}

void applyBilateralFilter(BYTE *src, BYTE *dst, int imageWidth, int imageHeight, int diameter, float sigmaI, float sigmaS) {
    int half = diameter / 2;

    for (int i = half; i < imageHeight - half; i++) {
        for (int j = half; j < imageWidth - half; j++) {
            applyBilateralFilterPixel(src, dst, imageWidth, imageHeight, diameter, sigmaI, sigmaS, i, j);
        }
    }
}

// https://github.com/brunokeymolen/canny/blob/master/canny.cpp
void applySobelFilter(BYTE *src, BYTE *gradient, float *direction, int imageWidth, int imageHeight) {
    int xFilter[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    int yFilter[3][3] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};
    float maxGradient = 0.f;
    float *gradientTemp = new float[imageWidth * imageHeight];
    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            float xGradient = 0;
            float yGradient = 0;
            for (int ii = 0; ii < 3; ii++) {
                for (int jj = 0; jj < 3; jj++) {
                    int neighbor_i = i + ii - 1;
                    int neighbor_j = j + jj - 1;
                    xGradient += xFilter[ii][jj] * src[neighbor_i * imageWidth + neighbor_j];
                    yGradient += yFilter[ii][jj] * src[neighbor_i * imageWidth + neighbor_j];
                }
            }
            float grad = sqrtf(xGradient * xGradient + yGradient * yGradient);

            gradientTemp[i * imageWidth + j] = grad;
            maxGradient = std::max(grad, maxGradient);

            direction[i * imageWidth + j] = atan2f(yGradient, xGradient) * 180.f / M_PI;
            if (direction[i * imageWidth + j] < 0) {
                direction[i * imageWidth + j] += 180.f;
            }
        }
    }

    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            gradient[i * imageWidth + j] = gradientTemp[i * imageWidth + j] / maxGradient * 255;
        }
    }
}

void applyNonMaxSuppression(BYTE *src, BYTE *dst, float *direction, int imageWidth, int imageHeight) {
    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            int q = 255;
            int r = 255;

            if ((0 <= direction[i * imageWidth + j] && direction[i * imageWidth + j] < 22.5) || (157.5 <= direction[i * imageWidth + j] && direction[i * imageWidth + j] <= 180)) {
                q = src[i * imageWidth + j + 1];
                r = src[i * imageWidth + j - 1];
            } else if (22.5 <= direction[i * imageWidth + j] && direction[i * imageWidth + j] < 67.5) {
                q = src[(i + 1) * imageWidth + j - 1];
                r = src[(i - 1) * imageWidth + j + 1];
            } else if (67.5 <= direction[i * imageWidth + j] && direction[i * imageWidth + j] < 112.5) {
                q = src[(i + 1) * imageWidth + j];
                r = src[(i - 1) * imageWidth + j];
            } else if (112.5 <= direction[i * imageWidth + j] && direction[i * imageWidth + j] < 157.5) {
                q = src[(i - 1) * imageWidth + j - 1];
                r = src[(i + 1) * imageWidth + j + 1];
            }
            if ((src[i * imageWidth + j] >= q) && (src[i * imageWidth + j] >= r)) {
                dst[i * imageWidth + j] = src[i * imageWidth + j];
            } else {
                dst[i * imageWidth + j] = 0;
            }
        }
    }
}

int WEAK = 0;
int STRONG = 255;

void applyThreshold(BYTE *image, int imageWidth, int imageHeight) {
    float lowThresholdRatio = 0.05f;
    float highThresholdRatio = 0.09f;
    int maxPixel = 0;
    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            if (image[i * imageWidth + j] > maxPixel) {
                maxPixel = image[i * imageWidth + j];
            }
        }
    }
    int highThreshold = maxPixel * highThresholdRatio;
    int lowThreshold = highThreshold * lowThresholdRatio;
    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            if (image[i * imageWidth + j] >= highThreshold) {
                image[i * imageWidth + j] = STRONG;
            } else if (image[i * imageWidth + j] >= lowThreshold) {
                image[i * imageWidth + j] = WEAK;
            } else {
                image[i * imageWidth + j] = 0;
            }
        }
    }
}

void applyHysteresis(BYTE *image, int imageWidth, int imageHeight) {
    for (int i = 1; i < imageHeight - 1; i++) {
        for (int j = 1; j < imageWidth - 1; j++) {
            if (image[i * imageWidth+ j] == WEAK) {
                if ((image[(i + 1) * imageWidth + j - 1] == STRONG) || (image[(i + 1) * imageWidth + j] == STRONG) 
                || (image[(i + 1) * imageWidth + j + 1] == STRONG) || (image[i * imageWidth + j - 1] == STRONG) 
                || (image[i * imageWidth + j + 1] == STRONG) || (image[(i - 1) * imageWidth + j - 1] == STRONG) 
                || (image[(i - 1) * imageWidth + j] == STRONG) || (image[(i - 1) * imageWidth + j + 1] == STRONG)) {
                    image[i* imageWidth+ j] = STRONG;
                } else {
                    image[i* imageWidth+ j] = 0;
                }
            }
        }
    }
}

void printImageInfo(FIBITMAP *image) {
    std::cout << "BPP: " << bpp << " Width: " << width << " Height: " << height << " Pitch: " << pitch << " Red mask: " << FreeImage_GetRedMask(image) << std::endl;
}

void saveImage(const char *filepath, BYTE *dst, const char *stageName){
    FIBITMAP *image = FreeImage_ConvertFromRawBits(dst, width, height, pitch, bpp, 0, 0, 0);
    std::string outpath;
    outpath = std::string(filepath); 
    outpath.insert(outpath.length() - 4, stageName);
    FreeImage_Save(FIF_PNG, image, outpath.c_str());
}

int main(int argc, char* argv[]) {
    const char* filepath = "";

    if(argc < 3){
        std::cout << "Usage: ./main -f <image file>" << std::endl;
        exit(-1);
    }
    if(strcmp(argv[1], "-f") == 0){
        filepath = argv[2];
    }
    else{
        std::cout << "Usage: ./main -f <image file>" << std::endl;
        exit(-1);
    }
    //"images/valve.png"

    FIBITMAP *image = FreeImage_Load(FIF_PNG, filepath);

    image = FreeImage_ConvertToGreyscale(image);
    image = FreeImage_ConvertTo8Bits(image);

    bpp = FreeImage_GetBPP(image);
    width = FreeImage_GetWidth(image);
    height = FreeImage_GetHeight(image);
    pitch = FreeImage_GetPitch(image);

    printImageInfo(image);

    BYTE *src = new BYTE[width * height];
    BYTE *dst = new BYTE[width * height]();

    FreeImage_ConvertToRawBits(src, image, pitch, bpp, 0, 0, 0);

    applyBilateralFilter(src, dst, width, height, 8, 30.f, 20.f);
    saveImage(filepath, dst, "-1-bilateral");

    src = dst;
    dst = new BYTE[width * height]();
    float *direction = new float[width * height];

    applySobelFilter(src, dst, direction, width, height);
    saveImage(filepath, dst, "-2-sobel");

    src = dst;
    dst = new BYTE[width * height]();

    applyNonMaxSuppression(src, dst, direction, width, height);
    saveImage(filepath, dst, "-3-nonmax");

    applyThreshold(dst, width, height);
    saveImage(filepath, dst, "-4-thres");

    // src = dst;
    // dst = new BYTE[width * height]();

    // applyHysteresis(dst, width,height);
    // saveImage(filepath, dst, "-5-hyster");
    return 0;
}
