#pragma once

#include <image.h>
#include <huffman.h>
#include <fstream>
#include "huffman.h"

void histoGrey(Image &img, unsigned int **data){
    allocation_tableau(*data, unsigned int, 256);

    for(int i=0; i<img.totalSize; i++){
        (*data)[img.data[i]] ++;
    }
}


void histoColor(Image &img, unsigned int **data){
    allocation_tableau(*data, unsigned int, 256 * 3);

    for(int i=0; i<img.totalSize; i+=3){
        (*data)[img.data[i] * 3] ++;
        (*data)[img.data[i + 1] * 3 + 1] ++;
        (*data)[img.data[i + 2] * 3 + 2] ++;
    }
}

void storeHistoGrey(unsigned int *data, const char * name){
    std::ofstream file(name);
    
    if(!file.is_open()) return;


    for(int i=0; i<256; i++){
        file << i << " " << data[i] << "\n";
    }
}

void storeHistoColor(unsigned int *dataR, unsigned int *dataG, unsigned int *dataB, const char *name){
    std::ofstream file(name);
    
    if(!file.is_open()) return;


    for(int i=0; i<256; i++){
        file << i << " " << dataR[i] << " " << dataG[i] << " " << dataB[i] << "\n";
    }
}

void displayHistoColor(const char *name){
    // Recupérer de mes tps
    //system(
    //   "gnuplot -persist -e \"set term png; set output 'histoCouleur.png'; plot 'histocoul.dat' using 1:2 with lines title 'Red' linecolor 'red', \
    //       'histocoul.dat' using 1:3 with lines title 'Green' linecolor 'green', \
    //       'histocoul.dat' using 1:4 with lines title 'Blue' linecolor 'blue';\"");
    std::string command = "gnuplot -persist -e \"set term png; set output '";
    command += name;
    command += ".png'; plot '";
    command += name;
    command += "' using 1:2 with lines title 'Red' linecolor 'red', \
          '";
    command += name;
    command += "' using 1:3 with lines title 'Green' linecolor 'green', \
          '";
    command += name;
    command += "' using 1:4 with lines title 'Blue' linecolor 'blue';\"";
    std::cout << command << std::endl;
    system(command.c_str());
}

void histoColor(Image &img, const char *name){
    unsigned int *dataR, *dataG, *dataB;
    allocation_tableau(dataR, unsigned int, 256);
    allocation_tableau(dataG, unsigned int, 256);
    allocation_tableau(dataB, unsigned int, 256);

    for(int i = 0; i < img.totalSize; i += 3){
        dataR[img.data[i]]++;
        dataG[img.data[i + 1]]++;
        dataB[img.data[i + 2]]++;
    }

    // Store data
    std::ofstream file(name);
    if(!file.is_open()) return;
    for(int i=0; i<256; i++){
        file << i << " " << dataR[i] << " " << dataG[i] << " " << dataB[i] << "\n";
    }
    file.close();
    // Display data
    displayHistoColor(name);
}

double PSNR(Image &original, Image &compressed){
    double total = 0;  
    for(size_t i=0; i<original.totalSize; i++){
        double difference = original.data[i] - compressed.data[i];
        difference *= difference;
        total += difference;
    }

    double EQM = total / (double) original.totalSize;

    return 10. * log10((255. * 255.) / EQM);
}

/**
 * @brief Calculates the Boundary Recall between a segmented image and a reference segmentation.
 * Detects the contours of the superpixels via an RGB gradient, then compares with the ground truth
 * to assess the quality of adhesion to real contours. Returns a score between 0 and 1.
 */
double BoundaryRecall(Image &Compressed, Image &BoundaryRecall){
    int xR, yR, normeR;
    int xG, yG, normeG;
    int xB, yB, normeB;
    int seuil = 0;

    Image CompressedBoundary(Compressed.width, Compressed.height, false);
    int k = 0;

    // Contour of the image using RGB gradient detection
    for (int i = 0; i < Compressed.height; i++) {
        for (int j = 0; j < Compressed.width * 3; j += 3) {
            int R = Compressed.data[i * Compressed.width * 3 + j];
            int G = Compressed.data[i * Compressed.width * 3 + j + 1];
            int B = Compressed.data[i * Compressed.width * 3 + j + 2];

            // Handling of the edges
            if (i == Compressed.height - 1 || j / 3 == Compressed.width - 1) {
                CompressedBoundary.data[k] = 255; // Bords blancs
            } else {
                // Compute the gradient for each channel
                xR = Compressed.data[(i + 1) * Compressed.width * 3 + j] - R;
                yR = Compressed.data[i * Compressed.width * 3 + j + 3] - R;
                normeR = sqrt(xR * xR + yR * yR);

                xG = Compressed.data[(i + 1) * Compressed.width * 3 + j + 1] - G;
                yG = Compressed.data[i * Compressed.width * 3 + j + 3 + 1] - G;
                normeG = sqrt(xG * xG + yG * yG);

                xB = Compressed.data[(i + 1) * Compressed.width * 3 + j + 2] - B;
                yB = Compressed.data[i * Compressed.width * 3 + j + 3 + 2] - B;
                normeB = sqrt(xB * xB + yB * yB);

                // Edge detection if at least one channel exceeds the threshold
                CompressedBoundary.data[k] = (normeR > seuil || normeG > seuil || normeB > seuil) ? 255 : 0;
            }
            k++;
        }
    }

    //Save the image of the detected contours
    CompressedBoundary.write("output/CompressedBoundary.ppm");

    // Boundary Recall calculation
    double TP = 0, FP = 0, FN = 0;
    for (int i = 0; i < BoundaryRecall.height * BoundaryRecall.width; i++) {
        if (BoundaryRecall.data[i] == 255) { // Real contour detected
                                             // TODO : Modify to accept a tolerance because my current tests
                                             //  have been done with an image already thresholded
            bool found = false;

            // Eight-neighborhood search
            for (int di = -1; di <= 1; di++) {
                for (int dj = -1; dj <= 1; dj++) {
                    int ni = i / BoundaryRecall.width + di;
                    int nj = i % BoundaryRecall.width + dj;

                    if (ni >= 0 && ni < BoundaryRecall.height && nj >= 0 && nj < BoundaryRecall.width) {
                        int neighborIndex = ni * BoundaryRecall.width + nj;
                        if (CompressedBoundary.data[neighborIndex] == 255) {
                            found = true;
                            break;
                        }
                    }
                }
                if (found) break;
            }

            if (found) {
              TP++; // Correctly detected contour
              BoundaryRecall.data[i] = 0; // Mark as detected
            }
            else{
              FN++; // Contour not detected
            }
        }
        else if (CompressedBoundary.data[i] == 255) {
            FP++; // False positive (detected contour but absent from the ground truth) but it doesn't matter for the recall
        }
    }

	// Save the boundary recall map (useful for visualizing errors)
    BoundaryRecall.write("output/BoundaryRecall.ppm");

    return TP / (TP + FN); // Return the Boundary Recall
}


// void compressFile(bool compress, char *srcPath, char *dstPath)
// {
//     // FILE *src, *dst, *frq;

//     // if ((src=fopen(srcPath, "rb"))==NULL) perror("fopen");
//     // if ((dst=fopen(dstPath, "wb"))==NULL) perror("fopen");

//     // frq=NULL;

//     // if ((arbre_d=(struct arbre_data *)malloc(512*sizeof(struct arbre_data)))==NULL) perror("malloc");

//     // if ((arbre=(struct arbre *)malloc(512*sizeof(struct arbre)))==NULL) {
//     //     free(arbre_d);
//     //     perror("malloc");
//     // }

//     // if (compress)
//     //     huffman_compacter(src, dst, frq);
//     // else
//     //     huffman_decompacter(src, dst, frq);


//     // free(arbre_d);
//     // free(arbre);

//     // fclose(src);
//     // fclose(dst);
//     // if (frq!=NULL)
//     //     fclose(frq);


//     if(compress){
//         compressFile(srcPath, dstPath);
//     } else {
//         decompressFile(srcPath, dstPath);
//     }
// }

uintmax_t getFileSize(char *path){
    return filesystem::file_size(path);
}


struct Pattern {
    float *data;
    int sideLength;

    Pattern(float *values, int sideLength){
        data = values;
        this->sideLength = sideLength;
    }

    int apply(int idx, Image &image){
        int offset = sideLength/2;

        int topLeftIdx = idx - image.width * offset - offset;

        float total = 0;

        for(int i=0; i<sideLength; i++){
            for(int j=0; j<sideLength; j++){
                int correctIdx = topLeftIdx + i * image.width + j;
                if(correctIdx >= 0 && correctIdx < image.nbPixel){
                    total += data[i * sideLength + j] * (float) image.data[correctIdx];
                }
            }
        }

        return total;
    }
};

void processGradient(Image &imageIn, Image &imageOut){
    imageOut = Image(imageIn.width, imageIn.height, false);
    

    float sobelDataX[9] = {-1, 0, 1,
                           -2, 0, 2,
                           -1, 0, 1};
    float sobelDataY[9] = {-1,-2,-1,
                            0, 0, 0,
                            1, 2, 1};
    
    Pattern Gx(sobelDataX, 3);
    Pattern Gy(sobelDataY, 3);

    for(int i=1; i<imageIn.height; i++){
        for(int j=1; j<imageIn.width; j++){
            float x = Gx.apply(i*imageIn.width + j, imageIn);
            float y = Gy.apply(i*imageIn.width + j, imageIn);
            float g = sqrt(x*x + y*y);
            imageOut[i*imageIn.width + j] = g + 127;
        }
    }
}



