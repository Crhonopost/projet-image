#pragma once
#include <image_ppm.h>
#include <vector>
#include "Vec3.h"
#include <regex>
#include <cstring>
#include <stddef.h>

struct Image {
    int width, height;
    OCTET *data;
    size_t totalSize;
    size_t nbPixel;
    bool isColor;

    OCTET& operator[](size_t i) {return data[i];}
    const OCTET& operator[](size_t i) const {return data[i];}
    
    Image(){};

    Image(Image &cpy){
        this->height = cpy.height;
        this->width = cpy.width;
        this->totalSize = cpy.totalSize;
        this->nbPixel = cpy.nbPixel;
        this->isColor = cpy.isColor;

        allocation_tableau(this->data, OCTET, this->totalSize);

        std::memcpy(this->data, cpy.data, this->totalSize);
    }
    Image(int width, int height, bool isColor = false){
        this->height = height;
        this->width = width;
        this->isColor = isColor;
        nbPixel = width * height;
        totalSize = isColor ? nbPixel * 3 : nbPixel;
    
        allocation_tableau(this->data, OCTET, this->totalSize);
    }

    // ~Image() {
    //     if (data) {
    //         free(data);
    //         data = nullptr;
    //     }
    // }

    void read(char *imagePath){
        std::regex ppm_regex(".*\\.ppm$", std::regex::icase);
        bool isPPM = std::regex_match(imagePath, ppm_regex);
        if(isPPM){
            isColor = true;
            lire_nb_lignes_colonnes_image_ppm(imagePath, &height, &width);

            nbPixel = width * height;
            totalSize = nbPixel * 3;
            allocation_tableau(data, OCTET, totalSize);
            lire_image_ppm(imagePath, data, nbPixel);
        } else {
            isColor = false;
            lire_nb_lignes_colonnes_image_pgm(imagePath, &height, &width);

            nbPixel = width * height;
            totalSize = nbPixel;
            allocation_tableau(data, OCTET, totalSize);
            lire_image_pgm(imagePath, data, nbPixel);
        }
    }
    void write(const char *imagePath){
        if(isColor){
            ecrire_image_ppm(imagePath, data,  height, width);
        } else {
            ecrire_image_pgm(imagePath, data,  height, width);
        }
    }
    void split(Image &r, Image &g, Image &b){
        allocation_tableau(r.data, OCTET, nbPixel);
        allocation_tableau(g.data, OCTET, nbPixel);
        allocation_tableau(b.data, OCTET, nbPixel);
        
        r.width = g.width = b.width = width;
        r.height = g.height = b.height = height;
        r.nbPixel = g.nbPixel = b.nbPixel = nbPixel;
        r.totalSize = g.totalSize = b.totalSize = nbPixel;
        r.isColor = g.isColor = b.isColor = false;

        for(size_t i=0; i<nbPixel; i++){
            size_t totalI = i*3;
            r.data[i] = data[totalI];
            g.data[i] = data[totalI + 1];
            b.data[i] = data[totalI + 2];
        }
    }
    void fuse(Image &r, Image &g, Image &b){
        height = r.height;
        width = r.width;
        nbPixel = r.nbPixel;
        totalSize = nbPixel * 3;
        isColor = true;

        allocation_tableau(data, OCTET, totalSize);

        for(int i=0; i<nbPixel; i++){
            int colorIdx = i*3;
            data[colorIdx] = r.data[i];
            data[colorIdx + 1] = g.data[i];
            data[colorIdx + 2] = b.data[i];
        }
    }
    void ppmToVec3(std::vector<Vec3> &res){
        for(size_t i=0; i<nbPixel; i++){
            size_t idx = i*3;
            Vec3 color(data[idx], data[idx+1], data[idx+2]);
            res.push_back(color);
        }
    }
    void vec3ToPPM(std::vector<Vec3> &dataVec){
        for(size_t i=0; i<nbPixel; i++){
            Vec3 color = dataVec[i];
            data[i*3] = color[0];
            data[i*3 + 1] = color[1];
            data[i*3 + 2] = color[2];
        }
    }
    void toPGM(Image &imgOut){
        imgOut = Image(width, height, false);

        for(int i=0; i<nbPixel; i++){
            float grey = 0.299f * (float) data[i*3] + 
                         0.587f * (float) data[i*3 + 1] + 
                         0.114f * (float) data[i*3 + 2];
            imgOut.data[i] = std::roundf(grey);
        }

    }
    void toPPM(Image &imgOut){
        imgOut = Image(width, height, true);

        for(int i=0; i<nbPixel; i++){
            imgOut.data[i*3] = data[i];
            imgOut.data[i*3 + 1] = data[i];
            imgOut.data[i*3 + 2] = data[i];
        }
    }

};
