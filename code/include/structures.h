#pragma once

#include <format.h>

struct coords {
    int x,y;
    coords(): x(0), y(0){};
    coords(coords &other){
        this->x = other.x;
        this->y = other.y;
    }
    coords(int x, int y): x(x), y(y){}
};

struct Pixel {
    int superpixel_id, temp_id;
    Lab lab;
    double distance, x, y;

    Pixel() : x(0), y(0), superpixel_id(-1), temp_id(-1), lab(), distance(MAXFLOAT) {}

    static double computeDistance(const Pixel& A, const Pixel& B, double s, double m) {
        double posDist = pow(A.x - B.x, 2) + pow(A.y - B.y, 2);
        double labDist = pow(A.lab.l - B.lab.l, 2) + pow(A.lab.a - B.lab.a, 2) + pow(A.lab.b - B.lab.b, 2);
        return sqrt(posDist / s + labDist / m);
    }

    static double computeDistanceSLIC(const Pixel& A, const Pixel& B, double s, double m) {
        double posDist = pow(A.x - B.x, 2) + pow(A.y - B.y, 2);
        double labDist = pow(A.lab.l - B.lab.l, 2) + pow(A.lab.a - B.lab.a, 2) + pow(A.lab.b - B.lab.b, 2);
        return sqrt(posDist + (m/s) * labDist);
    }


    static double computeSpaceDistance(const Pixel& A, const Pixel& B){
        return sqrt(pow(A.x - B.x, 2) + pow(A.y - B.y, 2));
    }
};

struct ComparePixels {
    bool operator()(const Pixel* a, const Pixel* b) const {
        return a->distance > b->distance; // Min-heap : plus petite distance en premier
    }
};


struct SuperPixel{
    std::vector<Pixel*> pixels;

    SuperPixel(){
        pixels = std::vector<Pixel*>();
    }

    Pixel getAveragePixel(){
        Pixel res = Pixel();
        res.distance = 0;
        for(Pixel *p : pixels){
            res.lab.add(p->lab);
            res.x += p->x;
            res.y += p->y;
        }
        res.lab.div(pixels.size());
        res.x /= pixels.size();
        res.y /= pixels.size();

        return res;
    }
};



// Structure pour les waterpixels
struct PixelNDG {
    int x, y;
    int value; // valeur extraite de l'image greg
};

struct ComparePixelNDG {
    bool operator()(const PixelNDG *a, const PixelNDG* b) const {
        return a->value > b->value; // Min-heap : plus petite valeur en premier
    }
};

struct minimaLocal
{
    int index;
    int value;
    int aire;
    int marqueur;
};
