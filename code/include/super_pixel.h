#include <cmath>
#include <queue>
#include <vector>
#include <algorithm>

#include <Image.h>

struct Rgb {
    OCTET r,g,b;
    Rgb(OCTET r, OCTET g, OCTET b){
        this->r = r;
        this->g = g;
        this->b = b;
    }
};

struct Xyz {
    double x,y,z;
};

struct Lab {
    double l, a, b;
    Lab(){
        l = a = b = 0;
    }

    void add(Lab &other){
        l += other.l;
        a += other.a;
        b += other.b;
    }

    void div(double value){
        l /= value;
        a /= value;
        b /= value;
    }
};

Xyz rgbToXyz(Rgb rgb){
    double R = rgb.r / 255.0;
    double G = rgb.g / 255.0;
    double B = rgb.b / 255.0;

    // Application de la transformation gamma inverse (standard de l'espace sRGB)
    R = (R > 0.04045) ? std::pow((R + 0.055) / 1.055, 2.4) : R / 12.92;
    G = (G > 0.04045) ? std::pow((G + 0.055) / 1.055, 2.4) : G / 12.92;
    B = (B > 0.04045) ? std::pow((B + 0.055) / 1.055, 2.4) : B / 12.92;


    // Transformation RGB → XYZ
    double X = R * 0.4124564 + G * 0.3575761 + B * 0.1804375;
    double Y = R * 0.2126729 + G * 0.7151522 + B * 0.0721750;
    double Z = R * 0.0193339 + G * 0.1191920 + B * 0.9503041;

    Xyz res;
    res.x = X;
    res.y = Y;
    res.z = Z;

    return res;
}

Lab rgbToLab(Rgb rgb) {
    Xyz xyz = rgbToXyz(rgb);

    double epsilon = 0.008856; //actual CIE standard
    double kappa = 903.3; //actual CIE standard


    // Xr = 0.950456 
    // Yr = 1.0 
    // Zr = 1.088754 

    xyz.x /= 0.95047;
    xyz.y /= 1.00000;
    xyz.z /= 1.08883;

    // Fonction de transformation f(t)
    auto f = [](double t) -> double {
        return (t > 0.008856) ? std::pow(t, 1.0 / 3.0) : (t * 903.3 + 16.0) / 116.0;
    };

    double X = f(xyz.x);
    double Y = f(xyz.y);
    double Z = f(xyz.z);

    // Calcul du CIELAB
    Lab lab;
    lab.l = 116.0 * Y - 16.0; // (Y > epsilon) ? (116.0 * std::pow(Y, 1.0 / 3.0) - 16.0) : (kappa * Y);
    lab.a = 500.0 * (X - Y);
    lab.b = 200.0 * (Y - Z);

    return lab;
}


struct Pixel {
    int x, y, superpixel_id;
    float distance;
    Lab lab;

    Pixel(){
        x = y = distance = 0;
        superpixel_id = -1;
        lab = Lab();
    }

    // Comparateur pour créer un min-heap (plus petite distance en premier)
    bool operator>(const Pixel& other) const {
        return distance > other.distance;
    }
};


struct SuperPixel{
    std::vector<Pixel> pixels;

    SuperPixel(){
        pixels = std::vector<Pixel>();
    }

    void getAveragePixel(Pixel &res){
        res = Pixel();
        for(Pixel p : pixels){
            res.lab.add(p.lab);
            res.x += p.x;
            res.y += p.y;
        }
        res.lab.div(pixels.size());
        res.x /= pixels.size();
        res.y /= pixels.size();
    }
};

double distance(Pixel &A, Pixel &B,
                double s, double m)
{
    double posDist = pow(A.x - B.x, 2) + pow(A.y - B.y, 2);
    double labDist = pow(A.lab.l - B.lab.l, 2) + pow(A.lab.a - B.lab.a, 2) + pow(B.lab.b - B.lab.b, 2);
    return sqrt(posDist / s + labDist / m);
}

// points clés de SNIC: n'utilise pas kmean, pas besoin de plusieurs itérations et meilleure connectivité dés le début
void SNIC(Image &imageIn, Image &imageOut){
    int k = 10; // Nb de centroid
    
    double s = sqrt((double) imageIn.nbPixel / (double) k ); // for squared super pixels
    double m = 10.0; // compactness factor


    std::vector<SuperPixel> superPixels;
    superPixels.resize(k);

    std::vector<Pixel> imagePixels;
    for(int i=0; i<imageIn.nbPixel; i++){
        int idx = i*3;
        Rgb rgb(imageIn[idx], imageIn[idx + 1], imageIn[idx + 2]);

        Lab lab = rgbToLab(rgb);
        int x, y = 0;
        x = i % imageIn.width;
        y = i - x;

        Pixel pix;
        pix.x = x;
        pix.y = y;
        pix.lab = lab;
        // pix.distance = lab.l;
        imagePixels.push_back(pix);
    }


    // init priority queue
    std::priority_queue<Pixel, std::vector<Pixel>, std::greater<Pixel>> Q;
    int iterationStep = imageIn.nbPixel / k;
    for(int i=0; i<k; i++){
        Q.push(imagePixels[(i + 1) * iterationStep - (iterationStep/2)]);
    }


    while (!Q.empty()) {
        Pixel p = Q.top();
        if(p.superpixel_id == -1){

        }
        Q.pop();
        std::cout << "Pixel (" << p.x << ", " << p.y << ") avec distance " << p.distance << "\n";
    }
}