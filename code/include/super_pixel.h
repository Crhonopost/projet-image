#include <cmath>
#include <queue>
#include <vector>
#include <algorithm>

#include <image.h>

struct Rgb {
    OCTET r,g,b;
    Rgb(): r(0), g(0), b(0){}
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

Rgb labToRgb(Lab lab) {
    // Inverse transformation f^-1(t)
    auto f_inv = [](double t) -> double {
        double t3 = t * t * t;
        return (t3 > 0.008856) ? t3 : (116.0 * t - 16.0) / 903.3;
    };
    
    // Calcul des valeurs XYZ
    double Y = (lab.l + 16.0) / 116.0;
    double X = lab.a / 500.0 + Y;
    double Z = Y - lab.b / 200.0;
    
    X = f_inv(X) * 0.95047;
    Y = f_inv(Y) * 1.00000;
    Z = f_inv(Z) * 1.08883;
    
    // Transformation XYZ → RGB
    double R = X *  3.2404542 + Y * -1.5371385 + Z * -0.4985314;
    double G = X * -0.9692660 + Y *  1.8760108 + Z *  0.0415560;
    double B = X *  0.0556434 + Y * -0.2040259 + Z *  1.0572252;
    
    // Application de la transformation gamma
    auto gamma_correction = [](double c) -> double {
        return (c > 0.0031308) ? (1.055 * std::pow(c, 1.0 / 2.4) - 0.055) : 12.92 * c;
    };
    
    R = gamma_correction(R);
    G = gamma_correction(G);
    B = gamma_correction(B);
    
    // Clamp et conversion en entiers [0, 255]
    auto clamp = [](double v) -> int {
        return static_cast<int>(std::round(std::max(0.0, std::min(1.0, v)) * 255.0));
    };
    
    Rgb rgb;
    rgb.r = clamp(R);
    rgb.g = clamp(G);
    rgb.b = clamp(B);
    
    return rgb;
}

struct Pixel {
    int superpixel_id, potential_sp_id;
    Lab lab;
    double distance, x, y;

    Pixel() : x(0), y(0), superpixel_id(-1), potential_sp_id(-1), lab(), distance(MAXFLOAT) {}

    static double computeDistance(const Pixel& A, const Pixel& B, double s, double m) {
        double posDist = pow(A.x - B.x, 2) + pow(A.y - B.y, 2);
        double labDist = pow(A.lab.l - B.lab.l, 2) + pow(A.lab.a - B.lab.a, 2) + pow(A.lab.b - B.lab.b, 2);
        return sqrt(posDist / s + labDist / m);
    }
};

struct ComparePixels {
    bool operator()(const Pixel* a, const Pixel* b) const {
        return a->distance > b->distance; // Min-heap : plus petite distance en premier
    }
};

void getNeighbors(Pixel p, Image &image, std::vector<Pixel*> &pixels, std::vector<Pixel*> &res){
    res.resize(0);
    if(p.x > 0){
        res.push_back(pixels[p.y * image.width + p.x - 1]);
    }
    if((p.y * image.width + p.x + 1) < image.nbPixel){
        res.push_back(pixels[p.y * image.width + p.x + 1]);
    }

    if(p.y > 0){
        res.push_back(pixels[(p.y - 1) * image.width + p.x]);
    }
    if((p.y + 1) * image.width + p.x < image.nbPixel){
        res.push_back(pixels[(p.y + 1) * image.width + p.x]);
    }
}


struct SuperPixel{
    std::vector<Pixel> pixels;

    SuperPixel(){
        pixels = std::vector<Pixel>();
    }

    Pixel getAveragePixel(){
        Pixel res = Pixel();
        res.distance = 0;
        for(Pixel p : pixels){
            res.lab.add(p.lab);
            res.x += p.x;
            res.y += p.y;
        }
        res.lab.div(pixels.size());
        res.x /= pixels.size();
        res.y /= pixels.size();
        
        return res;
    }
};


// ref: https://openaccess.thecvf.com/content_cvpr_2017/papers/Achanta_Superpixels_and_Polygons_CVPR_2017_paper.pdf
// points clés de SNIC: n'utilise pas kmean, pas besoin de plusieurs itérations et meilleure connectivité dés le début
void SNIC(Image &imageIn, Image &imageOut) {
    int k = 50000; // Nombre de superpixels
    double s = sqrt((double)imageIn.nbPixel / (double)k); // Taille d'un superpixel
    double m = 10.0; // Facteur de compacité

    std::vector<SuperPixel> superPixels;
    std::vector<Pixel*> imagePixels;

    for (int i = 0; i < imageIn.height; i++) {
        for (int j = 0; j < imageIn.width; j++) {
            int idx = (i * imageIn.width + j) * 3;

            Rgb rgb(imageIn[idx], imageIn[idx + 1], imageIn[idx + 2]);
            Lab lab = rgbToLab(rgb);

            Pixel *pix = new Pixel();
            (*pix).y = i;
            (*pix).x = j;
            (*pix).lab = lab;
            (*pix).superpixel_id = -1;
            (*pix).distance = std::numeric_limits<double>::max();
            imagePixels.push_back(pix);
        }
    }

    // Initialisation de la priority queue
    std::priority_queue<Pixel*, std::vector<Pixel*>, ComparePixels> Q;
    

    int gridSize = sqrt(k); // Nombre de superpixels par ligne et colonne
    double stepX = (double)imageIn.width / gridSize;
    double stepY = (double)imageIn.height / gridSize;

    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            int x = round((j + 0.5) * stepX); // Centré dans la cellule
            int y = round((i + 0.5) * stepY);

            x = std::min(x, imageIn.width - 1);  // Évite de sortir des bords
            y = std::min(y, imageIn.height - 1);

            int index = y * imageIn.width + x; // Convertit (x, y) en index de `imagePixels`
            
            Pixel &centroid = *imagePixels[index];
            centroid.superpixel_id = superPixels.size();
            centroid.potential_sp_id = superPixels.size();
            centroid.distance = 0;
            Q.push(&centroid);

            SuperPixel sp;
            sp.pixels.push_back(centroid);
            superPixels.push_back(sp);
        }
    }


    // int iterationStep = imageIn.nbPixel / k;
    
    // for (int i = 0; i < k; i++) {
    //     for (int j = 0; j < k; j++) {
    //         Pixel& centroid = imagePixels[i * imageIn.width iterationStep - (iterationStep / 2)];
    //         centroid.superpixel_id = i;
    //         centroid.distance = 0;
    //         Q.push(centroid);
    //         superPixels[i].pixels.push_back(centroid);
    //     }
    // }

    int loop = 0;

    while (!Q.empty()) {
        std::cout << "\033[2J\033[H";
 
        std::cout << "i - " << loop << "\n"
        << "Q size : " << Q.size() << "\n";

        loop ++;
        Pixel p = *Q.top();
        Q.pop();

        p.superpixel_id = p.potential_sp_id;

        if ((*imagePixels[p.y * imageIn.width + p.x]).superpixel_id == -1) {
            (*imagePixels[p.y * imageIn.width + p.x]).superpixel_id = p.superpixel_id;
            superPixels[p.superpixel_id].pixels.push_back(p);
        }

        
        std::vector<Pixel*> neighbors;
        getNeighbors(p, imageIn, imagePixels, neighbors);
        for (Pixel *neighborPtr : neighbors) {
            Pixel &neighbor = *neighborPtr;
            double dist = Pixel::computeDistance(neighbor, superPixels[p.superpixel_id].getAveragePixel(), s, m);
            
            if (neighbor.potential_sp_id == -1) {
                Q.push(&neighbor);
            }

            if(dist < neighbor.distance){
                neighbor.distance = dist;
                neighbor.potential_sp_id = p.superpixel_id;
            }

        }
    }

    
    std::vector<Rgb> colors;
    for(int i=0; i<superPixels.size(); i++){
        Lab lab = superPixels[i].getAveragePixel().lab;
        Rgb color = labToRgb(lab);
        colors.push_back(color);

        std::cout << (int)color.r << " / " << (int)color.g << " / " << (int)color.b << "\n";
    }
    
    
    imageOut = Image(imageIn.width, imageIn.height, true);
    // Mise à jour des centroïdes
    for(int i=0; i<imageOut.nbPixel; i++){
        int superPixelIdx = (*imagePixels[i]).superpixel_id;
        Rgb &color = colors[superPixelIdx];
        imageOut[i*3] = color.r;
        imageOut[i*3 + 1] = color.g;
        imageOut[i*3 + 2] = color.b;
    }

    imageOut.write("output/res.ppm");



    for(Pixel * pixel: imagePixels){
        free(pixel);
    }
}
