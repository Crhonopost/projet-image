#include <queue>

#include <image.h>
#include <format.h>
#include <fstream>
#include <util.h>
#include <huffman.h>

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


void storeSNIC(int width, int height, const std::vector<Rgb> &palette, const std::vector<int> &superpixelIds, std::vector<unsigned char> &buffer) {
    buffer.clear();

    int paletteSize = palette.size();

    buffer.insert(buffer.end(), reinterpret_cast<unsigned char*>(&width), reinterpret_cast<unsigned char*>(&width) + sizeof(int));
    buffer.insert(buffer.end(), reinterpret_cast<unsigned char*>(&height), reinterpret_cast<unsigned char*>(&height) + sizeof(int));
    buffer.insert(buffer.end(), reinterpret_cast<unsigned char*>(&paletteSize), reinterpret_cast<unsigned char*>(&paletteSize) + sizeof(int));

    buffer.insert(buffer.end(), reinterpret_cast<const unsigned char*>(palette.data()), reinterpret_cast<const unsigned char*>(palette.data()) + paletteSize * sizeof(Rgb));

    std::vector<short int> differentialIndices;
    int previousIdxVal = 0;
    int maxDiff = 0;

    for(int idx : superpixelIds){
        int diff = idx - previousIdxVal;
        maxDiff = std::max(maxDiff, std::abs(diff));

        differentialIndices.push_back(static_cast<short int>(diff));
        previousIdxVal = idx;
    }

    std::cout << "Différence max entre 2 pixels : " << maxDiff << "\n";

    buffer.insert(buffer.end(), reinterpret_cast<const unsigned char*>(differentialIndices.data()), reinterpret_cast<const unsigned char*>(differentialIndices.data()) + differentialIndices.size() * sizeof(short int));
}


void readSNIC(const std::vector<unsigned char> &buffer, Image &image) {
    if (buffer.empty()) {
        std::cerr << "Error: Empty buffer." << std::endl;
        return;
    }

    size_t pos = 0;

    // Lecture des métadonnées
    int width, height, paletteSize;
    std::memcpy(&width, &buffer[pos], sizeof(int)); pos += sizeof(int);
    std::memcpy(&height, &buffer[pos], sizeof(int)); pos += sizeof(int);
    std::memcpy(&paletteSize, &buffer[pos], sizeof(int)); pos += sizeof(int);

    int totalSize = width * height;

    // Lecture de la palette
    std::vector<Rgb> palette(paletteSize);
    std::memcpy(palette.data(), &buffer[pos], paletteSize * sizeof(Rgb));
    pos += paletteSize * sizeof(Rgb);

    // Lecture des indices différentiels
    std::vector<short int> differentialIds(totalSize);
    std::memcpy(differentialIds.data(), &buffer[pos], totalSize * sizeof(short int));
    pos += totalSize * sizeof(short int);

    // Reconstruction des indices absolus
    std::vector<int> ids(totalSize);
    ids[0] = differentialIds[0];
    for (int i = 1; i < totalSize; i++) {
        ids[i] = ids[i - 1] + differentialIds[i];
    }

    // Création de l'image
    image = Image(width, height, true);
    for (int i = 0; i < image.nbPixel; i++) {
        Rgb pixelColor = palette[ids[i]];
        image[i * 3] = pixelColor.r;
        image[i * 3 + 1] = pixelColor.g;
        image[i * 3 + 2] = pixelColor.b;
    }
}

// ref: https://openaccess.thecvf.com/content_cvpr_2017/papers/Achanta_Superpixels_and_Polygons_CVPR_2017_paper.pdf
// points clés de SNIC: n'utilise pas kmean, pas besoin de plusieurs itérations et meilleure connectivité dés le début
// imageIn : image d'entrée
// imageOut : image de sortie
// k : nombre de superpixels
// m : Facteur de compacité
void SNIC(Image &imageIn, Image &imageOut, int k = 5000, double m = 10.0) {
    double s = sqrt((double)imageIn.nbPixel / (double)k); // Taille d'un superpixel

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

    std::priority_queue<Pixel*, std::vector<Pixel*>, ComparePixels> Q;
    

    int gridSize = sqrt(k);
    double stepX = (double)imageIn.width / gridSize;
    double stepY = (double)imageIn.height / gridSize;

    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            int x = round((j + 0.5) * stepX);
            int y = round((i + 0.5) * stepY);

            x = std::min(x, imageIn.width - 1);
            y = std::min(y, imageIn.height - 1);

            int index = y * imageIn.width + x;
            
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
    }
    
    
    imageOut = Image(imageIn.width, imageIn.height, true);
    std::vector<int> ids;
    // Mise à jour des centroïdes
    for(int i=0; i<imageOut.nbPixel; i++){
        int superPixelIdx = (*imagePixels[i]).superpixel_id;
        ids.push_back(superPixelIdx);
        Rgb &color = colors[superPixelIdx];
        imageOut[i*3] = color.r;
        imageOut[i*3 + 1] = color.g;
        imageOut[i*3 + 2] = color.b;
    }

    // imageOut.write("output/SP_without_compression.ppm");

    std::vector<unsigned char> data;
    storeSNIC(imageIn.width, imageIn.height, colors, ids, data);// "output/taupe_SNIC.snic");
    compressFile("output/taupe.snic", data);

    for(Pixel * pixel: imagePixels){
        free(pixel);
    }
}
