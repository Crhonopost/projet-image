#pragma once

#include <queue>

#include <image.h>
#include <format.h>
#include <fstream>
#include <util.h>
#include <structures.h>
#include <huffman.h>
#include <logger.h>


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

void conversionImageVector(Image &imageIn, std::vector<Pixel*> &imagePixels){
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
}

void gridInitSLIC(Image &imageIn, std::vector<SuperPixel> &superPixels, std::vector<Pixel*> &imagePixels, int k){
    conversionImageVector(imageIn, imagePixels);

    int gridSize = sqrt(k);
    double stepX = (double)imageIn.width / gridSize;
    double stepY = (double)imageIn.height / gridSize;

    
    Image inGrey, imageGradient;
    imageIn.toPGM(inGrey);
    processGradient(inGrey, imageGradient);

    std::vector<coords> bestCentroidsCoords(gridSize * gridSize);
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            int x = round((j + 0.5) * stepX);
            int y = round((i + 0.5) * stepY);

            x = std::min(x, imageIn.width - 1);
            y = std::min(y, imageIn.height - 1);

            coords c(x,y);
            bestCentroidsCoords[i * gridSize + j] = c;
        }
    }

    for(coords &centroid : bestCentroidsCoords){
        int minGradient = abs(imageGradient[centroid.y * imageGradient.width + centroid.x] - 127);
        coords bestCoords(centroid);

        for(int x=-1; x<1; x++){
            for(int y=-1; y<1; y++){
                int neighborGradient = abs(imageGradient[(centroid.y + y) * imageGradient.width + centroid.x + x] - 127);
                if(neighborGradient < minGradient){
                    minGradient = neighborGradient;
                    bestCoords.x = centroid.x + x;
                    bestCoords.y = centroid.y + y;
                }
            }
        }

        centroid.x = bestCoords.x;
        centroid.y = bestCoords.y;
    }

    for(coords &centroidCoords : bestCentroidsCoords){
        int index = centroidCoords.y * imageIn.width + centroidCoords.x;
        
        Pixel *centroid = imagePixels[index];
        centroid->superpixel_id = superPixels.size();
        centroid->temp_id = superPixels.size();
        centroid->distance = 0;
    
        SuperPixel sp;
        sp.pixels.push_back(centroid);
        superPixels.push_back(sp);
    }


    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            int x = round((j + 0.5) * stepX);
            int y = round((i + 0.5) * stepY);

            x = std::min(x, imageIn.width - 1);
            y = std::min(y, imageIn.height - 1);

            int index = y * imageIn.width + x;
            
            Pixel *centroid = imagePixels[index];
            centroid->superpixel_id = superPixels.size();
            centroid->temp_id = superPixels.size();
            centroid->distance = 0;

            SuperPixel sp;
            sp.pixels.push_back(centroid);
            superPixels.push_back(sp);
        }
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

    gridInitSLIC(imageIn, superPixels, imagePixels, k);

    std::priority_queue<Pixel*, std::vector<Pixel*>, ComparePixels> Q;
    for(SuperPixel &p: superPixels){
        Q.push(p.pixels[0]);
    }

    int loop = 0;

    while (!Q.empty()) {
        std::cout << "\033[2J\033[H";
 
        std::cout << "i - " << loop << "\n"
        << "Q size : " << Q.size() << "\n";

        loop ++;
        Pixel *p = Q.top();
        Q.pop();

        p->superpixel_id = p->temp_id;

        if ((*imagePixels[p->y * imageIn.width + p->x]).superpixel_id == -1) {
            (*imagePixels[p->y * imageIn.width + p->x]).superpixel_id = p->superpixel_id;
            superPixels[p->superpixel_id].pixels.push_back(p);
        }

        
        std::vector<Pixel*> neighbors;
        getNeighbors(*p, imageIn, imagePixels, neighbors);
        for (Pixel *neighborPtr : neighbors) {
            Pixel &neighbor = *neighborPtr;
            double dist = Pixel::computeDistance(neighbor, superPixels[p->superpixel_id].getAveragePixel(), s, m);
            
            if (neighbor.temp_id == -1) {
                Q.push(&neighbor);
            }

            if(dist < neighbor.distance){
                neighbor.distance = dist;
                neighbor.temp_id = p->superpixel_id;
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

// https://vision.gel.ulaval.ca/~jflalonde/cours/4105/h17/tps/results/projet/111063028/index.html
// https://darshita1405.medium.com/superpixels-and-slic-6b2d8a6e4f08

void SLIC(Image &imageIn, Image &imageOut, int k = 5000, double m = 10.0) {
    ScopedLogger logger("SLIC");
    double s = sqrt((double)imageIn.nbPixel / (double)k); // Taille d'un superpixel

    std::vector<SuperPixel> superPixels;
    std::vector<Pixel*> imagePixels;

    gridInitSLIC(imageIn, superPixels, imagePixels, k);

    for(Pixel *pix: imagePixels){
        if(pix->superpixel_id == -1){
            pix->superpixel_id = 0;
            pix->temp_id = superPixels[0].pixels.size();
            superPixels[0].pixels.push_back(pix);
        } 
    }

    for(int iteration=0; iteration < 10; iteration++){
        ScopedLogger logger("SLIC iteration");
        auto startIte = std::chrono::high_resolution_clock::now();
        for(int spIdx=0; spIdx<superPixels.size(); spIdx++){
            ScopedLogger logger("SLIC super pixels");
            SuperPixel &sp = superPixels[spIdx];
            Pixel average = sp.getAveragePixel();
            int minX = max(average.x - s, 0.);
            int minY = max(average.y - s, 0.);

            for(int i=0; 
                i < 2*s && 
                i+minY < imageIn.height; 
                i++){
                for(int j=0; 
                    j < 2*s && 
                    j+minX < imageIn.width; 
                    j++){
                    ScopedLogger logger("SLIC pixels in super pixels");
                    Pixel *pix = imagePixels[(minY + i) * imageIn.width + j + minX];
                    double distance = Pixel::computeDistanceSLIC(*pix, average, s, m);
                    if(distance < pix->distance){
                        SuperPixel &concernedSP = superPixels[pix->superpixel_id];
                        concernedSP.pixels.erase(concernedSP.pixels.begin() + pix->temp_id);

                        for(int i=pix->temp_id; i<concernedSP.pixels.size(); i++){
                            concernedSP.pixels[i]->temp_id -= 1;
                        }

                        pix->temp_id = sp.pixels.size();
                        pix->superpixel_id = spIdx;
                        pix->distance = distance;
                        sp.pixels.push_back(pix);
                    }
                }
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

    GlobalLogger::getInstance().printAverages();
}

// Reg sert a savoir si on veut plus de régularité ou plus suivre les contours
// Plus reg est grand, plus on aura de régularité (grille carré)
void Waterpixel(Image& imageIn, Image& imageOut, int k, float percentageRho=0.2f, double reg=50., bool debugMode = false){
    int nH = imageIn.height;
    int nW = imageIn.width;
    int nTaille = nH * nW;

    // Conversion en structure pixels
    std::vector<Pixel*> imagePixels;
    conversionImageVector(imageIn, imagePixels);

    // ---------------------
    // Calcul du gradient
    // ---------------------
    Image ImgNormeGradient;
    ImgNormeGradient = Image(nW, nH, false);
    std::vector<PixelNDG*> ImgNormeGradientVector;
    processGradientLab(imagePixels, ImgNormeGradientVector, ImgNormeGradient);

    // Debug
    if (debugMode){
        std::string normeGradient = std::string("WATERPIXEL-DEBUG_ImgNormeGradient.pgm");
        ImgNormeGradient.write(normeGradient.c_str());
    }

    struct centersSquare{
        int x, y;
    };

    // ---------------------
    // Partie grille génération des centres
    // ---------------------
    std::vector<centersSquare> centers;
    // Paramètres
    float sigma = std::sqrt((nW * nH) / static_cast<float>(k));
    float sigma_y = sigma * std::sqrt(3.f) / 2.f; // Sqrt pour avoir un petit décalage en haut et en bas
    // On va rajouter une marge autour de nos hexagones/carrés
    float rho = sigma * (percentageRho);
    ;
    // ---------------------
    // Génération de la grille
    // ---------------------
    for (float y = 0.f; y < nH + sigma_y; y += sigma_y){
        for (float x = 0.f; x < nW; x += sigma){
            int x_center = static_cast<int>(x);
            int y_center = static_cast<int>(y);

            // On peut vérifier les bords
            if (x_center >= 0 && x_center < nW && y_center >= 0){
                centers.push_back({x_center, y_center});
            }
        }
    }


    std::vector<int> bestIndices;

    Image ImgInterGrilleColor = Image(nW, nH, true);
    // ---------------------
    // Affectation de chaque pixel au centre le plus proche
    // ---------------------
    for (int y = 0; y < nH; y++)
    {
        for (int x = 0; x < nW; x++)
        {
            int minDist = 1e9, bestIndex = 0;

            for (size_t i = 0; i < centers.size(); i++)
            {
                int dx = x - centers[i].x;
                int dy = y - centers[i].y;
                int dist = dx * dx + dy * dy;
                if (dist < minDist && (x < centers[i].x + rho && x > centers[i].x - rho) && (y < centers[i].y + rho && y
                    > centers[i].y - rho))
                {
                    minDist = dist;
                    bestIndex = i + 1;
                }
            }

            int index = (y * nW + x) * 3;
            int color = bestIndex % 256;
            ImgInterGrilleColor[index + 0] = (color * 3) % 256;
            ImgInterGrilleColor[index + 1] = (color * 7) % 256;
            ImgInterGrilleColor[index + 2] = (color * 5) % 256;
            bestIndices.push_back(bestIndex);
        }
    }
    if (debugMode){
        std::string InterGrilleColor = std::string( "WATERPIXEL-DEBUG_ImgInterGrilleColor.ppm");
        ImgInterGrilleColor.write(InterGrilleColor.c_str());
    }

    // ---------------------
    // Extraction des marqueurs
    // ---------------------
    std::vector<int> indicesCentresSuperPixels;
    for (int i = 0; i < nH; i++)
    {
        for (int j = 0; j < nW; j++)
        {
            int index = (i * nW + j);
            // Si je suis le pixel en haut à droite de la petite zone
            if ((i == 0 && j == 0 && bestIndices[index] != 0) ||
                (i == 0 && bestIndices[(i) * nW + j] != 0 && bestIndices[(i) * nW + j - 1] == 0) ||
                (j == 0 && bestIndices[(i) * nW + j] != 0 && bestIndices[(i - 1) * nW + j] == 0) ||
                (i > 0 && j > 0 && bestIndices[(i) * nW + j] != 0 && bestIndices[(i - 1) * nW + j] == 0 && bestIndices[i
                    * nW + j - 1] == 0)
            )
            {
                indicesCentresSuperPixels.push_back(getTheMarker(index, rho * 2, bestIndices, ImgNormeGradient));
            }
        }
    }


    // On a maintenant les indices des centres des superpixels dans indicesCentresSuperPixels et on va les utiliser
    // pour créer les ACM-waterpixels
    std::vector<SuperPixel> superPixels;
    for (int i = 0; i < indicesCentresSuperPixels.size(); i++)
    {
        Pixel *centroid = imagePixels[indicesCentresSuperPixels[i]];
        centroid->superpixel_id = superPixels.size()+1;// On commence à 1 pour éviter les confusions avec les pixels non assignés
        centroid->temp_id = superPixels.size()+1; // On commence à 1 pour éviter les confusions avec les pixels non assignés
        centroid->distance = 0;
        SuperPixel sp;
        sp.pixels.push_back(centroid);
        superPixels.push_back(sp);
    }
    Image ImgInterGrilleColorAfterChangeCenters = Image(nW, nH, true);


    // Boucle pour assigner chaque pixel à un superpixel en fonction de la distance en une passe
    // (formule fournis : dQ(p)=2/σ min(qi∈Q) d(p,qi) )

    // Si on prends les mots de l'article en disant que les m-waterpixels sont
    // quand on prends les bons marqueurs en prenant bien en compte les minima locaux
    // Et les c-waterpixels sont les pixels qui sont les centres de cellules
    // Alors les notres on peut les nommer les acm-waterpixels car on prends
    // pas exactement les minima locaux, ni les centres de cellules donc "average center marker"
    for (int i = 0; i < nH; i++)
    {
        for (int j = 0; j < nW; j++)
        {
            int index = i * nW + j;
            Pixel *pix = imagePixels[index];
            double minDist = std::numeric_limits<double>::max();
            int bestSuperPixel = -1;
            for (int spIdx = 0; spIdx < superPixels.size(); spIdx++)
            {
                SuperPixel &sp = superPixels[spIdx];
                Pixel &centroid = *sp.pixels[0];
                double distance = Pixel::computeSpaceDistance(*pix, centroid);
                distance = 2 * distance / sigma;
                if (distance < minDist){
                    minDist = distance;
                    bestSuperPixel = spIdx;
                }
            }
            pix->distance = minDist;
            pix->superpixel_id = bestSuperPixel+1; // On commence à 1 pour éviter les confusions avec les pixels non assignés

            // si le pixel est le centroid, on l'ajoute pas une seconde fois
            if (superPixels[bestSuperPixel].pixels[0]->x != pix->x || superPixels[bestSuperPixel].pixels[0]->y != pix->y)
                superPixels[bestSuperPixel].pixels.push_back(pix);

            ImgInterGrilleColorAfterChangeCenters[index*3 + 0] = labToRgb(superPixels[bestSuperPixel].pixels[0]->lab).r;
            ImgInterGrilleColorAfterChangeCenters[index*3 + 1] = labToRgb(superPixels[bestSuperPixel].pixels[0]->lab).g;
            ImgInterGrilleColorAfterChangeCenters[index*3 + 2] = labToRgb(superPixels[bestSuperPixel].pixels[0]->lab).b;
        }
    }

    // Ecriutre
    if (debugMode){
        std::string InterGrilleColorAfterChangeCenters = std::string("WATERPIXEL-DEBUG_ImgInterGrilleColorAfterChangeCenters.ppm");
        ImgInterGrilleColorAfterChangeCenters.write(InterGrilleColorAfterChangeCenters.c_str());
    }


    // Etape 5 : Régulariser spatialement le gradient
    // on va construire notre image g_reg qui est une image de gradient régularisé en utilisant les ACM-waterpixels
    // et en utilisant la formule fournis : g_reg(p) = g + reg*d(p, q) avec reg fournis par utilisateur
    // Plus reg est grand et plus il y a de la régularisation
    Image ImgNormeGradientReg = Image(nW, nH, false);
    std::vector<PixelNDG*> ImgNormeGradientRegVector = std::vector<PixelNDG*>(nW*nH); // On va stocker les pixels de l'image régularisée

    for (int spIdx = 0; spIdx < superPixels.size(); spIdx++){
        SuperPixel &sp = superPixels[spIdx];
        for (int pIdx = 0; pIdx < sp.pixels.size(); pIdx++){
            Pixel &pix = *sp.pixels[pIdx];
            double dQ = pix.distance;
            ImgNormeGradientReg[pix.y*nW+pix.x] = ImgNormeGradientVector[pix.y*nW+pix.x]->value + reg * dQ;
            if (ImgNormeGradientReg[pix.y*nW+pix.x]>255){
                ImgNormeGradientReg[pix.y*nW+pix.x] = 255;
            }
            else if (ImgNormeGradientReg[pix.y*nW+pix.x]<0)
            {
                ImgNormeGradientReg[pix.y*nW+pix.x] = 0;
            }

            PixelNDG *pixNDG = new PixelNDG();
            pixNDG->y = pix.y;
            pixNDG->x = pix.x;
            pixNDG->value = ImgNormeGradientReg[pix.y*nW+pix.x];
            ImgNormeGradientRegVector[pix.y*nW+pix.x] = pixNDG;
        }
    }

    if (debugMode){
        std::string NormeGradientReg = std::string("WATERPIXEL-DEBUG_ImgNormeGradientReg.pgm");
        ImgNormeGradientReg.write(NormeGradientReg.c_str());
    }


    // Création d'un vecteur labels qui contient les labels de chaque superpixel
    // Pour l'instant on va initialiser les labels à 0 pour les pixels qui ne sont pas les centroids
    // et à l'index du superpixel pour les pixels qui sont les centroids
    Image imgMarqueurCell2 = Image(nW, nH, false);
    std::vector<int> labels(nH*nW, 0);
    for (int spIdx = 0; spIdx < superPixels.size(); spIdx++){
        SuperPixel &sp = superPixels[spIdx];
        // On met le label du superpixel pour le pixel qui est le centroid
        labels[sp.pixels[0]->y*nW+sp.pixels[0]->x] = sp.pixels[0]->superpixel_id;
        imgMarqueurCell2[sp.pixels[0]->y*nW+sp.pixels[0]->x] = 255;
    }
    // Ecriture de marqueur cell 2
    if (debugMode){
        std::string MarqueurCell2 = std::string("WATERPIXEL-DEBUG_ImgMarqueurCell2.pgm");
        imgMarqueurCell2.write(MarqueurCell2.c_str());
    }



    // PARTIE WATERSHED : Pour cette partie, on va utiliser la norme du gradient régularisé calculé précédemment
    // Les marqueurs de notre watershed sont représentés par les centroide des ACM-waterpixels (superpixels temporaires)
    // Labels qui est un vecteur qui contient les labels de chaque pixel, qui a partir de la valeur 1 correspond à un superpixel
    // On va utiliser une file de priorité pour parcourir les pixels en fonction de leur valeur de gradient régularisé
    std::priority_queue<PixelNDG*, std::vector<PixelNDG*>, ComparePixelNDG> pq;

    // On ajoute les marqueurs à la file de priorité
    for (int spIdx = 0; spIdx < superPixels.size(); spIdx++){
        SuperPixel &sp = superPixels[spIdx];
        PixelNDG *pix = ImgNormeGradientRegVector[sp.pixels[0]->y*nW+sp.pixels[0]->x];
        ImgNormeGradientReg[sp.pixels[0]->y*nW+sp.pixels[0]->x] = 255;
        pq.push(pix);
    }

    // Processus d'inondation
    while (!pq.empty()){
        PixelNDG* p = pq.top();
        pq.pop();


        // Parcourir les voisins de p
        for (int di = -1; di <= 1; di++){
            for (int dj = -1; dj <= 1; dj++){
                if (di == 0 && dj == 0) continue; // On ignore le pixel courant
                int ni = p->y + di, nj = p->x + dj;
                // Si le voisin est dans l'image
                if (ni >= 0 && ni < nH && nj >= 0 && nj < nW){
                    PixelNDG *neighbor = ImgNormeGradientRegVector[ni*nW+nj];
                    if (labels[ni*nW+nj] == 0){
                        labels[ni*nW+nj] = labels[p->y*nW+p->x];
                        pq.push(neighbor);
                    }
                    // Si le voisin est déjà marqué par un autre superpixel et que je ne suis pas un watershed
                    else if (labels[ni*nW+nj] != labels[p->y*nW+p->x] && labels[p->y*nW+p->x] != -1){
                        int lab = labels[ni*nW+nj];
                        int labCompare = labels[p->y*nW+p->x];
                        //labels[ni*nW+nj] = -1; // Marquer comme watershed
                    }
                }
            }
        }
    }

    // On va réintialiser la liste des superpixels pour avoir une couleur adaptée à l'image de sortie
    for (int i = 0; i < superPixels.size(); i++){
        superPixels[i].pixels.clear();
    }
    // Grâce à labels, on va pouvoir reconstruire les superpixels
    for (int i = 0; i < labels.size(); i++){
        superPixels[labels[i]-1].pixels.push_back(imagePixels[i]);
    }

    // colors sert à stocker les couleurs moyennes des superpixels pour l'écriture de l'image de sortie
    std::vector<Rgb> colors;
    for(int i=0; i<superPixels.size(); i++){
        Lab lab = superPixels[i].getAveragePixel().lab;
        Rgb color = labToRgb(lab);
        colors.push_back(color);
    }
    // Ecriture de l'image de sortie
    // Ici on a les labels de chaque pixel, qui a partir de la valeur 1 correspond à un superpixel, on a plus qu'à reconstruire l'image
    // en mettant les couleurs moyennes des superpixels à la place des pixels
    imageOut = Image(imageIn.width, imageIn.height, true);
    for (int index = 0; index < labels.size(); index++){
        if (labels[index] == -1){
            imageOut[index*3 + 0] = 255;
            imageOut[index*3 + 1] = 255;
            imageOut[index*3 + 2] = 255;
        }
        else if (labels[index] == 0){
            imageOut[index*3 + 0] = 0;
            imageOut[index*3 + 1] = 0;
            imageOut[index*3 + 2] = 0;
        }
        else{
            int spIdx = labels[index] - 1;
            Rgb color = colors[spIdx];
            imageOut[index*3 + 0] = color.r;
            imageOut[index*3 + 1] = color.g;
            imageOut[index*3 + 2] = color.b;
        }
    }

    // Ecriture de l'image de sortie
    if (debugMode){
        std::string imgOutName = std::string("WATERPIXEL-DEBUG_ImgOut.ppm");
        imageOut.write(imgOutName.c_str());
    }

    // Libération de la mémoire
    for (Pixel * pixel: imagePixels){
        free(pixel);
    }
    for (PixelNDG * pixel: ImgNormeGradientRegVector){
        free(pixel);
    }

}