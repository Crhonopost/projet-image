#include <image.h>
#include <super_pixel.h>
#include <util.h>
int main(int argc, char* argv[])
{
    if (argc == 5) {
        Image imageIn, imageOut;
        char filename[256];
        sscanf(argv[1], "%s", filename);
        imageIn.read(filename);
        int k = atoi(argv[3]);
        double m = atof(argv[4]);
        std::cout << "Lecture de l'image" << std::endl;
        SNIC(imageIn, imageOut, k, m);
        std::cout << "Image superpixelisée" << std::endl;
        imageOut.write(argv[2]);
        // Ecrire le PSNR dans un fichier argv[2] privé de .ppm remplacé par .psnr
        double psnr = PSNR(imageIn, imageOut);
        std::string psnrFile = std::string(argv[2]);
        psnrFile = psnrFile.substr(0, psnrFile.size() - 4) + ".psnr";
        std::ofstream file(psnrFile);
        file << psnr << std::endl;
        file.close();
        std::cout << "PSNR : " << psnr << std::endl;

        // Ecrire les histogrammes de couleurs dans des fichiers argv[2] privé de .ppm remplacé par .histo
        std::cout << "Histogramme de couleurs : " << std::endl;
        std::string histoFile = std::string(argv[2]);
        histoFile = histoFile.substr(0, histoFile.size() - 4) + ".histo";
        histoColor(imageOut, histoFile.c_str());
        return 0;
    }

    Image imageIn, imageOut;
    imageIn.read("images/taupe.ppm");
    
    SNIC(imageIn, imageOut);
    
    return 1;
}