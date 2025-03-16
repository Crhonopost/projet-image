#include <image.h>
#include <super_pixel.h>
#include <util.h>
#include <huffman.h>

int main(int argc, char* argv[])
{
    if (argc == 5 || argc == 6) {
        Image imageIn, imageOut;

        std::cout << "Lecture de l'image" << std::endl;
        char filename[256];
        sscanf(argv[1], "%s", filename);
        imageIn.read(filename);

        int k = atoi(argv[3]);
        double m = atof(argv[4]);

        std::cout << "Image superpixelisée" << std::endl;
        SNIC(imageIn, imageOut, k, m);
        imageOut.write(argv[2]);

        // Ecrire le PSNR dans un fichier argv[2] privé de .ppm remplacé par .psnr
        double psnr = PSNR(imageIn, imageOut);
        std::string psnrFile = std::string(argv[2]);
        psnrFile = psnrFile.substr(0, psnrFile.size() - 4) + ".psnr";
        std::ofstream file(psnrFile);
        file << psnr << std::endl;
        file.close();
        std::cout << "PSNR : " << psnr << std::endl;

        if (argc == 6){
            Image imageBoundaryRecall;
            imageBoundaryRecall.read(argv[5]);
            // Ecrire le calcul du boundary recall dans un fichier argv[2] privé de .ppm remplacé par .recall
            double recall = BoundaryRecall(imageOut, imageBoundaryRecall);
            std::string recallFile = std::string(argv[2]);
            recallFile = recallFile.substr(0, recallFile.size() - 4) + ".recall";
            std::ofstream file2(recallFile);
            file2 << recall << std::endl;
            file2.close();
            std::cout << "Boundary Recall : " << recall << std::endl;
        }



        // Ecrire les histogrammes de couleurs dans des fichiers argv[2] privé de .ppm remplacé par .histo
        std::cout << "Histogramme de couleurs : " << std::endl;
        std::string histoFile = std::string(argv[2]);
        histoFile = histoFile.substr(0, histoFile.size() - 4) + ".histo";
        histoColor(imageOut, histoFile.c_str());
        return 0;
    } else if(argc>1){
        std::cerr << "Usage: " << argv[0] << " |OU| " << argv[0] << " imageIn.ppm imageOut.ppm k m" << " [imageBoundaryRecall.pgm (test de boundary recall)]"  << std::endl;
        return 1;
    }

    Image imageIn, imageOut;
    imageIn.read("images/taupe.ppm");
    SNIC(imageIn, imageOut, 100000);

    /////// read back from compressed file ///////
    std::vector<unsigned char> data;
    decompressFile("output/taupe.snic", data);
    readSNIC(data, imageOut);
    imageOut.write("output/res.ppm");

    long double inSize = getFileSize("images/taupe.ppm");
    long double compressedSize = getFileSize("output/taupe.snic");

    auto compressionRate = inSize / compressedSize;

    std::cout << "taille d'origine: " << inSize / 1000. << "kb  /  taille compressée: " << compressedSize / 1000. << "kb\n";
    std::cout << "taux de compression: " << compressionRate << "\n";
    //////////////////////////////////////////////

    // double recall = BoundaryRecall(imageOut, imageOut);
    // std::cout << "Boundary Recall : " << recall << std::endl;
    
    return 1;
}