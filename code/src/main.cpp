#include <image.h>
#include <super_pixel.h>
#include <util.h>
#include <huffman.h>

int main(int argc, char* argv[])
{
    if (argc == 5 || argc == 6) {
        Image imageIn,imageOutSNIC, imageOutSLIC,imageOutWaterpixel;

        std::cout << "Lecture de l'image" << std::endl;
        char filename[256];
        sscanf(argv[1], "%s", filename);
        imageIn.read(filename);

        int k = atoi(argv[3]);
        double m = atof(argv[4]);

        // -------------------PARTIE SNIC--------------------
        // Chemin de l'image de compression
        std::string pathCompressedSNIC = std::string(argv[2]);
        pathCompressedSNIC = pathCompressedSNIC.substr(0, pathCompressedSNIC.size() - 4) + "Compressed_SNIC";

        char *pathCharSNIC = new char[pathCompressedSNIC.length() + 1];
        std::strcpy(pathCharSNIC, pathCompressedSNIC.c_str());

        // Executer la fonction SNIC
        SNIC(imageIn, imageOutSNIC, k, m,pathCharSNIC);
        // Chemin de l'image de sortie
        std::string SNICPath = std::string(argv[2]);
        SNICPath = SNICPath.substr(0, SNICPath.size() - 4) + "SNIC.ppm";
        std::cout << "SNIC écriture ici :"<< SNICPath << std::endl;
        // Ecrire l'image de sortie
        imageOutSNIC.write(SNICPath.c_str());
        // Recuperer la taille de l'image compressée
        long double compressedSizeSNIC = getFileSize(pathCharSNIC);

        // -------------------PARTIE SLIC--------------------
        // Chemin de l'image de compression
        std::string pathCompressedSLIC = std::string(argv[2]);
        pathCompressedSLIC = pathCompressedSLIC.substr(0, pathCompressedSLIC.size() - 4) + "Compressed_SLIC";
        char *pathCharSLIC = new char[pathCompressedSLIC.length() + 1];
        std::strcpy(pathCharSLIC, pathCompressedSLIC.c_str());

        // Executer la fonction SLIC
        SLIC(imageIn, imageOutSLIC, k, m, pathCharSLIC);
        // Chemin de l'image de sortie
        std::string SLICPath = std::string(argv[2]);
        SLICPath = SLICPath.substr(0, SLICPath.size() - 4) + "SLIC.ppm";
        std::cout << "SLIC ecriture ici :"<< SLICPath << std::endl;
        // Ecrire l'image de sortie
        imageOutSLIC.write(SLICPath.c_str());
        // Recuperer la taille de l'image compressée
        long double compressedSizeSLIC = getFileSize(pathCharSLIC);

        // -------------------PARTIE WATERPIXEL--------------------
        // Chemin de l'image de compression
        std::string pathCompressedWaterpixel = std::string(argv[2]);
        pathCompressedWaterpixel = pathCompressedWaterpixel.substr(0, pathCompressedWaterpixel.size() - 4) + "Compressed_Waterpixel";
        char *pathCharWaterpixel = new char[pathCompressedWaterpixel.length() + 1];
        std::strcpy(pathCharWaterpixel, pathCompressedWaterpixel.c_str());
        // Executer la fonction Waterpixel
        Waterpixel(imageIn,imageOutWaterpixel,k,0.2f,m,false, pathCharWaterpixel);
        // Chemin de l'image de sortie
        std::string WaterpixelPath = std::string(argv[2]);
        WaterpixelPath = WaterpixelPath.substr(0, WaterpixelPath.size() - 4) + "Waterpixel.ppm";
        std::cout << "Waterpixel : ecriture ici : "<< WaterpixelPath << std::endl;
        imageOutWaterpixel.write(WaterpixelPath.c_str());
        // Recuperer la taille de l'image compressée
        long double compressedSizeWaterpixel = getFileSize(pathCharWaterpixel);

        // Ecrire le PSNR dans un fichier argv[2] privé de .ppm remplacé par .data
        double psnrSNIC = PSNR(imageIn, imageOutSNIC);
        double psnrSLIC = PSNR(imageIn, imageOutSLIC);
        double psnrWaterpixel = PSNR(imageIn, imageOutWaterpixel);

        // Ouvrir le fichier de sortie
        std::string dataFile = std::string(argv[2]);
        dataFile = dataFile.substr(0, dataFile.size() - 4) + ".data";
        std::ofstream file(dataFile);
        if (argc == 6){
            Image imageBoundaryRecallSNIC, imageBoundaryRecallSLIC, imageBoundaryRecallWaterpixel;
            imageBoundaryRecallSNIC.read(argv[5]);
            imageBoundaryRecallSLIC.read(argv[5]);
            imageBoundaryRecallWaterpixel.read(argv[5]);

            // Ecrire le calcul du boundary recall dans un fichier argv[2] privé de .ppm remplacé par .data
            double recallSNIC = BoundaryRecall(imageOutSNIC, imageBoundaryRecallSNIC,
                    SNICPath);
            double recallSLIC = BoundaryRecall(imageOutSLIC, imageBoundaryRecallSLIC,SLICPath);
            double recallWaterpixel = BoundaryRecall(imageOutWaterpixel, imageBoundaryRecallWaterpixel,WaterpixelPath);
            file <<"\nSNIC : k :" << k << " m :"<< m << " psnr :"<< psnrSNIC << " SizeCompressed : "<<compressedSizeSNIC <<" Recall :"<< recallSNIC<<
                    "\n" << "SLIC : k :" << k << " m :"<< m <<" psnr :"<< psnrSLIC<< " SizeCompressed : "<<compressedSizeSLIC << " Recall :"<< recallSLIC<<
                    "\n" << "WATERPIXEL : k :" << k << " m :"<< m <<" psnr :"<< psnrWaterpixel<< " SizeCompressed : "<<compressedSizeWaterpixel << " Recall :"<< recallWaterpixel
            << std::endl;
            file.close();
        }
        else{
            file <<"\nSNIC : k :" << k << " m :"<< m << " psnr :"<< psnrSNIC << "\n" << "SLIC : k :" << k << " m :"<< m <<" psnr :"<< psnrSLIC << "\n" << "WATERPIXEL : k :" << k << " m :"<< m <<"psnr :"<< psnrWaterpixel << std::endl;
            std::cout << "PSNR SNIC :"<< psnrSNIC<< "\nPSNR SLIC :"<< psnrSLIC << "\nPSNR WATERPIXEL : "<< psnrWaterpixel << std::endl;
        }
        file.close();

        // Ecrire les histogrammes de couleurs dans des fichiers argv[2] privé de .ppm remplacé par .histo
        /*    std::cout << "Histogramme de couleurs : " << std::endl;
            std::string histoFileSNIC = std::string(SNICPath);
            //histoFileSNIC = histoFileSNIC.substr(0, histoFileSNIC.size() - 4) + ".histo";
            //histoColor(imageOutSNIC, histoFileSNIC.c_str());

            std::string histoFileSLIC = std::string(SLICPath);
            histoFileSLIC = histoFileSLIC.substr(0, histoFileSLIC.size() - 4) + ".histo";
            //histoColor(imageOutSLIC, histoFileSLIC.c_str());

            std::string histoFileWaterpixel = std::string(WaterpixelPath);
            //histoFileWaterpixel = histoFileWaterpixel.substr(0, histoFileWaterpixel.size() - 4) + ".histo";
            //histoColor(imageOutWaterpixel, histoFileWaterpixel.c_str());*/
        return 0;
    } else if(argc>1){
        std::cerr << "Usage: " << argv[0] << " |OU| " << argv[0] << " imageIn.ppm imageOut.ppm k m" << " [imageBoundaryRecall.pgm (test de boundary recall)]"  << std::endl;
        return 1;
    } else {
        Image imageIn, imageOut, imageOutSLIC, imageOutWaterpixel;
        imageIn.read("images/taupe.ppm");
        
        std::ofstream file("output/taupeComparison.dat");
        long double inSize = getFileSize("images/taupe.ppm");
        /*for(int k=1; k <= 150001; k += 10000){
            SNIC(imageIn, imageOut, k);
            SLIC(imageIn, imageOutSLIC, k);
            float psnrSNIC = PSNR(imageIn, imageOut);
            float psnrSLIC = PSNR(imageIn, imageOutSLIC);
            
            
            long double compressedSize = getFileSize("output/taupe.snic");
            auto compressionRate = inSize / compressedSize;
            
            file << k << " " << compressionRate << " " << psnrSNIC << " " << psnrSLIC << std::endl;
        }
        file.close();*/
        
        SNIC(imageIn, imageOut, 1000);
        SLIC(imageIn, imageOutSLIC, 1000);
        Waterpixel(imageIn, imageOutWaterpixel, 1000,0.2f, 5.,  true);
        imageOut.write("output/taupe_SNIC.ppm");
        imageOutSLIC.write("output/taupe_SLIC.ppm");
        imageOutWaterpixel.write("output/taupe_waterpixel.ppm");
    }

    
    //////////////////////////////////////////////

    // double recall = BoundaryRecall(imageOut, imageOut);
    // std::cout << "Boundary Recall : " << recall << std::endl;
    
    return 0;
}