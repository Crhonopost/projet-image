#pragma once

#include <image.h>
#include <huffman.h>
#include <fstream>

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

int compressFile(bool compress, char *srcPath, char *dstPath)
{
   FILE *src, *dst, *frq;

   /* Ouverture du fichier source en lecture */
   if ((src=fopen(srcPath, "rb"))==NULL)
   {
      perror("fopen");
      return -1;
   }

   /* Ouverture du fichier cible en écriture */
   if ((dst=fopen(dstPath, "wb"))==NULL)
   {
      perror("fopen");
      return -1;
   }

   frq=NULL;


   /* Allocation mémoire pour les données diverses necessaires à la construction de l'arbre */
   if ((arbre_d=(struct arbre_data *)malloc(512*sizeof(struct arbre_data)))==NULL)
   {
      perror("malloc");
      return -1;
   }

   /* Allocation d'une zone mémoire pour l'arbre */
   if ((arbre=(struct arbre *)malloc(512*sizeof(struct arbre)))==NULL)
   {
      free(arbre_d);
      perror("malloc");
      return -1;
   }

   /* Compression ou décompression ? */
   if (compress)
      huffman_compacter(src, dst, frq);
   else
      huffman_decompacter(src, dst, frq);
//    else if (*argv[1]=='f')
//       huffman_creer_fichier_frequences(src, dst);

   /* Libération de la mémoire */
   free(arbre_d);
   free(arbre);

   /* Fermeture des fichiers */
   fclose(src);
   fclose(dst);
   if (frq!=NULL)
      fclose(frq);
#ifdef LINUX_COMPIL
   printf("Linux rules!\n");
#endif
   return 0;
}



