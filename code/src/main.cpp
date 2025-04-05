#include <image.h>
#include <super_pixel.h>
#include <util.h>
#include <huffman.h>
#include "imgui/imgui.h"
#include "imgui/backend/imgui_impl_glfw.h"
#include "imgui/backend/imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include <imgui/imfilebrowser.h>
#include <interface/helper.hpp>


InterfaceState currentState;

int main() {
    glfwInit();
    GLFWwindow* window = glfwCreateWindow(800, 600, "Mon App", NULL, NULL);
    glfwMakeContextCurrent(window);
    
    ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init();

    ImGui::FileBrowser fileDialog;
    
    // fileDialog.SetTitle("title");
    // fileDialog.SetTypeFilters({ ".h", ".cpp" });

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
    
        if(ImGui::Begin("Menu")) {
            if (ImGui::Button("File explo")) fileDialog.Open();
        }
        ImGui::End();
    
        fileDialog.Display();
        
        if(fileDialog.HasSelected()) {
            currentState.setInputImage(fileDialog.GetSelected().string().data());
            fileDialog.ClearSelected();
        }

        currentState.interfaceUpdate();
    
        // Clear le buffer avant de dessiner
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f); // Couleur de fond (gris foncé)
        glClear(GL_COLOR_BUFFER_BIT);
    
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();


    // if (argc == 5 || argc == 6) {
    //     Image imageIn,imageOutSNIC, imageOutSLIC,imageOutWaterpixel;

    //     std::cout << "Lecture de l'image" << std::endl;
    //     char filename[256];
    //     sscanf(argv[1], "%s", filename);
    //     imageIn.read(filename);

    //     int k = atoi(argv[3]);
    //     double m = atof(argv[4]);

    //     std::string SNICPath = std::string(argv[2]);
    //     SNICPath = SNICPath.substr(0, SNICPath.size() - 4) + "SNIC.ppm";
    //     std::cout << "Image superpixelisée \n SNIC écriture ici :"<< SNICPath << std::endl;
    //     SNIC(imageIn, imageOutSNIC, k, m);
    //     imageOutSNIC.write(SNICPath.c_str());


    //     std::string SLICPath = std::string(argv[2]);
    //     SLICPath = SLICPath.substr(0, SLICPath.size() - 4) + "SLIC.ppm";
    //     std::cout << "SLIC ecriture ici :"<< SLICPath << std::endl;

    //     SLIC(imageIn, imageOutSLIC, k, m);
    //     imageOutSLIC.write(SLICPath.c_str());


    //     std::string WaterpixelPath = std::string(argv[2]);
    //     WaterpixelPath = WaterpixelPath.substr(0, WaterpixelPath.size() - 4) + "Waterpixel.ppm";
    //     std::cout << "Waterpixel : ecriture ici : "<< WaterpixelPath << std::endl;

    //     Waterpixel(imageIn,imageOutWaterpixel,k,0.2f,0.);
    //     imageOutWaterpixel.write(WaterpixelPath.c_str());

    //     // Ecrire le PSNR dans un fichier argv[2] privé de .ppm remplacé par .psnr
    //     double psnrSNIC = PSNR(imageIn, imageOutSNIC);
    //     double psnrSLIC = PSNR(imageIn, imageOutSLIC);
    //     double psnrWaterpixel = PSNR(imageIn, imageOutWaterpixel);
    //     std::string psnrFile = std::string(argv[2]);
    //     psnrFile = psnrFile.substr(0, psnrFile.size() - 4) + ".psnr";
    //     std::ofstream file(psnrFile);
    //     file << "PSNR SNIC :"<< psnrSNIC<< "\nPSNR SLIC :"<< psnrSLIC << "\nPSNR WATERPIXEL : "<< psnrWaterpixel << std::endl;
    //     file.close();
    //     std::cout << "PSNR SNIC :"<< psnrSNIC<< "\nPSNR SLIC :"<< psnrSLIC << "\nPSNR WATERPIXEL : "<< psnrWaterpixel << std::endl;

    //     if (argc == 6){
    //         Image imageBoundaryRecall;
    //         imageBoundaryRecall.read(argv[5]);
    //         // Ecrire le calcul du boundary recall dans un fichier argv[2] privé de .ppm remplacé par .recall
    //         double recallSNIC = BoundaryRecall(imageOutSNIC, imageBoundaryRecall);
    //         double recallSLIC = BoundaryRecall(imageOutSLIC, imageBoundaryRecall);
    //         double recallWaterpixel = BoundaryRecall(imageOutWaterpixel, imageBoundaryRecall);
    //         std::string recallFile = std::string(argv[2]);
    //         recallFile = recallFile.substr(0, recallFile.size() - 4) + ".recall";
    //         std::ofstream file2(recallFile);
    //         file2 << "Recall SNIC :"<< recallSNIC << "\nRecall SLIC :"<< recallSLIC << "\nRecall Waterpixel:"<< recallWaterpixel << std::endl;
    //         file2.close();
    //         std::cout << "Boundary Recall SNIC :"<< recallSNIC << "\nRecall SLIC :"<< recallSLIC << "\nRecall Waterpixel:"<< recallWaterpixel << std::endl;
    //     }

    //     // Ecrire les histogrammes de couleurs dans des fichiers argv[2] privé de .ppm remplacé par .histo
    //     std::cout << "Histogramme de couleurs : " << std::endl;
    //     std::string histoFileSNIC = std::string(SNICPath);
    //     histoFileSNIC = histoFileSNIC.substr(0, histoFileSNIC.size() - 4) + ".histo";
    //     histoColor(imageOutSNIC, histoFileSNIC.c_str());

    //     std::string histoFileSLIC = std::string(SLICPath);
    //     histoFileSLIC = histoFileSLIC.substr(0, histoFileSLIC.size() - 4) + ".histo";
    //     histoColor(imageOutSLIC, histoFileSLIC.c_str());

    //     std::string histoFileWaterpixel = std::string(WaterpixelPath);
    //     histoFileWaterpixel = histoFileWaterpixel.substr(0, histoFileWaterpixel.size() - 4) + ".histo";
    //     histoColor(imageOutWaterpixel, histoFileWaterpixel.c_str());
    //     return 0;
    // } else if(argc>1){
    //     std::cerr << "Usage: " << argv[0] << " |OU| " << argv[0] << " imageIn.ppm imageOut.ppm k m" << " [imageBoundaryRecall.pgm (test de boundary recall)]"  << std::endl;
    //     return 1;
    // } else {
    //     Image imageIn, imageOut, imageOutSLIC, imageOutWaterpixel;
    //     imageIn.read("images/taupe.ppm");
        
    //     std::ofstream file("output/taupeComparison.dat");
    //     long double inSize = getFileSize("images/taupe.ppm");
    //     /*for(int k=1; k <= 150001; k += 10000){
    //         SNIC(imageIn, imageOut, k);
    //         SLIC(imageIn, imageOutSLIC, k);
    //         float psnrSNIC = PSNR(imageIn, imageOut);
    //         float psnrSLIC = PSNR(imageIn, imageOutSLIC);
            
            
    //         long double compressedSize = getFileSize("output/taupe.snic");
    //         auto compressionRate = inSize / compressedSize;
            
    //         file << k << " " << compressionRate << " " << psnrSNIC << " " << psnrSLIC << std::endl;
    //     }
    //     file.close();*/
        
    //     SNIC(imageIn, imageOut, 1000);
    //     SLIC(imageIn, imageOutSLIC, 1000);
    //     Waterpixel(imageIn, imageOutWaterpixel, 1000,0.2f, 5.,  true);
    //     imageOut.write("output/taupe_SNIC.ppm");
    //     imageOutSLIC.write("output/taupe_SLIC.ppm");
    //     imageOutWaterpixel.write("output/taupe_waterpixel.ppm");
    // }

    
    //////////////////////////////////////////////

    // double recall = BoundaryRecall(imageOut, imageOut);
    // std::cout << "Boundary Recall : " << recall << std::endl;
    
    return 0;
}