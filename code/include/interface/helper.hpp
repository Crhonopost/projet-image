#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_IMPLEMENTATION
#include "image.h"
#include "imgui/imgui.h"
#include "super_pixel.h"


struct InterfaceState {
    int compressionMode = 0;
    Image imgInput, imgOutput;
    GLuint inputId, outputId;
    bool inputDisplay, outputDisplay;

    int k=5000;
    float m=5;

    float percentageRho = 0.2f;
    double reg = 50.0;

    char *outputPath = "output/compressed.spr";
    float psnr;
    double inSize, compressionRate;

    bool LoadTextureFromMemory(Image &image, GLuint* out_texture)
    {    
        // Create a OpenGL texture identifier
        GLuint image_texture;
        glGenTextures(1, &image_texture);
        glBindTexture(GL_TEXTURE_2D, image_texture);

        // Setup filtering parameters for display
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        // Upload pixels into texture
        glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
        GLenum format = image.isColor ? GL_RGB : GL_RED;
        glTexImage2D(GL_TEXTURE_2D, 0, format, image.width, image.height, 0, format, GL_UNSIGNED_BYTE, image.data);

        *out_texture = image_texture;

        return true;
    }

    void displayImageInput(Image &image, GLuint textureId) {
        ImGui::Text("dimmensions = %d x %d", image.width, image.height);
        
        ImGui::RadioButton("SLIC", &compressionMode, 0); ImGui::SameLine();
        ImGui::RadioButton("SNIC", &compressionMode, 1); ImGui::SameLine();
        ImGui::RadioButton("WATERPIXEL", &compressionMode, 2);

        if(compressionMode == 0 || compressionMode == 1){
            ImGui::Text("Paramètres");
            ImGui::InputInt("valeur de k", &k);
            ImGui::InputFloat("valeur de m", &m);
        } else if(compressionMode == 2){
            ImGui::InputInt("valeur de k", &k);
            ImGui::SliderFloat("pourcentage rho", &percentageRho, 0.01f, 1.f);
            ImGui::InputDouble("reg", &reg);
        }
    
        if(ImGui::Button("Compresser")){
            Image imageOut;
            switch(compressionMode){
                case 0:
                    SLIC(image, imageOut, k, m, outputPath);
                    break;
                case 1:
                    SNIC(image, imageOut, k, m, outputPath);
                    break;
                case 2:
                    Waterpixel(image, imageOut, k, percentageRho, reg, false, outputPath);
                    break;
            }
            setOutputImage(imageOut);
        }
        
    
        ImGui::Image((ImTextureID)(intptr_t)textureId, ImVec2(image.width, image.height));
    }

    void interfaceUpdate(){
        if(inputDisplay){
            displayImageInput(imgInput, inputId);
        }

        if(outputDisplay){
            ImGui::SameLine();
            ImGui::Image((ImTextureID)(intptr_t)outputId, ImVec2(imgOutput.width, imgOutput.height));

            ImGui::Text("Image enregistrée sous '%s", outputPath);

            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << psnr;  // 2 décimales
            ImGui::Text("PSNR = %sdB", oss.str().c_str());

            std::ostringstream crs;
            crs << std::fixed << std::setprecision(2) << compressionRate;  // 2 décimales
            ImGui::Text("Taux de compression = %s", crs.str().c_str());


        }
    }

    void setInputImage(char *path){
        inSize = getFileSize(path);
        imgInput.read(path);
        inputDisplay = true;
        LoadTextureFromMemory(imgInput, &inputId);
    }

    void setOutputImage(Image image){
        imgOutput = image;
        outputDisplay = true;
        psnr = PSNR(imgInput, imgOutput);
        double outSize = getFileSize(outputPath);
        compressionRate = inSize / outSize;
        LoadTextureFromMemory(imgOutput, &outputId);
    }
};