#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_IMPLEMENTATION
#include "image.h"
#include "imgui/imgui.h"
#include "super_pixel.h"

// Simple helper function to load an image into a OpenGL texture with common settings
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

int compressionMode = 0;
void displayImage(Image &image, GLuint textureId) {
    ImGui::Begin("OpenGL Texture Text");
    ImGui::Text("pointer = %x", textureId);
    ImGui::Text("size = %d x %d", image.width, image.height);
    
    ImGui::RadioButton("SLIC", &compressionMode, 0); ImGui::SameLine();
    ImGui::RadioButton("SNIC", &compressionMode, 1); ImGui::SameLine();
    ImGui::RadioButton("WATERPIXEL", &compressionMode, 2);

    if(ImGui::Button("Compresser")){
        Image imageOut;
        switch(compressionMode){
            case 0:
                SLIC(image, imageOut);
                break;
            case 1:
                SNIC(image, imageOut);
                break;
            case 2:
                Waterpixel(image, imageOut, 5000);
                break;
        }
    }
    

    ImGui::Image((ImTextureID)(intptr_t)textureId, ImVec2(image.width, image.height));
    ImGui::End();
}