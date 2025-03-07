#include <Image.h>
#include <super_pixel.h>

int main(int argc, char* argv[])
{
    Image imageIn, imageOut;
    imageIn.read("images/taupe.ppm");
    
    SNIC(imageIn, imageOut);
    
    return 1;
}