#include <Image.h>

int main(int argc, char* argv[])
{
    Image test;
    test.read("images/taupe.ppm");
    
    OCTET premierPixel[3] = {test[0], test[1], test[2]};
    std::cout << "premier pixel: " << (int) premierPixel[0]<< " " << (int) premierPixel[1] << " " << (int) premierPixel[2] << "\n";
    
    test.write("output/test.ppm");
    return 1;
}