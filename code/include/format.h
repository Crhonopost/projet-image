#pragma once

#include <cmath>

#include <image.h>


struct Rgb {
    OCTET r,g,b;
    Rgb(): r(0), g(0), b(0){}
    Rgb(OCTET r, OCTET g, OCTET b){
        this->r = r;
        this->g = g;
        this->b = b;
    }
};

struct Xyz {
    double x,y,z;
};

struct Lab {
    double l, a, b;
    Lab(){
        l = a = b = 0;
    }

    void add(Lab &other){
        l += other.l;
        a += other.a;
        b += other.b;
    }

    void div(double value){
        l /= value;
        a /= value;
        b /= value;
    }
};

Xyz rgbToXyz(Rgb rgb){
    double R = rgb.r / 255.0;
    double G = rgb.g / 255.0;
    double B = rgb.b / 255.0;

    // Application de la transformation gamma inverse (standard de l'espace sRGB)
    R = (R > 0.04045) ? std::pow((R + 0.055) / 1.055, 2.4) : R / 12.92;
    G = (G > 0.04045) ? std::pow((G + 0.055) / 1.055, 2.4) : G / 12.92;
    B = (B > 0.04045) ? std::pow((B + 0.055) / 1.055, 2.4) : B / 12.92;


    // Transformation RGB → XYZ
    double X = R * 0.4124564 + G * 0.3575761 + B * 0.1804375;
    double Y = R * 0.2126729 + G * 0.7151522 + B * 0.0721750;
    double Z = R * 0.0193339 + G * 0.1191920 + B * 0.9503041;

    Xyz res;
    res.x = X;
    res.y = Y;
    res.z = Z;

    return res;
}

Lab rgbToLab(Rgb rgb) {
    Xyz xyz = rgbToXyz(rgb);

    double epsilon = 0.008856; //actual CIE standard
    double kappa = 903.3; //actual CIE standard


    // Xr = 0.950456 
    // Yr = 1.0 
    // Zr = 1.088754 

    xyz.x /= 0.95047;
    xyz.y /= 1.00000;
    xyz.z /= 1.08883;

    // Fonction de transformation f(t)
    auto f = [](double t) -> double {
        return (t > 0.008856) ? std::pow(t, 1.0 / 3.0) : (t * 903.3 + 16.0) / 116.0;
    };

    double X = f(xyz.x);
    double Y = f(xyz.y);
    double Z = f(xyz.z);

    // Calcul du CIELAB
    Lab lab;
    lab.l = 116.0 * Y - 16.0; // (Y > epsilon) ? (116.0 * std::pow(Y, 1.0 / 3.0) - 16.0) : (kappa * Y);
    lab.a = 500.0 * (X - Y);
    lab.b = 200.0 * (Y - Z);

    return lab;
}

Xyz labToXyz(Lab lab) {
    auto f_inv = [](double t) -> double {
        double t3 = t * t * t;
        return (t3 > 0.008856) ? t3 : (116.0 * t - 16.0) / 903.3;
    };
    
    double Y = (lab.l + 16.0) / 116.0;
    double X = lab.a / 500.0 + Y;
    double Z = Y - lab.b / 200.0;
    
    Xyz xyz;
    xyz.x = f_inv(X) * 0.95047;
    xyz.y = f_inv(Y) * 1.00000;
    xyz.z = f_inv(Z) * 1.08883;
    
    return xyz;
}

Rgb labToRgb(Lab lab) {
    Xyz xyz = labToXyz(lab);
    
    double R = xyz.x *  3.2404542 + xyz.y * -1.5371385 + xyz.z * -0.4985314;
    double G = xyz.x * -0.9692660 + xyz.y *  1.8760108 + xyz.z *  0.0415560;
    double B = xyz.x *  0.0556434 + xyz.y * -0.2040259 + xyz.z *  1.0572252;
    
    auto gamma_correction = [](double c) -> double {
        return (c > 0.0031308) ? (1.055 * std::pow(c, 1.0 / 2.4) - 0.055) : 12.92 * c;
    };
    
    R = gamma_correction(R);
    G = gamma_correction(G);
    B = gamma_correction(B);
    
    auto clamp = [](double v) -> int {
        return static_cast<int>(std::round(std::max(0.0, std::min(1.0, v)) * 255.0));
    };
    
    Rgb rgb;
    rgb.r = clamp(R);
    rgb.g = clamp(G);
    rgb.b = clamp(B);
    
    return rgb;
}



