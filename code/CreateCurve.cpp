#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cstdio>    // pour remove()
#include <cstring>

struct DataPoint {
    int k;
    double psnr;
    double sizeCompressed;
    double recall;
};

void plotCurvesByM_All(const char* inputFile,
                       const char* outputPlotPsnr,
                       const char* outputPlotSize,
                       const char* outputPlotRecall,
                       const char* outputPlotPsnrSize) {
    // Map : clé = valeur de m, valeur = vecteur de DataPoint
    std::map<double, std::vector<DataPoint>> groups;

    // Lecture des fichiers
    std::ifstream infile(inputFile);
    if (!infile.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << inputFile << std::endl;
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty())
            continue;
        std::istringstream iss(line);
        int k;
        double m, psnr;
        double sizeCompressed;
        double recall;
        if (!(iss >> k >> m >> psnr >> sizeCompressed >> recall)) {
            std::cerr << "Erreur de parsing sur la ligne : " << line << std::endl;
            continue;
        }
        groups[m].push_back({k, psnr, sizeCompressed, recall});
    }
    infile.close();

    // Pour chaque groupe (valeur de m), trier par k et sauvegarder dans des fichiers temporaires
    // Les noms des fichiers seront de la forme : plot_psnr_m_<m>.dat, plot_size_m_<m>.dat, plot_recall_m_<m>.dat
    for (auto &entry : groups) {
        double mVal = entry.first;
        auto &points = entry.second;
        std::sort(points.begin(), points.end(), [](const DataPoint &a, const DataPoint &b) {
            return a.k < b.k;
        });

        // Construire les noms de fichiers
        std::ostringstream ossPsnr, ossSize, ossRecall, ossPsnrSize;
        ossPsnr << "plot_psnr_m_" << mVal << ".dat";
        ossSize << "plot_size_m_" << mVal << ".dat";
        ossRecall << "plot_recall_m_" << mVal << ".dat";
        ossPsnrSize << "plot_psnr_size_m_" << mVal << ".dat";
        std::string filePsnr = ossPsnr.str();
        std::string fileSize = ossSize.str();
        std::string fileRecall = ossRecall.str();
        std::string filePsnrSize = ossPsnrSize.str();
        {
            std::ofstream outPsnr(filePsnr);
            if (!outPsnr.is_open()) {
                std::cerr << "Erreur lors de la création de " << filePsnr << std::endl;
                continue;
            }
            for (const auto &pt : points)
                outPsnr << pt.k << " " << pt.psnr << "\n";
        }
        {
            std::ofstream outSize(fileSize);
            if (!outSize.is_open()) {
                std::cerr << "Erreur lors de la création de " << fileSize << std::endl;
                continue;
            }
            for (const auto &pt : points)
                outSize << pt.k << " " << pt.sizeCompressed << "\n";
        }
        {
            std::ofstream outRecall(fileRecall);
            if (!outRecall.is_open()) {
                std::cerr << "Erreur lors de la création de " << fileRecall << std::endl;
                continue;
            }
            for (const auto &pt : points)
                outRecall << pt.k << " " << pt.recall << "\n";
        }
        {
            std::ofstream outPsnrSize(filePsnrSize);
            if (!outPsnrSize.is_open()) {
                std::cerr << "Erreur lors de la création de " << filePsnrSize << std::endl;
                continue;
            }
            for (const auto &pt : points)
                outPsnrSize << pt.psnr << " " << pt.sizeCompressed << "\n";
        }
    }

    // Construction des commandes gnuplot pour les quatres graphes
    std::string gnuplotCmdPsnr = "gnuplot -persist -e \"set term png; set output '";
    gnuplotCmdPsnr += outputPlotPsnr;
    gnuplotCmdPsnr += ".png'; plot ";

    std::string gnuplotCmdSize = "gnuplot -persist -e \"set term png; set output '";
    gnuplotCmdSize += outputPlotSize;
    gnuplotCmdSize += ".png'; plot ";

    std::string gnuplotCmdRecall = "gnuplot -persist -e \"set term png; set output '";
    gnuplotCmdRecall += outputPlotRecall;
    gnuplotCmdRecall += ".png'; plot ";

    std::string gnuplotCmdPsnrSize = "gnuplot -persist -e \"set term png; set output '";
    gnuplotCmdPsnrSize += outputPlotPsnrSize;
    gnuplotCmdPsnrSize += ".png'; plot";

    bool first = true;
    for (auto &entry : groups) {
        double mVal = entry.first;
        std::ostringstream ossPsnr, ossSize, ossRecall,ossPsnrSize;
        ossPsnr << "plot_psnr_m_" << mVal << ".dat";
        ossSize << "plot_size_m_" << mVal << ".dat";
        ossRecall << "plot_recall_m_" << mVal << ".dat";
        ossPsnrSize << "plot_psnr_size_m_"<< mVal << ".dat";
        std::string filePsnr = ossPsnr.str();
        std::string fileSize = ossSize.str();
        std::string fileRecall = ossRecall.str();
        std::string filePsnrSize = ossPsnrSize.str();

        if (!first) {
            gnuplotCmdPsnr += ", ";
            gnuplotCmdSize += ", ";
            gnuplotCmdRecall += ", ";
            gnuplotCmdPsnrSize += ", ";
        } else {
            first = false;
        }
        // Pour psnr
        gnuplotCmdPsnr += "'" + filePsnr + "' using 1:2 with lines title 'm = " + std::to_string(mVal) + "'";
        // Pour taille
        gnuplotCmdSize += "'" + fileSize + "' using 1:2 with lines title 'm = " + std::to_string(mVal) + "'";
        // Pour recall
        gnuplotCmdRecall += "'" + fileRecall + "' using 1:2 with lines title 'm = " + std::to_string(mVal) + "'";
        // Pour psnr size
        gnuplotCmdPsnrSize += "'" + filePsnrSize + "' using 1:2 with lines title 'm = " + std::to_string(mVal) + "'";

    }
    gnuplotCmdPsnr += "\"";
    gnuplotCmdSize += "\"";
    gnuplotCmdRecall += "\"";
    gnuplotCmdPsnrSize += "\"";

    // Exécution des commandes gnuplot
    std::cout << "Commande gnuplot PSNR : " << gnuplotCmdPsnr << std::endl;
    system(gnuplotCmdPsnr.c_str());

    std::cout << "Commande gnuplot Taille : " << gnuplotCmdSize << std::endl;
    system(gnuplotCmdSize.c_str());

    std::cout << "Commande gnuplot Recall : " << gnuplotCmdRecall << std::endl;
    system(gnuplotCmdRecall.c_str());

    std::cout << "Commande gnuplot Recall : " << gnuplotCmdPsnrSize << std::endl;
    system(gnuplotCmdPsnrSize.c_str());

    // Supprime des fichiers qui ont permis de plot les graphes
    for (auto &entry : groups) {
        double mVal = entry.first;
        std::ostringstream ossPsnr, ossSize, ossRecall, ossPsnrSize;
        ossPsnr << "plot_psnr_m_" << mVal << ".dat";
        ossSize << "plot_size_m_" << mVal << ".dat";
        ossRecall << "plot_recall_m_" << mVal << ".dat";
        ossPsnrSize << "plot_psnr_size_m_"<< mVal << ".dat";
        std::remove(ossPsnr.str().c_str());
        std::remove(ossSize.str().c_str());
        std::remove(ossRecall.str().c_str());
        std::remove(ossPsnrSize.str().c_str());
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " input_file" << std::endl;
        return 1;
    }

    const char* inputFile = argv[1];
    std::string inputFileStr = argv[1];
    std::string outputPlotPsnr = std::string(inputFile) + "_psnr";
    std::string outputPlotSize = std::string(inputFile) + "_size";
    std::string outputPlotRecall = std::string(inputFile) + "_recall";
    std::string outputPlotPsnrSize = std::string(inputFile) + "_psnr_size";

    char *charPathPsnr = new char[outputPlotPsnr.size() + 1];
    std::strcpy(charPathPsnr, outputPlotPsnr.c_str());

    char *charPathSize = new char[outputPlotSize.size() + 1];
    std::strcpy(charPathSize, outputPlotSize.c_str());

    char *charPathRecall = new char[outputPlotRecall.size() + 1];
    std::strcpy(charPathRecall, outputPlotRecall.c_str());

    char *charPathPsnrSize = new char[outputPlotPsnrSize.size() +1];
    std::strcpy(charPathPsnrSize, outputPlotPsnrSize.c_str());

    plotCurvesByM_All(inputFile, charPathPsnr, charPathSize, charPathRecall,charPathPsnrSize);

    return 0;
}