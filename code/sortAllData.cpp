#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib> // pour system()

struct Record {
    int k;
    double m;
    double psnr;
    int sizeCompressed;
    double recall;
};

struct Accumulator {
    double sumPsnr = 0.0;
    double sumSize = 0.0;
    double sumRecall = 0.0;
    int count = 0;
};

int main(int argc, char* argv[]) {
    if(argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input_file.txt" << std::endl;
        return 1;
    }

    std::ifstream infile(argv[1]);
    if(!infile) {
        std::cerr << "Erreur d'ouverture du fichier " << argv[1] << std::endl;
        return 1;
    }

    // Map : clé = pair (algorithme, number), valeur = vecteur de Record
    std::map<std::pair<std::string, std::string>, std::vector<Record>> data;
    std::string line;
    std::string currentNumber; // extrait du chemin (ex: "42049")

    // Regex pour extraire le numéro dans le chemin (ex: "./output_stat/42049/42049_100_2.0.psnr :")
    std::regex fileRegex(R"(\./output_stat/(\d+)/)");

    // Regex pour extraire une ligne contenant toutes les informations
    // Exemple attendu :
    // SNIC : k :100 m :0.5 psnr :21.4715 SizeCompressed : 50927 Recall :0.719965
    std::regex algoRegex(R"(^\s*(\w+)\s*:\s*k\s*:\s*([\d\.]+)\s*m\s*:\s*([\d\.]+)\s*psnr\s*:\s*([\d\.]+)\s*SizeCompressed\s*:\s*([\d\.]+)\s*Recall\s*:\s*([\d\.]+))");

    while(std::getline(infile, line)) {
        std::smatch match;
        if(std::regex_search(line, match, fileRegex)) {
            currentNumber = match[1];
            continue;
        }
        if(std::regex_search(line, match, algoRegex)) {
            std::string algo = match[1];
            int k = std::stoi(match[2]);
            double m = std::stod(match[3]);
            double psnr = std::stod(match[4]);
            int sizeCompressed = std::stoi(match[5]);
            double recall = std::stod(match[6]);

            Record rec { k, m, psnr, sizeCompressed, recall };
            data[{algo, currentNumber}].push_back(rec);
        }
    }
    infile.close();

    // Pour calculer la moyenne par algorithme pour chaque paire (k, m)
    // Map : clé = algorithme, valeur = map avec clé (k, m) et valeur Accumulator
    std::map<std::string, std::map<std::pair<int, double>, Accumulator>> avgData;

    // Pour chaque couple (algorithme, number), trier par k puis par m et écrire dans un fichier.
    for(auto& kv : data) {
        std::string algo = kv.first.first;
        std::string number = kv.first.second;
        auto& records = kv.second;

        // Tri par k croissant, puis par m décroissant (selon votre version initiale)
        std::sort(records.begin(), records.end(), [](const Record& a, const Record& b) {
            return (a.k == b.k) ? (a.m > b.m) : (a.k < b.k);
        });

        // On crée le dossier de sortie pour cet algorithme
        std::string dirName = "all_data_values_sorted_" + algo;
        std::string mkdirCmd = "mkdir -p " + dirName;
        system(mkdirCmd.c_str());

        std::string outFileName = dirName + "/all_data_sorted_" + number;
        std::ofstream outfile(outFileName);
        if(!outfile) {
            std::cerr << "Erreur lors de la création du fichier " << outFileName << std::endl;
            continue;
        }

        // Pour chaque enregistrement, on écrit dans le fichier individuel
        // et on accumule les valeurs pour le calcul de la moyenne.
        for(const auto& rec : records) {
            outfile << rec.k << " "
                    << rec.m << " "
                    << rec.psnr << " "
                    << rec.sizeCompressed << " "
                    << rec.recall << "\n";

            std::pair<int, double> key = {rec.k, rec.m};
            avgData[algo][key].sumPsnr   += rec.psnr;
            avgData[algo][key].sumSize   += rec.sizeCompressed;
            avgData[algo][key].sumRecall += rec.recall;
            avgData[algo][key].count++;
        }
        outfile.close();
        std::cout << "Fichier généré : " << outFileName << std::endl;
    }

    // Calcul et écriture des moyennes pour chaque algorithme.
    // Pour chaque algorithme, on parcourt toutes les paires (k, m) et on calcule la moyenne.
    for(auto& algoPair : avgData) {
        std::string algo = algoPair.first;
        std::string avgFileName = "all_data_values_sorted_" + algo + "/all_data_sorted_avg_" + algo;
        std::ofstream avgFile(avgFileName);
        if(!avgFile) {
            std::cerr << "Erreur lors de la création du fichier " << avgFileName << std::endl;
            continue;
        }

        // Pour obtenir un ordre cohérent, on transfère les données dans un vecteur que l'on trie par k (ascendant) puis par m (ascendant)
        std::vector<std::pair<std::pair<int,double>, Accumulator>> entries(algoPair.second.begin(), algoPair.second.end());
        std::sort(entries.begin(), entries.end(), [](const auto &a, const auto &b) {
            if(a.first.first == b.first.first)
                return a.first.second > b.first.second;
            return a.first.first < b.first.first;
        });

        // Écriture de la moyenne pour chaque paire (k, m)
        for(auto &entry : entries) {
            int k = entry.first.first;
            double m = entry.first.second;
            Accumulator acc = entry.second;
            double avgPsnr = acc.sumPsnr / acc.count;
            double avgSize = acc.sumSize / acc.count;
            double avgRecall = acc.sumRecall / acc.count;

            avgFile << k << " " << m << " "
                    << avgPsnr << " "
                    << avgSize << " "
                    << avgRecall << "\n";
        }
        avgFile.close();
        std::cout << "Fichier moyen généré : " << avgFileName << std::endl;
    }

    return 0;
}
