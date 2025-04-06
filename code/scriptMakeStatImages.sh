#!/bin/bash

# Définition des valeurs de k et m
k_values=(100 500 2500 5000 10000 15000 20000 30000)
m_values=(0.5 1.0 2.0 4.0 8.0 10.0)

# Création du dossier output s'il n'existe pas
mkdir -p output_stat

# Préparer une liste de commandes à exécuter
commands_file="commands.txt"
> "$commands_file"  # vide le fichier s'il existe

for image in images_test_orginal/*.ppm; do
    number=$(basename "$image" .ppm)
    mkdir -p "output_stat/$number"
    for k in "${k_values[@]}"; do
        for m in "${m_values[@]}"; do
            output_file="output_stat/${number}/${number}_${k}_${m}.ppm"
            if [ ! -f "output_stat/${number}/${number}_${k}_${m}SLIC.ppm" ]; then
                echo "./main images_test_orginal/${number}.ppm $output_file $k $m contour_humain/${number}.pgm" >> "$commands_file"
            fi
        done
    done
done

# Exécuter les commandes en parallèle (par exemple avec 4 jobs en parallèle)
parallel -j 10 < "$commands_file"
echo "Traitement terminé."
