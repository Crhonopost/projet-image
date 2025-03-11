#!/bin/bash

# Définition des valeurs de k et m
k_values=(100 500 2500 5000 10000)
m_values=(0.5 1.0 2.0 4.0 8.0 10.0)

# Création du dossier output s'il n'existe pas
mkdir -p output

# Fonction pour traiter une image
process_image() {
    image="$1"
    filename=$(basename -- "$image")
    name="${filename%.*}"
    mkdir -p "output/$name"

    for k in "${k_values[@]}"; do
        for m in "${m_values[@]}"; do
            mkdir -p "output/$name/${name}_${k}_${m}"
            output_path="output/$name/${name}_${k}_${m}/${name}_${k}_${m}.ppm"
            if [ ! -f "$output_path" ]; then
                echo "Processing $image with k=$k and m=$m -> $output_path"
                ./main "$image" "$output_path" "$k" "$m"
            else
                echo "Skipping $output_path, already exists."
            fi
        done
    done
}

export -f process_image

# Trouver tous les fichiers .ppm et les traiter en parallèle
find images -name '*.ppm' | while read -r img; do
    process_image "$img" &
done

# Attendre que tous les processus en arrière-plan terminent
wait
