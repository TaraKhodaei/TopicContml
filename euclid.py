#!/usr/bin/env python
# calculates distance matrix from infile generated in topiccontml
# when using the option --neighbor this will then be used in the neighbor
# program 
import sys
import numpy as np

def calculate_euclidean_distance_matrix(frequencies):
    num_species = len(frequencies)
    distance_matrix = np.zeros((num_species, num_species))

    for i in range(num_species):
        for j in range(i + 1, num_species):
            distance = np.linalg.norm(np.array(frequencies[i]) - np.array(frequencies[j]))
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance

    return distance_matrix

def print_distance_matrix(species_names, distance_matrix, filename):
    num_species = len(species_names)
    with open(filename,'w') as f:
        f.write(f'{num_species}\n')
        for i in range(num_species):
            f.write(f"{species_names[i]:<15}")
            for j in range(num_species):
                if j == num_species - 1:
                    f.write(f"{distance_matrix[i][j]:.6f}")
                else:
                    f.write(f"{distance_matrix[i][j]:.6f} ")
            f.write('\n')

def distances(infile,outfile):
    #input_data = sys.stdin.read().splitlines()
    with open(infile,'r') as f:
        input_data = f.readlines()
    num_species, num_loci = map(int, input_data[0].split())
    alleles_per_locus = list(map(int, input_data[1].split()))

    species_names = []
    frequencies = []

    for line in input_data[2:]:
        parts = line.split()
        species_names.append(parts[0])
        freqs = list(map(float, parts[1:]))
        
        # Include the additional frequency for each locus
        extended_freqs = []
        start = 0
        for num_alleles in alleles_per_locus:
            locus_freqs = freqs[start:start + num_alleles]
            extended_freqs.extend(locus_freqs)
            extended_freqs.append(1.0 - sum(locus_freqs))
            start += num_alleles
        
        frequencies.append(extended_freqs)

    distance_matrix = calculate_euclidean_distance_matrix(frequencies)
    print_distance_matrix(species_names, distance_matrix, outfile)

if __name__ == "__main__":
    # give infile filename as argument
    distances(sys.argv[1])
