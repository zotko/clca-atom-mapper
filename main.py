"""CLCA atom mapper.

This script uses the Canonical Labeling for Clique Approximation (CLCA) algorithm to identify atom-atom mappings between two molecules.

This script accepts two XYZ ﬁles with chemical elements and their cartesian coordinates as input.

The `plotly` package is required for visualisation of molecular graphs.
"""
from itertools import chain
from sys import argv

from clca import clca
from molecular_graph import Graph
from plot import plot, gen_trace

# The algorithm takes two molecular graphs as input.
script, file1, file2 = argv

molecule1, molecule2 = Graph(), Graph()

molecule1.read_file(file1)
molecule2.read_file(file2)

primes1, primes2 = clca(molecule1, molecule2)

# Convert prime numbers to natural numbers and set all values '1' to 's' (singular)
all_primes = sorted(list(set(prime for prime, mapped in chain(primes1, primes2))))
if 1 in all_primes: all_primes.remove(1)
numbers1 = [(all_primes.index(prime), mapped) if prime is not 1 else ('s', False) for prime, mapped in primes1]
numbers2 = [(all_primes.index(prime), mapped) if prime is not 1 else ('s', False) for prime, mapped in primes2]

# Identify only unambiguous atom-atom mappings
mapping = {i: primes2.index((prime, mapped)) for i, (prime, mapped) in enumerate(primes1) if mapped}
print('Mapping:')
print(mapping)

labels_and_color1 = [(label, 'red') if mapped else (label, 'black') for label, mapped in numbers1]
labels_and_color2 = [(label, 'red') if mapped else (label, 'black') for label, mapped in numbers2]

# Generate atom and bond traces for the plot
trace1 = gen_trace(adjacency_list=molecule1.adjacency_list, elements=molecule1.elements,
                   x_coordinates=molecule1.x_coordinates, y_coordinates=molecule1.y_coordinates,
                   z_coordinates=molecule1.z_coordinates)

trace2 = gen_trace(adjacency_list=molecule2.adjacency_list, elements=molecule2.elements,
                   x_coordinates=molecule2.x_coordinates, y_coordinates=molecule2.y_coordinates,
                   z_coordinates=molecule2.z_coordinates)

# Generate equal labels for equal atoms for the plot
element_labels1 = [dict(text=f'{el} ({label})', x=x, y=y, z=z, showarrow=False, yshift=15, font=dict(color=color))
                   for (label, color), (el, (x, y, z)) in zip(labels_and_color1, molecule1)]
element_labels2 = [dict(text=f'{el} ({label})', x=x, y=y, z=z, showarrow=False, yshift=15, font=dict(color=color))
                   for (label, color), (el, (x, y, z)) in zip(labels_and_color2, molecule2)]

plot(trace1, trace2, annotations1=element_labels1, annotations2=element_labels2)