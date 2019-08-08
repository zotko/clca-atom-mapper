import functools
from collections import Counter
from itertools import chain
from operator import itemgetter, mul


def clca(molecule1, molecule2):

    # Create atom string for each atom
    orig_strings1 = [molecule1.atom_string(i) for i in range(len(molecule1))]
    orig_strings2 = [molecule2.atom_string(i) for i in range(len(molecule2))]

    # Set maximum prime number to 0
    max_prime = 0

    # Get adjacency list for each molecule
    adj_list1 = molecule1.adjacency_list
    adj_list2 = molecule2.adjacency_list

    # Start atom strings
    strings_and_primes1 = orig_strings1
    strings_and_primes2 = orig_strings2

    while True:
        # Get only strings from all values
        strings1 = [i for i in strings_and_primes1 if isinstance(i, str)]
        strings2 = [i for i in strings_and_primes2 if isinstance(i, str)]

        # Find singular atom strings
        string_counter1 = Counter(strings1)
        string_counter2 = Counter(strings2)
        string_counter_all = string_counter1 + string_counter2
        singular_strings = [string for string, number in string_counter_all.items() if number is 1]

        # Find atom strings unique in either graph, but common to both graphs
        unique_strings1 = set(str_ for str_, number in string_counter1.items() if number is 1)
        unique_strings2 = set(str_ for str_, number in string_counter2.items() if number is 1)
        unique_strings = unique_strings1 & unique_strings2
        print(strings1)
        print(strings2)
        print(unique_strings1)
        print(unique_strings2)
        print(unique_strings)

        prime_generator = gen_primes()

        # Check which prime numbers are already taken
        used_primes1 = set([i for i in strings_and_primes1 if isinstance(i, int) and i != 1])
        used_primes2 = set([i for i in strings_and_primes2 if isinstance(i, int) and i != 1])
        assert used_primes1 == used_primes2

        for prime in used_primes1: next(prime_generator) # Clean the generator from used prime numbers

        # Get prime numbers of unique atom strings
        primes_of_uniques = {string: next(prime_generator) for string in unique_strings}

        # Assign the value '1' to each singular atom string
        singular1 = (1 if i in singular_strings else i for i in strings_and_primes1)
        singular2 = (1 if i in singular_strings else i for i in strings_and_primes2)

        # Assign next prime number to atom strings unique in either graph, but common to both graphs
        singular_and_unique1 = [primes_of_uniques.get(i, i) for i in singular1]
        singular_and_unique2 = [primes_of_uniques.get(i, i) for i in singular2]

        # Get prime numbers of the rest of atom strings
        strings_rest = set([i for i in chain(singular_and_unique1, singular_and_unique2) if isinstance(i, str)])
        primes_of_rest = {str_: next(prime_generator) for str_ in strings_rest}

        # Assign next higher prime numbers to the rest of atom strings
        primes1 = [primes_of_rest.get(i, i) for i in singular_and_unique1]
        primes2 = [primes_of_rest.get(i, i) for i in singular_and_unique2]

        # Change on Github!!!
        # Calculate prime numbers product of adjacent atoms
        groups_of_neighbs1 = ([node] + list(neighbs) for node, neighbs in sorted(adj_list1.items()))
        groups_of_neighbs2 = ([node] + list(neighbs) for node, neighbs in sorted(adj_list2.items()))

        primes_of_neighbs1 = (itemgetter(*neighbs)(primes1) for neighbs in groups_of_neighbs1)
        primes_of_neighbs2 = (itemgetter(*neighbs)(primes2) for neighbs in groups_of_neighbs2)

        product_of_pimes1 = [functools.reduce(mul, primes) for primes in primes_of_neighbs1]
        product_of_pimes2 = [functools.reduce(mul, primes) for primes in primes_of_neighbs2]

        primes_of_mapped1 = [i if isinstance(i, int) else None for i in singular_and_unique1]
        primes_of_mapped2 = [i if isinstance(i, int) else None for i in singular_and_unique2]
        print(primes_of_mapped1)
        print(primes_of_mapped2)
        print()

        # Generate new atom strings by concatenation of original atom strings with prime numbers products
        strings_and_primes1 = [string + str(product) if not prime else prime
                               for string, product, prime in
                               zip(orig_strings1, product_of_pimes1, primes_of_mapped1)]
        strings_and_primes2 = [string + str(product) if not prime else prime
                               for string, product, prime in
                               zip(orig_strings2, product_of_pimes2, primes_of_mapped2)]
                                   
        # Evaluate new maximum prime number
        new_max_prime = max(primes1 + primes2, default=0)        

        # If the highest assigned prime number has not changed,
        # no new information can be retrieved and the algorithm stops
        if new_max_prime == max_prime:
            break
        else:
            max_prime = new_max_prime


    # Add information about mapped atoms
    primes1 = [(prime, bool(mapped_prime)) for prime, mapped_prime in zip(primes1, primes_of_mapped1)]
    primes2 = [(prime, bool(mapped_prime)) for prime, mapped_prime in zip(primes2, primes_of_mapped2)]

    return primes1, primes2
