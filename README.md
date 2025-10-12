## Complex
This portion thus far explores the generation of Butson-type complex Hadamard matrices, along with their more specific relatives, permutation-core Butson-type complex Hadamard matrices.

The `.sage` program is the current iteration of my work, and the `.ipynb` is severely out of date.
The sage program is a bit harder to read, since a good portion of the logic is parallelized and therefore a bit less straightforward to read.

### Output

This contains the outputs of the program. It works essentially as intended, but generates more matrices than are given by the aalto.fi database. For example, butson-6-6.txt contains 6 CHMs, but the database suggests there should be only four, up to monomial equivalence. This is, however, a drastic improvement over the previous 20 that were generated before significant pruning measures were implemented.

## Real
At present, this is more of a shot-in-the-dark hobby project than anything deserving substantial recognition. `minimal-real-generator.py` is, in theory, an attempt to generate the lexicographically-minimal real Hadamard matrix in a given dimension.

The logic seems to be sound, but the combinatorial explosion seems to occur for a 4n x 4n real Hadamard matrix for n values greater than or equal to 5. Until and unless I come up with some flash of inspiration to radically optimize the program, it's essentially a moot project.