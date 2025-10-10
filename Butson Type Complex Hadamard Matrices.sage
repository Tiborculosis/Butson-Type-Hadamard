#type: ignore

from itertools import product
import multiprocessing as mp
from sage.all import Integer
from functools import lru_cache
import time

@lru_cache(maxsize=None)  # unlimited cache
def cached_permutations(tail_tuple):
    """Return all permutations of a tail (tuple), caching results."""
    return list(Permutations(tail_tuple))

@lru_cache(maxsize=None)
def rows_orthogonal_cached(u_tuple, v_tuple):
    # Use global zeta
    return sum(zeta**(u_i - v_i) for u_i, v_i in zip(u_tuple, v_tuple)) == 0

def cycle_exponents(k, d, shift=0):
	"""
	Generate a cycle of length k of dth roots of unity, where k is a prime divisor of d, rotated by shift.
	"""
	step = d // k
	return [(shift + j * step) % d for j in range(k)]

def partitions_with_parts(n, parts):
	"""
	Generate all integer combinations of given parts that sum to n.

	Args:
	- n (int): the total to reach (L1 norm)
	- parts (list[int]): allowed part sizes, e.g. [p, q]

	Returns:
	- list of tuples (a1, a2, ..., ak) where sum(ai * parts[i]) == n
	"""
	results = []
	def helper(idx, remaining, coeffs):
		"""
		Generate the actual integer combinations of parts recursively. Currently only valid for up to 2 primes.
		
		Args:
		- idx (int): the index of which part is actually being incremented
		- remaining (int): the remainder; it tells us when we've got a valid combination or have missed our window to find one
		- coeffs (list[tuple[int]]): coefficients of the parts from the parent function
		"""
		if idx == len(parts) - 1:
			if remaining % parts[idx] == 0:
				results.append(tuple(coeffs + [remaining // parts[idx]]))
			return
		for count in range(remaining // parts[idx] + 1):
			helper(idx + 1, remaining - count * parts[idx], coeffs + [count])
	helper(0, n, [])
	return results

def generate_L1_vectors(n, d, require_zero=True):
	"""
	Generate a list of vectors of dth roots of unity which have L1 norm n.

	Args:
	- n (int): The L1 norm of the vectors
	- d (int): The primitive root of unity out of which these are constructed
	- require_zero (boolean): If true, forces the vectors to contain at least one 1 (exponent 0) to account for dephased matrix

	Returns:
	- List of all vectors satisfying the above conditions
	"""

	# Make sure n and d are integers
	if not isinstance(n, (int, Integer)):
		raise TypeError(f"n must be an integer, not {type(n).__name__}")
	if not isinstance(d, (int, Integer)):
		raise TypeError(f"d must be an integer, not {type(d).__name__}")
    
	# Value checks
	if n < 0:
		raise ValueError(f"n must be nonnegative (L1 norm), got {n}")
	if d < 2:
		raise ValueError(f"d must be a natural number ≥ 2, got {d}")

	primes = prime_divisors(d)

	if len(primes) > 2:
		raise ValueError(f"d must have at most 2 prime factors (at least for now), got {len(primes)}")

	partitions = partitions_with_parts(n, primes)
	all_vectors = set()

	for coeffs in partitions:
		prime_counts = list(zip(primes, coeffs))
		# Build all combinations of shifts for all cycles
		shift_space = []
		for k, count in prime_counts:
			# each cycle of size k can be rotated independently
			shift_space.append(product(range(d // k), repeat=count))
		# Cartesian product over all primes’ shift sets
		for combo in product(*shift_space):
			vector = []
			for (k, count), shifts in zip(prime_counts, combo):
				for shift in shifts:
					vector.extend(cycle_exponents(k, d, shift))
			vector.sort()
			if require_zero and 0 not in vector:
				continue
			all_vectors.add(tuple(vector))

	return [list(v) for v in sorted(all_vectors)]

def valid_permutation(candidate_row, existing_rows):
	"""
	Return True if candidate_row respects the "don't swap identical positions in previous rows" rule.
	"""
	n = len(candidate_row)
	for i in range(n):
		for j in range(i+1, n):
			if candidate_row[i] > candidate_row[j]:
				# Check if positions i and j are identical in all previous rows
				if all(row[i] == row[j] for row in existing_rows):
					return False
	return True

def build_permutation_core_CHMs(core_coeffs, n, d):
	"""
	core_coeffs: list of integers
	n: size of the CHM
	d: dth root of unity, i.e. e^{2pi i/d}
	Returns: list of matrices (each matrix is a list of rows, each row is a list)
	"""
	first_row = [0]*n
	second_row = core_coeffs
	fixed_first = core_coeffs[0]
	tail = core_coeffs[1:]
	results = []

	def backtrack(existing_rows):
		if len(existing_rows) == n:
	    	# Deep copy of matrix
			results.append([row[:] for row in existing_rows])
			return

    	# Generate only unique permutations of tail
		for perm in cached_permutations(tuple(tail)):
			candidate_row = [fixed_first] + list(perm)

			# Canonical ordering: only consider rows >= last row
			if candidate_row > existing_rows[-1] and \
				all(rows_orthogonal_cached(tuple(candidate_row), tuple(row)) for row in existing_rows) and \
					valid_permutation(candidate_row, existing_rows):
				existing_rows.append(candidate_row)
				backtrack(existing_rows)
				existing_rows.pop()

	backtrack([first_row, second_row])
	return results

def write_chms_to_file(chms, filename="output.txt"):
	with open(f"output/{filename}", "w") as f:
		if chms:
			for i, mat in enumerate(chms, 1):
				f.write(f"Matrix {i}:\n")
				for row in mat:
					f.write(" ".join(map(str, row)) + "\n")
				f.write("\n")
		else:
			f.write(f"No matrices found.")

# --- Helper for chunking ---
def chunk_list(lst, chunk_size):
    """Split a list into smaller chunks for better CPU utilization."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i+chunk_size]

# --- Top-level function for permutation-core CHMs ---
def process_chunk(chunk, n, d):
    """
    Compute permutation-core CHMs for a chunk of rows.
    Returns a list of matrices.
    """
    results = []
    for row in chunk:
        results.extend(build_permutation_core_CHMs(row, int(n), int(d)))
    return results

# --- Top-level function for Butson CHMs ---
def process_butson(_unused, vector_list, n, d):
    """
    Compute Butson CHMs.
    _unused is just a dummy argument for starmap compatibility.
    """
    return build_Butson_CHMs(vector_list, int(n), int(d))

# --- Main parallel runner ---
def run_parallel(n, d):
    # Ensure inputs are Sage Integers
    n = Integer(n)
    d = Integer(d)
    t0 = time.time()

    # Generate the L1 vectors (single-threaded)
	
    vector_list = generate_L1_vectors(n, d)

    # Convert all numbers to plain Python ints to reduce pickling overhead
    vector_list = [[int(x) for x in row] for row in vector_list]
	
    vector_list.sort() # Minimize canonical checks

    t1 = time.time()
    print(f"L1 vector generation ({n}, {d}): {t1 - t0:.3f} seconds")

    # Determine chunk size: more tasks than cores
    cpu_count = mp.cpu_count()
    chunk_size = max(1, len(vector_list) // (cpu_count * 4))  # ~4 tasks per core
    chunks = list(chunk_list(vector_list, chunk_size))

    # Use spawn context for WSL/Linux
    ctx = mp.get_context('spawn')
    t2 = time.time()

    # --- Parallel permutation-core CHMs ---
    permutation_chms = []
    with ctx.Pool(cpu_count, initializer=init_worker, initargs=(d,)) as pool:
        for res in pool.starmap(process_chunk, [(chunk, n, d) for chunk in chunks]):
            permutation_chms.extend(res)

    t3 = time.time()
    print(f"Permutation-core CHMs ({n}, {d}): {t3 - t2:.3f} seconds")
    # Write results
    write_chms_to_file(permutation_chms, f"Permutation-core H({int(n)}, {int(d)}).txt")

    # --- Parallel Butson CHMs ---
    safe_vector_list = [tuple(row) for row in vector_list]  # tuples are pickle-safe
    t4 = time.time()
    butson_chms = []

    with ctx.Pool(cpu_count, initializer=init_worker, initargs=(d,)) as pool:
        # map each starting row chunk to a process
        for res in pool.starmap(process_butson, [(chunk, vector_list, n, d) for chunk in chunks]):
            butson_chms.extend(res)

    t5 = time.time()
    print(f"Butson-type CHMs ({n}, {d}): {t5 - t4:.3f} seconds")
    # write to file
    write_chms_to_file(butson_chms, f"Butson-type H({int(n)}, {int(d)}).txt")

    print(f"Total runtime ({n}, {d}): {t5 - t0: .3f} seconds")

# --- Top-level function for Butson CHMs with a starting row ---
def build_Butson_CHMs(candidate_rows, n, d, start_row):
    """
    Generate all n x n Butson-type CHMs whose rows are drawn from candidate_rows,
    starting with start_row. Comparisons still consider all rows.
    """
    results = []
    first_row = [0]*n  # standard dephased first row
	
    def backtrack(existing_rows, recent_row):
        if len(existing_rows) == n:
            results.append([row[:] for row in existing_rows])
            return

        for base_row in candidate_rows:
            if base_row < recent_row:
                continue
			
            fixed_first = base_row[0]
            tail = base_row[1:]
			
            for perm in cached_permutations(tuple(tail)):
                row = [fixed_first] + list(perm)
				
                if base_row == recent_row and row < existing_rows[-1]:
                    continue
				
                # skip duplicate row
                if row in existing_rows:
                    continue

                # ensure orthogonality with all previous rows
                if not all(rows_orthogonal_cached(tuple(row), tuple(prev)) for prev in existing_rows):
                    continue

                # symmetry breaking: don't swap identical columns of previous rows
                if not valid_permutation(row, existing_rows):
                    continue

                # canonical ordering: because candidate_rows is pre-sorted, 
                # we can skip this check
                # if row < existing_rows[-1]:
                #     continue
				
                recent_row = base_row
                existing_rows.append(row)
                backtrack(existing_rows, recent_row)
                existing_rows.pop()

    backtrack([first_row, start_row], start_row)
    return results

# --- Worker function for multiprocessing ---
def process_butson(start_rows_chunk, full_vector_list, n, d):
    """
    Each process takes a chunk of starting rows, builds all CHMs starting from them.
    """
    results = []
    for start_row in start_rows_chunk:
        results.extend(build_Butson_CHMs(full_vector_list, n, d, start_row))
    return results

def init_worker(d_value):
    global zeta
    from sage.all import CyclotomicField
    zeta = CyclotomicField(d_value).gen()

if __name__ == "__main__":
    start_time = time.time()
    for d in range(2, 13):
        R = CyclotomicField(d)
        zeta = R.gen()
        rows_orthogonal_cached.cache_clear() # Avoid conflicts that may arise due to using the same tuples with different zeta values
        for n in range(2, 12):
            run_parallel(n, d)