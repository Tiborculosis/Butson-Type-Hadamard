import time

def build_row(bit_index, ones_left, current_row, prev_rows, overlaps_needed, m, n):
    """
    Incrementally build a row in minimal lexicographic order,
    pruning branches that cannot satisfy n-1 overlaps with previous rows.
    """
    if bit_index == m:
        if ones_left == 0 and all(x == 0 for x in overlaps_needed):
            yield current_row
        return

    remaining_positions = m - bit_index

    # Try placing 1 at this position
    if ones_left > 0:
        new_overlaps = []
        prune = False
        for prev, needed in zip(prev_rows, overlaps_needed):
            if (prev >> (m - 1 - bit_index)) & 1:
                needed -= 1
            if needed < 0 or needed > remaining_positions - 1:
                prune = True
                break
            new_overlaps.append(needed)
        if not prune:
            yield from build_row(
                bit_index + 1, ones_left - 1,
                current_row | (1 << (m - 1 - bit_index)),
                prev_rows, new_overlaps, m, n
            )

    # Try placing 0 at this position
    new_overlaps = []
    prune = False
    for prev, needed in zip(prev_rows, overlaps_needed):
        if (prev >> (m - 1 - bit_index)) & 1:
            if needed > remaining_positions - 1:
                prune = True
                break
        new_overlaps.append(needed)
    if not prune:
        yield from build_row(bit_index + 1, ones_left, current_row, prev_rows, new_overlaps, m, n)

def find_minimal_hadamard(n):
    size = 4 * n
    total_ones = 2 * n - 1
    m = size - 1  # exclude first column (normalized)

    # first row: leftmost ones
    first_row = sum(1 << (m - 1 - i) for i in range(total_ones))
    matrix = [first_row]

    start_time = time.perf_counter()

    def backtrack(rows):
        if len(rows) == size - 1:
            return rows
        overlaps_needed = [n - 1] * len(rows)
        for candidate in build_row(0, total_ones, 0, rows, overlaps_needed, m, n):
            rows.append(candidate)
            result = backtrack(rows)
            if result:
                return result
            rows.pop()
        return None

    result = backtrack(matrix)
    end_time = time.perf_counter()
    print(f"Elapsed time: {end_time - start_time:.6f} seconds")
    return result

def print_hadamard_01(n, rows_bitmask):
    size = 4 * n
    m = size - 1
    print("Row 0:", [1]*size)
    for i, row in enumerate(rows_bitmask, start=1):
        row_list = [1] + [1 if (row >> (m - 1 - j)) & 1 else 0 for j in range(m)]
        print(f"Row {i}:", row_list)

if __name__ == "__main__":
    n = 4  # adjust as needed
    rows_bitmask = find_minimal_hadamard(n)
    if rows_bitmask:
        print_hadamard_01(n, rows_bitmask)