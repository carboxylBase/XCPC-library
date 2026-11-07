"""Compile solution.cpp and run a local interactor; no third-party packages.

Example: CXX=g++-13 python3 check.py
Each input batch obeys T <= 250 and sum(n) <= 500.
"""

import argparse
from functools import lru_cache
import itertools
import os
from pathlib import Path
import random
import select
import shlex
import subprocess
import tempfile


def check_batch(executable, cases):
    process = subprocess.Popen(
        [str(executable)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, bufsize=0,
    )
    pending = bytearray()
    maxima = {}

    def send(value):
        process.stdin.write((str(value) + "\n").encode())

    def line():
        while b"\n" not in pending:
            ready, _, _ = select.select([process.stdout], [], [], 5)
            assert ready, "Timed out waiting for flushed output"
            chunk = os.read(process.stdout.fileno(), 4096)
            assert chunk, "Program exited before giving an answer"
            pending.extend(chunk)
        end = pending.index(b"\n")
        result = bytes(pending[:end]).decode().split()
        del pending[:end + 1]
        return result

    try:
        send(len(cases))
        for values in cases:
            n = len(values)
            send(n)
            count = 0
            while True:
                fields = line()
                assert fields, "Empty output line"
                if fields[0] == "?":
                    assert len(fields) == 3, fields
                    i, j = map(int, fields[1:])
                    assert 1 <= i <= n and 1 <= j <= n and i != j, fields
                    count += 1
                    assert count <= 2 * n - 3, (values, count)
                    send(int(values[i - 1] > values[j - 1]))
                else:
                    assert fields[0] == "!", fields
                    order = [int(v) - 1 for v in fields[1:]]
                    assert sorted(order) == list(range(n)), (values, order)
                    assert all(values[order[i]] <= values[order[i + 1]]
                               for i in range(n - 1)), (values, order)
                    maxima[n] = max(maxima.get(n, 0), count)
                    break
        process.stdin.close()
        assert process.wait(timeout=5) == 0, process.stderr.read().decode()
        assert not pending.strip() and not process.stdout.read().strip(), "Extra output"
        return maxima
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
        for stream in (process.stdin, process.stdout, process.stderr):
            if not stream.closed:
                stream.close()


def check_termination(executable, reply):
    """After the first query, -1 or EOF must terminate without more output."""
    result = subprocess.run(
        [str(executable)], input=b"1\n2\n" + reply,
        capture_output=True, timeout=5,
    )
    assert result.returncode == 0, result.stderr.decode()
    assert result.stdout.decode().splitlines() == ["? 1 2"], result.stdout


def audit_rejected_strategy():
    """Search every continuation of the n=4 transcript discussed in notes.md."""
    arrays = list(itertools.product(range(3), repeat=4))
    orders = list(itertools.permutations(range(4)))
    valid = []
    for values in arrays:
        mask = 0
        for index, order in enumerate(orders):
            if all(values[order[i]] <= values[order[i + 1]] for i in range(3)):
                mask |= 1 << index
        valid.append(mask)
    state = tuple(k for k, a in enumerate(arrays)
                  if a[0] > a[1] and a[2] <= a[0] and a[3] <= a[0])

    @lru_cache(None)
    def can_finish(possible, depth):
        common = (1 << len(orders)) - 1
        for index in possible:
            common &= valid[index]
        if common:
            return True
        if depth == 0:
            return False
        for i in range(4):
            for j in range(4):
                if i == j:
                    continue
                yes = tuple(k for k in possible if arrays[k][i] > arrays[k][j])
                no = tuple(k for k in possible if arrays[k][i] <= arrays[k][j])
                if yes and no and can_finish(yes, depth - 1) and can_finish(no, depth - 1):
                    return True
        return False

    assert len(state) == 22
    assert not can_finish(state, 2)
    assert can_finish(state, 3)
    print("Rejected scan strategy: 22 consistent arrays; 2 more queries fail, 3 suffice.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--random-cases", type=int, default=100)
    args = parser.parse_args()
    assert 2 <= args.max_n <= 11 and args.random_cases >= 0
    rng = random.Random(20261027)
    cases = [a for n in range(2, args.max_n + 1)
             for a in itertools.product(range(3), repeat=n)]
    cases += [tuple([value] * 500) for value in range(3)]
    cases += [tuple(sorted(i % 3 for i in range(500))),
              tuple(sorted((i % 3 for i in range(500)), reverse=True)),
              tuple(i % 3 for i in range(500)),
              tuple(2 - i % 3 for i in range(500)),
              tuple([2] * 499 + [0]), tuple([0] + [2] * 499)]
    cases += [tuple(rng.randrange(3) for _ in range(500))
              for _ in range(args.random_cases)]
    maxima = {}
    batches = 0
    with tempfile.TemporaryDirectory(prefix="hidden-order-check-") as folder:
        executable = Path(folder) / "solution"
        source = Path(__file__).resolve().with_name("solution.cpp")
        subprocess.run(shlex.split(args.cxx) + ["-std=c++17", "-O2", str(source),
                       "-o", str(executable)], check=True)
        batch, size = [], 0
        for values in cases + [None]:
            if batch and (values is None or len(batch) == 250 or size + len(values) > 500):
                for n, count in check_batch(executable, batch).items():
                    maxima[n] = max(maxima.get(n, 0), count)
                batches += 1
                batch, size = [], 0
            if values is not None:
                batch.append(values)
                size += len(values)
        check_termination(executable, b"-1\n")
        check_termination(executable, b"")
    print(f"C++17 compiled; {len(cases)} cases passed in {batches} valid input batches.")
    print("Max queries by n:", maxima)
    print("Read failure and -1 termination checks passed.")
    audit_rejected_strategy()


if __name__ == "__main__":
    main()
