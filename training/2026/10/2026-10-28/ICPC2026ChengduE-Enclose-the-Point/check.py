"""Compile the archived code and compare its output with exact brute force.

Run: CXX=g++-13 python3 check.py
The generated input obeys the statement's total n, total q and T limits.
"""

import argparse
import itertools
import os
from pathlib import Path
import random
import shlex
import subprocess
import tempfile


def cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def hull(points):
    points = sorted(points)
    lower, upper = [], []
    for point in points:
        while len(lower) > 1 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    for point in reversed(points):
        while len(upper) > 1 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def in_convex(polygon, point):
    return len(polygon) >= 3 and all(
        cross(polygon[i], polygon[(i + 1) % len(polygon)], point) > 0
        for i in range(len(polygon))
    )


def in_triangle(triangle, point):
    values = [cross(triangle[i], triangle[(i + 1) % 3], point) for i in range(3)]
    return all(v > 0 for v in values) or all(v < 0 for v in values)


def on_segment(a, b, point):
    return cross(a, b, point) == 0 and (
        min(a[0], b[0]) <= point[0] <= max(a[0], b[0]) and
        min(a[1], b[1]) <= point[1] <= max(a[1], b[1])
    )


def intersects(a, b, c, d):
    if any((on_segment(a, b, p) for p in (c, d))):
        return True
    if any((on_segment(c, d, p) for p in (a, b))):
        return True
    return (cross(a, b, c) * cross(a, b, d) < 0 and
            cross(c, d, a) * cross(c, d, b) < 0)


def validate_polygon(polygon, point):
    n = len(polygon)
    inside = False
    for i in range(n):
        a, b = polygon[i], polygon[(i + 1) % n]
        assert not on_segment(a, b, point), (polygon, point)
        if (a[1] > point[1]) != (b[1] > point[1]):
            if (cross(a, b, point) > 0) == (b[1] > a[1]):
                inside = not inside
        for j in range(i + 1, n):
            if j == i + 1 or (i == 0 and j == n - 1):
                continue
            assert not intersects(a, b, polygon[j], polygon[(j + 1) % n]), polygon
    assert inside, (polygon, point)


def optimum(points, convex, point):
    if not in_convex(convex, point):
        return -1
    if any(in_triangle(t, point) for t in itertools.combinations(points, 3)):
        return 3
    # Independently find a four-point witness rather than assuming it exists.
    assert any(in_convex(hull(t), point) for t in itertools.combinations(points, 4))
    return 4


def make_cases():
    grid = list(itertools.product(range(-2, 3, 2), repeat=2))
    queries = list(itertools.product(range(-2, 3), repeat=2))
    cases = [(list(a), queries) for n in range(3, 10)
             for a in itertools.combinations(grid, n)]
    rng = random.Random(20261028)
    pool = list(itertools.product(range(-8, 9, 2), repeat=2))
    for _ in range(1000):
        points = rng.sample(pool, rng.randrange(3, 13))
        qs = [(rng.randrange(-9, 10), rng.randrange(-9, 10)) for _ in range(10)]
        cases.append((points, qs))
    limit = 10 ** 9
    cases.append((
        [(-limit, -limit), (limit, -limit), (limit, limit), (-limit, limit), (0, 1)],
        [(0, 0), (0, 1), (limit, 0), (limit, limit), (-limit + 1, -limit + 1)],
    ))
    square = [(-2, -2), (2, -2), (2, 2), (-2, 2)]
    cases.append((square + [(0, 1)], [(0, 0)]))
    cases.append((square + [(1, 1), (-1, -1), (0, 0)], [(0, 0)]))
    # The valid witness has index >= 4 after sorting a; indexing b is invalid.
    cases.append((
        [(-10, -10), (10, -10), (10, 10), (-10, 10),
         (-9, -9), (-8, -8), (-7, -7), (0, 1)],
        [(0, 0), (0, 1), (-9, -9)],
    ))
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    args = parser.parse_args()
    cases = make_cases()
    assert len(cases) <= 20000
    assert sum(len(a) for a, _ in cases) <= 100000
    assert sum(len(q) for _, q in cases) <= 100000
    data = [str(len(cases))]
    for points, queries in cases:
        data.append(f"{len(points)} {len(queries)}")
        data.extend(f"{x} {y}" for x, y in points)
        data.extend(f"{x} {y}" for x, y in queries)
    with tempfile.TemporaryDirectory(prefix="enclose-check-") as folder:
        executable = Path(folder) / "solution"
        source = Path(__file__).resolve().with_name("solution.cpp")
        subprocess.run(shlex.split(args.cxx) + ["-std=c++17", "-O2", str(source),
                       "-o", str(executable)], check=True)
        result = subprocess.run([str(executable)], input="\n".join(data) + "\n",
                                text=True, capture_output=True, check=True, timeout=30)
    tokens = iter(result.stdout.split())
    counts = {-1: 0, 3: 0, 4: 0}
    for points, queries in cases:
        convex = hull(points)
        for point in queries:
            expected = optimum(points, convex, point)
            answer = int(next(tokens))
            assert answer == expected, (points, point, answer, expected)
            counts[answer] += 1
            if answer != -1:
                indices = [int(next(tokens)) - 1 for _ in range(answer)]
                assert len(set(indices)) == answer, indices
                assert all(0 <= i < len(points) for i in indices), indices
                validate_polygon([points[i] for i in indices], point)
    assert next(tokens, None) is None, "Unexpected extra output"
    print("C++17 compiled; point sets:", len(cases))
    print("Queries checked:", sum(counts.values()), "answer counts:", counts)
    print("All indices, simple polygons, strict containment and minima passed.")


if __name__ == "__main__":
    main()
