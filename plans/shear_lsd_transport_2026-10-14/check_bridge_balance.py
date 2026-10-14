"""Enumerate bridge paths and verify detailed balance on a finite target."""
import itertools
import math


def check(lengths):
    states = list(itertools.product((0, 1), repeat=2))
    prior = (0.35, 0.65)
    factors = (((0.08, 0.92), (0.75, 0.25)),
               ((0.86, 0.14), (0.18, 0.82)))

    def logf(group, subject, value):
        return math.log(factors[group][subject][value])

    def paths(old, new, subject, start, length):
        for tail in itertools.product((0, 1), repeat=length - 1):
            path = (start,) + tail
            probability = 1.0
            for i in range(1, length):
                beta = i / length
                a, b = path[i - 1:i + 1]

                def density(x):
                    return ((1 - beta) * logf(old, subject, x)
                            + beta * logf(new, subject, x))

                flip = 0.6 * min(1, math.exp(density(1 - a) - density(a)))
                probability *= flip if b != a else 1 - flip
            work = sum(logf(new, subject, x) - logf(old, subject, x)
                       for x in path) / length
            yield path[-1], probability, work

    joint = list(itertools.product((0, 1), states))
    transition = [[0.0 for _ in joint] for _ in joint]
    target = [prior[g] * math.prod(factors[g][s][x[s]] for s in range(2))
              for g, x in joint]
    for row, (old, x) in enumerate(joint):
        new = 1 - old
        alternatives = [list(paths(old, new, s, x[s], lengths[s])) for s in range(2)]
        for candidate in itertools.product(*alternatives):
            y = tuple(c[0] for c in candidate)
            probability = math.prod(c[1] for c in candidate)
            work = math.log(prior[new] / prior[old]) + sum(c[2] for c in candidate)
            accepted = probability * min(1, math.exp(work))
            transition[row][joint.index((new, y))] += accepted
            transition[row][row] += probability - accepted
    error = max(abs(target[i] * transition[i][j] - target[j] * transition[j][i])
                for i in range(len(joint)) for j in range(len(joint)))
    row_error = max(abs(sum(row) - 1) for row in transition)
    assert error < 1e-13, error
    assert row_error < 1e-13, row_error
    print(f"K={lengths}: detailed-balance error {error:.3g}; row-sum error {row_error:.3g}")


if __name__ == "__main__":
    for lengths in ((1, 1), (2, 2), (4, 4), (1, 4)):
        check(lengths)
