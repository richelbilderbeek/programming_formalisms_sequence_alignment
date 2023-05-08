def add_one(x):
    return x + 1

from collections import deque

def alignment_score(x, y, alignment):
    """Score an alignment.

    x, y -- sequences.
    alignment -- an alignment of x and y.
    """
    score_gap = -1
    score_same = +1
    score_different = -1

    score = 0
    for i, j in alignment:
        if (i is None) or (j is None):
            score += score_gap
        elif x[i] == y[j]:
            score += score_same
        elif x[i] != y[j]:
            score += score_different

    return score

def all_alignments(x, y):
    """Return an iterable of all alignments of two
    sequences.

    x, y -- Sequences.
    """

    def F(x, y):
        """A helper function that recursively builds the
        alignments.

        x, y -- Sequence indices for the original x and y.
        """
        if len(x) == 0 and len(y) == 0:
            yield deque()

        scenarios = []
        if len(x) > 0 and len(y) > 0:
            scenarios.append((x[0], x[1:], y[0], y[1:]))
        if len(x) > 0:
            scenarios.append((x[0], x[1:], None, y))
        if len(y) > 0:
            scenarios.append((None, x, y[0], y[1:]))

        # NOTE: "xh" and "xt" stand for "x-head" and "x-tail",
        # with "head" being the front of the sequence, and
        # "tail" being the rest of the sequence. Similarly for
        # "yh" and "yt".
        for xh, xt, yh, yt in scenarios:
            for alignment in F(xt, yt):
                alignment.appendleft((xh, yh))
                yield alignment

    alignments = F(range(len(x)), range(len(y)))
    return map(list, alignments)
