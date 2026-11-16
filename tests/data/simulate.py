"""Simulate a 101 bp TruSeq paired-end isolate library contaminated with adapter dimers.

Usage: simulate.py N_GENOMIC_PAIRS N_DIMER_PAIRS  (writes sim_R1.fastq.gz, sim_R2.fastq.gz)

Genomic pairs come from a random 400 kb genome (10% with short inserts, so the adapter is
read through). Dimer pairs mimic the Nepal V. cholerae data: reads start at GATCGG...,
R1 runs off into Q2 poly-A after the adapter, R2 mates are low quality throughout.
"""
import gzip
import random
import sys

random.seed(7)
L = 101
A1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" + "ATCACG" + "ATCTCGTATGCCGTCTTCTGCTTG"
A2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


G = "".join(random.choice("ACGT") if random.random() > 0.05 else random.choice("AT")
            for _ in range(400000))


def qual(n, start, end, sd=3):
    return [max(2, min(40, int(random.gauss(start + (end - start) * i / n, sd)))) for i in range(n)]


def mutate(s, q):
    return "".join(random.choice("ACGT") if random.random() < 10 ** (-qq / 10) else b
                   for b, qq in zip(s, q, strict=True))


def rec(name, s, q):
    return f"@{name}\n{s}\n+\n{''.join(chr(33 + x) for x in q)}\n"


n_gen, n_dim = int(sys.argv[1]), int(sys.argv[2])
order = ["g"] * n_gen + ["d"] * n_dim
random.shuffle(order)
with gzip.open("sim_R1.fastq.gz", "wt") as f1, gzip.open("sim_R2.fastq.gz", "wt") as f2:
    for i, t in enumerate(order):
        if t == "g":
            ins = random.randint(60, 140) if random.random() < 0.1 else random.randint(150, 600)
            p = random.randint(0, len(G) - ins)
            frag = G[p:p + ins]
            s1 = (frag + A1 + "A" * L)[:L]
            s2 = (rc(frag) + A2 + "A" * L)[:L]
            q1 = qual(L, 38, 30)
            q2 = qual(L, 36, 26, 4)
        else:
            s1 = (A1[1:] + "A" * L)[:L]
            s2 = (A2[1:] + "A" * L)[:L]
            q1 = qual(62, 38, 34) + [2] * (L - 62)
            q2 = qual(L, 18, 10, 5)
        name = f"r{i}_{'gen' if t == 'g' else 'dimer'}"
        f1.write(rec(f"{name}/1", mutate(s1, q1), q1))
        f2.write(rec(f"{name}/2", mutate(s2, q2), q2))
