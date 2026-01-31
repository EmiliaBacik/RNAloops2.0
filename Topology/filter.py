#!/usr/bin/env python3
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-vtopol", required=True, help="wartość topol w formacie X-X-X")
args = parser.parse_args()

topol = args.vtopol
ta = topol.split("-")
ta_len = len(ta)

for line_number, line in enumerate(sys.stdin, start=1):
    if line_number <= 2:
        continue

    line = line.rstrip("\n")
    fields = line.split(";")

    if len(fields) < 4:
        continue

    a = fields[3].split("-")

    if len(a) != ta_len:
        continue

    b_parts = []
    for segment in a:
        count = segment.count(".")
        b_parts.append(str(count))

    b = "-".join(b_parts)

    if b == topol:
        print(f"{line};{b}")
        continue

    for j in range(1, ta_len):
        rotated = []
        for i in range(ta_len):
            idx = (j + i) % ta_len
            rotated.append(ta[idx])
        d = "-".join(rotated)

        if d == b:
            print(f"{line};{b}")
            break
