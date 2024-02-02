#!/usr/bin/env python3
import os
import sys

CS = set()
with open("test.prot", "r") as prot:
    for line in prot:
        line.strip()
        if line[0] == '#':
            continue
        l = line.split()
        CS.add(l[1])

with open("test.peaks", "r") as pk:
    for line in pk:
        line.strip()
        if line[0] == '#':
            continue
        l = line.split()
        if len(l) >= 4:
            c1 = l[1]
            c2 = l[2]
            c3 = l[3]
            if c1 in CS and c2 in CS and c3 in CS:
                print(line)

