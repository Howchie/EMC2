#!/usr/bin/env python3
"""Summarise sampled tree memory and rebuild evaluation logs.

Sampling can miss short-lived children. Tree totals are approximate: samples
are grouped into 0.1 s intervals and a process is counted at most once per bin.
"""
import collections
import csv
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
for arm in ('base', 'new'):
    path = root / f'EXG_{arm}_memory.csv'
    bins = collections.defaultdict(dict)
    with path.open() as stream:
        for row in csv.DictReader(stream):
            key = round(float(row['elapsed_s']) * 10)
            bins[key][row['pid']] = (int(row['pss_kb']), int(row['private_kb']))
    print(arm, 'sampled tree peak PSS KiB', max(sum(x[0] for x in b.values()) for b in bins.values()),
          'private KiB', max(sum(x[1] for x in b.values()) for b in bins.values()),
          'max sampled processes', max(len(b) for b in bins.values()))
    records = collections.defaultdict(list)
    for p in (root / f'EXG_{arm}_counts').glob('*.log'):
        for line in p.read_text().splitlines():
            items = line.split()
            if len(items) == 2 and items[1] != 'NA':
                records[items[0]].append(float(items[1]))
    for name, values in records.items():
        print(' ', name, 'calls', len(values), 'sum', sum(values))
