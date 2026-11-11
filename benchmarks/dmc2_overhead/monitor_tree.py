#!/usr/bin/env python3
"""Sample a command's Linux process tree; fast-lived processes may be missed."""
import csv
import os
import subprocess
import sys
import time

if len(sys.argv) < 4 or sys.argv[2] != '--':
    raise SystemExit('usage: monitor_tree.py output.csv -- command [args...]')
child = subprocess.Popen(sys.argv[3:])
with open(sys.argv[1], 'w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['elapsed_s', 'pid', 'ppid', 'rss_kb', 'pss_kb', 'private_kb', 'command'])
    start = time.monotonic()
    while child.poll() is None:
        pending = [child.pid]
        seen = set()
        while pending:
            pid = pending.pop()
            if pid in seen:
                continue
            seen.add(pid)
            try:
                with open(f'/proc/{pid}/stat') as f:
                    stat = f.read()
                tail = stat[stat.rfind(')') + 2:].split()
                ppid = int(tail[1]); rss = int(tail[21]) * os.sysconf('SC_PAGE_SIZE') // 1024
                with open(f'/proc/{pid}/smaps_rollup') as f:
                    fields = dict(line.split(':', 1) for line in f if ':' in line)
                value = lambda k: int(fields.get(k, '0 kB').split()[0])
                pss = value('Pss'); private = value('Private_Clean') + value('Private_Dirty')
                with open(f'/proc/{pid}/cmdline', 'rb') as f:
                    command = f.read().replace(b'\0', b' ').decode(errors='replace').strip()
                writer.writerow([round(time.monotonic() - start, 3), pid, ppid, rss, pss, private, command])
                with open(f'/proc/{pid}/task/{pid}/children') as f:
                    pending.extend(map(int, f.read().split()))
            except (OSError, ValueError, IndexError):
                continue
        stream.flush()
        time.sleep(.1)
sys.exit(child.returncode)
