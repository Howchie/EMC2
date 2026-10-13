"""Record per-process RSS during a benchmark command (Linux)."""
import csv
from pathlib import Path
import subprocess
import sys
import time

output, *command = sys.argv[1:]
process = subprocess.Popen(command)
with open(output, "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["elapsed", "pid", "ppid", "rss_kb"])
    start = time.monotonic()
    while process.poll() is None:
        records = {}
        for directory in Path("/proc").iterdir():
            if not directory.name.isdigit():
                continue
            try:
                status = (directory / "status").read_text()
                fields = dict(line.split(":", 1) for line in status.splitlines() if ":" in line)
                records[int(directory.name)] = (int(fields["PPid"]),
                    int(fields.get("VmRSS", "0 kB").split()[0]))
            except (OSError, ValueError, KeyError):
                continue
        descendants = {process.pid}
        while True:
            found = {pid for pid, (parent, _) in records.items() if parent in descendants}
            if found <= descendants:
                break
            descendants |= found
        # Sampling templates can be reparented to the caller's subreaper.
        for pid in sorted(descendants):
            if pid in records:
                parent, rss = records[pid]
                writer.writerow([time.monotonic() - start, pid, parent, rss])
        handle.flush()
        time.sleep(0.2)
sys.exit(process.returncode)
