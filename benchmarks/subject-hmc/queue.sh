#!/usr/bin/env bash
# Run benchmark jobs one at a time, pinned to 27 cores (3 chains x 9): one
# fit scales better with cores than several smaller fits side by side, and
# the allocator prices moves by wall-clock, so concurrent load changes plans,
# not just CPU.
#   queue.sh <jobs.txt>
# Each line of jobs.txt is a shell command; CORES is substituted with the slot's
# core list. Lines starting with # are skipped.
jobs=$1
slots=("0-26")
pids=("")
run_line() {
  local line=$1 slot=$2
  local cmd=${line//CORES/${slots[$slot]}}
  bash -c "$cmd" &
  pids[$slot]=$!
}
while IFS= read -r line || [ -n "$line" ]; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  while true; do
    for s in 0; do
      if [ -z "${pids[$s]}" ] || ! kill -0 "${pids[$s]}" 2>/dev/null; then
        run_line "$line" $s
        continue 3
      fi
    done
    sleep 30
  done
done < "$jobs"
wait
