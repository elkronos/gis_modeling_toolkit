#!/bin/bash
# Runs one command with a time limit. If it is still running at the limit,
# records each process's state and native stack, then kills it.
# Usage: probe.sh <name> <limit in seconds> <directory> <command...>
name=$1; limit=$2; dir=$3; shift 3
logs="${FEDORA_LOGS:?}"; mkdir -p "$logs"
log="$logs/$name.log"; stacks="$logs/$name.stacks.txt"

tree() {  # a process and all its descendants
  local p
  echo "$1"
  for p in $(pgrep -P "$1" 2>/dev/null); do tree "$p"; done
}

start=$SECONDS
( cd "$dir" && exec "$@" ) > "$log" 2>&1 &
pid=$!
while kill -0 "$pid" 2>/dev/null && (( SECONDS - start < limit )); do sleep 2; done

if kill -0 "$pid" 2>/dev/null; then
  result=HANG; rc=NA
  procs=$(tree "$pid")
  {
    echo "== process states, 10 seconds apart"
    for round in 1 2; do
      for p in $procs; do
        ps -L -o pid,tid,pcpu,cputime,stat,wchan:24,comm -p "$p" 2>/dev/null
      done
      [ "$round" = 1 ] && sleep 10
    done
    for p in $procs; do
      echo; echo "== native stacks of process $p"
      timeout 120 gdb -p "$p" -batch -nx -ex 'set pagination off' \
        -ex 'info threads' -ex 'thread apply all bt 45' 2>&1
    done
  } > "$stacks"
  for p in $procs; do kill -KILL "$p" 2>/dev/null; done
  wait "$pid" 2>/dev/null
else
  wait "$pid"; rc=$?
  if [ "$rc" = 0 ]; then result=ok; else result=FAILED; fi
fi

secs=$(( SECONDS - start ))
printf '%s\t%s\t%s\t%s\n' "$name" "$result" "$rc" "$secs" >> "$logs/summary.tsv"
echo "=== $name: $result (exit $rc) after $secs s"
tail -n 12 "$log"
if [ -f "$stacks" ]; then
  echo "--- process states and stacks"
  head -n 220 "$stacks"
fi
exit 0
