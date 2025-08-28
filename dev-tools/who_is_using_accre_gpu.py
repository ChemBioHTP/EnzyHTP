#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess, re
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

UNKNOWN = "unknown_model"

def run(cmd):
    return subprocess.check_output(cmd, text=True)

# 1) Get (JobID, User, State) from squeue
rows = []
for line in run(["squeue","-h","-t","R,PD","-o","%i %u %T"]).splitlines():
    parts = line.strip().split(maxsplit=2)
    if len(parts) < 3:
        continue
    jid, user, state = parts
    rows.append((jid, user, state))

# 2) Parse GPU request model and count from scontrol single-line output
pat_reqgres = re.compile(r"\bReqGRES=([^ ]+)")
pat_gres    = re.compile(r"\bGres=([^ ]+)")
pat_reqtres = re.compile(r"\bReqTRES=([^ ]+)")
pat_tpn     = re.compile(r"\bTresPerNode=([^ ]+)")

# Parse formats like "gpu:nvidia_titan_x:1,gpu:1"
def parse_from_reqgres(req: str):
    model_counts = defaultdict(int)
    for tok in req.split(","):
        tok = tok.strip()
        if not tok.startswith("gpu"):
            continue
        # Allow gpu:nvidia_rtx_a4000:2 / gpu:2
        parts = tok.split(":")
        # parts[0] == "gpu"
        if len(parts) == 3 and parts[2].isdigit():
            model = parts[1] or UNKNOWN
            model_counts[model] += int(parts[2])
        elif len(parts) == 2 and parts[1].isdigit():
            model_counts[UNKNOWN] += int(parts[1])
    return dict(model_counts)

# Parse ReqTRES= "...,gres/gpu=1,gres/gpu:nvidia_titan_x=1,..."
def parse_from_reqtres(tres: str):
    model_counts = defaultdict(int)
    generic = 0
    for tok in tres.split(","):
        tok = tok.strip()
        if tok.startswith("gres/gpu"):
            # Format: gres/gpu[:MODEL]=N
            m = re.match(r"gres/gpu(?::([^=,]+))?=(\d+)$", tok)
            if m:
                model, cnt = m.group(1), int(m.group(2))
                if model:
                    model_counts[model] += cnt
                else:
                    generic += cnt
    if model_counts:
        return dict(model_counts)  # Already have specific models, ignore generic entries to avoid duplication
    if generic:
        return {UNKNOWN: generic}
    return {}

# Parse TresPerNode="gres/gpu:nvidia_titan_x:1,..." (note: colon-separated count here)
def parse_from_tpn(tpn: str):
    model_counts = defaultdict(int)
    generic = 0
    for tok in tpn.split(","):
        tok = tok.strip()
        if tok.startswith("gres/gpu"):
            # Format: gres/gpu[:MODEL]:N
            m = re.match(r"gres/gpu(?::([^:,=]+))?:(\d+)$", tok)
            if m:
                model, cnt = m.group(1), int(m.group(2))
                if model:
                    model_counts[model] += cnt
                else:
                    generic += cnt
    if model_counts:
        return dict(model_counts)
    if generic:
        return {UNKNOWN: generic}
    return {}

def parse_models_counts(jobline: str):
    # 1) ReqGRES / Gres
    m = pat_reqgres.search(jobline) or pat_gres.search(jobline)
    if m:
        res = parse_from_reqgres(m.group(1))
        if res:
            return res
    # 2) ReqTRES
    m = pat_reqtres.search(jobline)
    if m:
        res = parse_from_reqtres(m.group(1))
        if res:
            return res
    # 3) TresPerNode
    m = pat_tpn.search(jobline)
    if m:
        res = parse_from_tpn(m.group(1))
        if res:
            return res
    return {}

def fetch(job_tuple):
    jid, user, state = job_tuple
    try:
        # -o ensures single-line key=value format for easy parsing
        out = run(["scontrol","show","job","-o", jid])
    except subprocess.CalledProcessError:
        return []
    models = parse_models_counts(out)
    if not models:
        return []
    items = []
    # Standardize State
    s = state
    if s in ("RUNNING", "R"):
        s = "R"
    elif s in ("PENDING", "PD"):
        s = "PD"
    else:
        return []
    for model, cnt in models.items():
        items.append((user, s, model, int(cnt)))
    return items

# 3) Concurrent fetch and aggregation
all_items = []
with ThreadPoolExecutor(max_workers=24) as ex:
    futs = [ex.submit(fetch, r) for r in rows]
    for f in as_completed(futs):
        all_items.extend(f.result())

by_model = defaultdict(lambda: {"R":0,"PD":0})
by_user_model = defaultdict(lambda: {"R":0,"PD":0})

for user, state, model, cnt in all_items:
    by_model[model][state] += cnt
    by_user_model[(user, model)][state] += cnt

# 4) Print results
print(f"{'GPU Model':30} {'Running(R)':>10} {'Pending(PD)':>12}")
for model in sorted(by_model):
    u = by_model[model]["R"]
    q = by_model[model]["PD"]
    print(f"{model:30} {u:10d} {q:12d}")

print("\n=== By User - Model ===")
print(f"{'User':20} {'GPU Model':30} {'R':>6} {'PD':>6}")
for (user, model) in sorted(by_user_model):
    u = by_user_model[(user,model)]["R"]
    q = by_user_model[(user,model)]["PD"]
    print(f"{user:20} {model:30} {u:6d} {q:6d}")
