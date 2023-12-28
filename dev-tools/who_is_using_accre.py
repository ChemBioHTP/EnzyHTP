"""ad hoc tool for checking ACCRE queue.

Usage:
1. squeue|grep productio > test.log
2. python who_is_using_accre.py"""
from collections import defaultdict 

top_running_user = defaultdict(int)
top_pending_user = defaultdict(int)
not_counted = 0

with open("test.log") as f:
    fc = f.readlines()
    for line in fc:
        lp = line.strip().split()
        if len(lp) == 8:
            job_id, partition, job_name, user, state, time, num_node, nodelist = lp
            if state == "R":
                top_running_user[user] += 1
            else:
                top_pending_user[user] += 1
        else:
            not_counted += 1
top_running_user = dict(sorted(top_running_user.items(), key=lambda item: item[1], reverse=1))
top_pending_user = dict(sorted(top_pending_user.items(), key=lambda item: item[1], reverse=1))

print("top_running_user\n")
print(*top_running_user.items(), sep="\n")
print("top_pending_user\n")
print(*top_pending_user.items(), sep="\n")
