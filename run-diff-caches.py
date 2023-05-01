#!/bin/python3
import os
import subprocess
import csv
config_path = "/home/am6429/research/Reverse_Time_Migration/workloads/bp_model/veloc_config.cfg"
res = []
times = []
for g in [x/10.0 for x in range(5, 61, 5)]:
    t = []
    for c in [0, 1]:
        print(f"Running for GPU cache: {g}, compress: {c}")
        os.system('killall veloc-backend')
        os.system('rm -rf /data/scratch/*')
        f = "false" if c == 1 else "true"
        os.system(f"sed -i 's/gpu_cache_size=.*/gpu_cache_size={g}/g' {config_path}")
        os.system(f"sed -i 's/force_no_compress=.*/force_no_compress={f}/g' {config_path}")
        process = subprocess.Popen(['./bin/Engine'],
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        temp = stdout.decode('ascii').split("<<<<<=====")[1]
        temp = temp.split("=====>>>>>")[0]
        temp = [float(num) for num in temp.strip().split(',')]
        tot_ckpt_res_time = temp[2]+temp[3]
        tot_comp_decomp_time = temp[4]+temp[5]
        ans = {
            "gpu_cache": g,
            "compress": c,
            "ckpt_res": tot_ckpt_res_time,
            "comp_decomp": tot_comp_decomp_time,
            "comp_to_ckpt": tot_comp_decomp_time/tot_ckpt_res_time,
            "log": temp,
        }
        print(ans)
        # times[g][c] = tot_ckpt_res_time
        t.append(tot_ckpt_res_time)
        res.append(ans)
    ans = {
        "gpu_cache": g,
        "no_compress": t[0],
        "compress": t[1],
    }
    times.append(ans)


with open('diff-cache-log.csv', 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames = res[0].keys())
    writer.writeheader()
    writer.writerows(res)

print(times)
with open('diff-cache-log.csv', 'a') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames = times[0].keys())
    writer.writeheader()
    writer.writerows(times)