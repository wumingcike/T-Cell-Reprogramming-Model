# -*- coding: utf-8 -*-
"""v2 实验套件：相图、单变量消融、安全性权衡、蒙特卡洛。

运行：python src/experiments.py
输出：results/*.csv 与 results/summary.json
"""

import os
import sys
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcell_model import run_one  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "results")
os.makedirs(OUT, exist_ok=True)

N_SEED = 30
BASE_R, BASE_A, BASE_K, BASE_INFLOW = 0.20, 0.10, 0.80, 0.05


def _worker(kw):
    kw = {k: v for k, v in kw.items() if not k.startswith("_")}
    return run_one(**kw)


def execute(tasks, workers=8):
    with ProcessPoolExecutor(max_workers=workers) as ex:
        out = list(ex.map(_worker, tasks, chunksize=4))
    return pd.DataFrame(out)


def phase_map(mode, seeds=N_SEED):
    tasks = []
    for r in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]:
        for a in [0.02, 0.05, 0.10, 0.20, 0.30]:
            for s in range(1, seeds + 1):
                kw = dict(r=r, alpha=a, kappa=BASE_K, seed=s, supply_mode=mode)
                if mode == "constant":
                    kw["inflow"] = BASE_INFLOW
                else:
                    kw["lam_scale"] = 0.002
                tasks.append(kw)
    return execute(tasks), tasks


def ablation(seeds=N_SEED):
    tasks = []
    # 只改诱导效率
    for a in [0.02, 0.05, 0.10, 0.20, 0.30]:
        for s in range(1, seeds + 1):
            tasks.append(dict(r=BASE_R, alpha=a, kappa=BASE_K, inflow=BASE_INFLOW,
                              seed=s, supply_mode="constant", _arm="alpha", _level=a))
    # 只改杀伤效率
    for k in [0.20, 0.40, 0.60, 0.80, 0.95]:
        for s in range(1, seeds + 1):
            tasks.append(dict(r=BASE_R, alpha=BASE_A, kappa=k, inflow=BASE_INFLOW,
                              seed=s, supply_mode="constant", _arm="kappa", _level=k))
    # 只改恒定供给
    for inf in [0.02, 0.05, 0.20, 0.50, 1.00, 2.00]:
        for s in range(1, seeds + 1):
            tasks.append(dict(r=BASE_R, alpha=BASE_A, kappa=BASE_K, inflow=inf,
                              seed=s, supply_mode="constant", _arm="inflow", _level=inf))
    # 只改比例供给系数
    for ls in [0.0005, 0.0010, 0.0020, 0.0040, 0.0080]:
        for s in range(1, seeds + 1):
            tasks.append(dict(r=BASE_R, alpha=BASE_A, kappa=BASE_K, lam_scale=ls,
                              seed=s, supply_mode="proportional", _arm="lam_scale", _level=ls))
    return execute(tasks), tasks


def safety(seeds=20):
    tasks = []
    for ko in [0.0, 0.02, 0.05, 0.10, 0.20, 0.30]:
        for s in range(1, seeds + 1):
            tasks.append(dict(r=BASE_R, alpha=0.20, kappa=BASE_K, inflow=BASE_INFLOW,
                              seed=s, supply_mode="constant", normal_frac=0.35,
                              kappa_off=ko, _level=ko))
    return execute(tasks), tasks


def monte_carlo(seeds=300):
    tasks = []
    for s in range(1, seeds + 1):
        tasks.append(dict(r=0.25, alpha=0.50, kappa=0.95, lam_scale=0.02,
                          supply_mode="proportional", chemo_days=(30, 60, 90),
                          chemo_kill=0.60, seed=s, _level=1.0))
    return execute(tasks), tasks


def summarize(df, group):
    g = df.groupby(group)
    out = pd.DataFrame({
        "n": g.size(),
        "clear_rate": g["cleared"].mean(),
        "relapse_rate": g["relapsed"].mean(),
        "median_final": g["final_cancer"].median(),
        "median_peak": g["peak_cancer"].median(),
        "median_clear_day": g["cleared_day"].median(),
    }).reset_index()
    return out


def main():
    t0 = time.time()
    workers = min(12, (os.cpu_count() or 4))
    print(f"使用 {workers} 个进程", flush=True)

    print("[1/5] 相图：恒定供给", flush=True)
    df1, _ = phase_map("constant")
    df1.to_csv(os.path.join(OUT, "E1_phase_constant.csv"), index=False)

    print("[2/5] 相图：比例供给", flush=True)
    df2, _ = phase_map("proportional")
    df2.to_csv(os.path.join(OUT, "E2_phase_proportional.csv"), index=False)

    print("[3/5] 单变量消融", flush=True)
    df3, tasks3 = ablation()
    meta = pd.DataFrame([{k: v for k, v in t.items() if k.startswith("_")} for t in tasks3])
    df3 = pd.concat([df3.reset_index(drop=True), meta], axis=1)
    df3.to_csv(os.path.join(OUT, "E3_ablation.csv"), index=False)

    print("[4/5] 安全性权衡", flush=True)
    df4, _ = safety()
    df4.to_csv(os.path.join(OUT, "E4_safety.csv"), index=False)

    print("[5/5] 蒙特卡洛", flush=True)
    df5, _ = monte_carlo()
    df5.to_csv(os.path.join(OUT, "E5_montecarlo.csv"), index=False)

    summary = {
        "n_runs": int(len(df1) + len(df2) + len(df3) + len(df4) + len(df5)),
        "phase_constant": summarize(df1, ["r", "alpha"]).to_dict("records"),
        "phase_proportional": summarize(df2, ["r", "alpha"]).to_dict("records"),
        "ablation": summarize(df3, ["_arm", "_level"]).to_dict("records"),
        "safety": summarize(df4, ["_level"]).to_dict("records"),
        "monte_carlo": {
            "n": int(len(df5)),
            "clear_rate": float(df5["cleared"].mean()),
            "relapse_rate": float(df5["relapsed"].mean()),
            "median_clear_day": float(df5["cleared_day"].median()),
            "mean_clear_day": float(df5["cleared_day"].mean()),
            "std_clear_day": float(df5["cleared_day"].std()),
            "p90_clear_day": float(df5["cleared_day"].quantile(0.9)),
        },
        "elapsed_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(json.dumps(summary["monte_carlo"], ensure_ascii=False, indent=2))
    print(f"完成，用时 {summary['elapsed_sec']} s")


if __name__ == "__main__":
    main()
