# -*- coding: utf-8 -*-
"""从已有 CSV 生成 summary.json（不重跑仿真）。"""
import os
import json
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")


def summarize(df, group):
    g = df.groupby(group)
    return pd.DataFrame({
        "n": g.size(),
        "clear_rate": g["cleared"].mean(),
        "relapse_rate": g["relapsed"].mean(),
        "median_final": g["final_cancer"].median(),
        "median_peak": g["peak_cancer"].median(),
        "median_clear_day": g["cleared_day"].median(),
    }).reset_index()


d1 = pd.read_csv(os.path.join(RES, "E1_phase_constant.csv"))
d2 = pd.read_csv(os.path.join(RES, "E2_phase_proportional.csv"))
d3 = pd.read_csv(os.path.join(RES, "E3_ablation.csv"))
d4 = pd.read_csv(os.path.join(RES, "E4_safety.csv"))
d5 = pd.read_csv(os.path.join(RES, "E5_montecarlo.csv"))

safety = summarize(d4, ["kappa_off"])
safety["mean_off_target"] = d4.groupby("kappa_off")["off_target"].mean().values

# 50% 清除阈值：在相邻 α 之间线性插值
def threshold(mat):
    out = {}
    for r in mat.columns:
        col = mat[r]
        xs, ys = list(col.index), list(col.values * 100)
        if max(ys) < 50:
            out[float(r)] = None
            continue
        if ys[0] >= 50:
            out[float(r)] = float(xs[0])
            continue
        for i in range(1, len(xs)):
            if ys[i] >= 50:
                x0, x1, y0, y1 = xs[i - 1], xs[i], ys[i - 1], ys[i]
                out[float(r)] = float(x0 + (50 - y0) * (x1 - x0) / (y1 - y0))
                break
    return out


m1 = d1.pivot_table(index="alpha", columns="r", values="cleared", aggfunc="mean").sort_index()
m2 = d2.pivot_table(index="alpha", columns="r", values="cleared", aggfunc="mean").sort_index()

summary = {
    "n_runs": int(len(d1) + len(d2) + len(d3) + len(d4) + len(d5)),
    "phase_constant": summarize(d1, ["r", "alpha"]).to_dict("records"),
    "phase_proportional": summarize(d2, ["r", "alpha"]).to_dict("records"),
    "ablation": summarize(d3, ["_arm", "_level"]).to_dict("records"),
    "safety": safety.to_dict("records"),
    "monte_carlo": {
        "n": int(len(d5)),
        "clear_rate": float(d5["cleared"].mean()),
        "relapse_rate": float(d5["relapsed"].mean()),
        "median_clear_day": float(d5["cleared_day"].median()),
        "mean_clear_day": float(d5["cleared_day"].mean()),
        "std_clear_day": float(d5["cleared_day"].std()),
        "p90_clear_day": float(d5["cleared_day"].quantile(0.9)),
    },
    "alpha50_constant": threshold(m1),
    "alpha50_proportional": threshold(m2),
}
with open(os.path.join(RES, "summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)
print(json.dumps({k: summary[k] for k in ("n_runs", "monte_carlo",
                                          "alpha50_constant", "alpha50_proportional")},
                 ensure_ascii=False, indent=2))
