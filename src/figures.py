# -*- coding: utf-8 -*-
"""根据 results/*.csv 生成论文插图。

运行：python src/figures.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tcell_model import Sim, SimConfig  # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
IMG = os.path.join(BASE, "img")
os.makedirs(IMG, exist_ok=True)

FONT = r"C:\Windows\Fonts\msyh.ttc"
if os.path.exists(FONT):
    font_manager.fontManager.addfont(FONT)
    plt.rcParams["font.family"] = "Microsoft YaHei"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["xtick.major.width"] = 0.8
plt.rcParams["ytick.major.width"] = 0.8

GREEN = "#1B7F4B"
LGREEN = "#6BBF8E"
GREY = "#8A8A8A"
DARK = "#222222"


def clear_matrix(df):
    p = df.pivot_table(index="alpha", columns="r", values="cleared", aggfunc="mean")
    return p.sort_index().sort_index(axis=1)


def draw_heat(ax, mat, title):
    r_vals = list(mat.columns)
    a_vals = list(mat.index)
    data = mat.values * 100
    cmap = LinearSegmentedColormap.from_list("clear", ["#F4F4F5", "#CDE8D8", "#6BBF8E", "#1B7F4B"])
    im = ax.imshow(data, origin="lower", aspect="auto", cmap=cmap, vmin=0, vmax=100)
    ax.set_xticks(range(len(r_vals)))
    ax.set_xticklabels([f"{v:.2f}" for v in r_vals], fontsize=8)
    ax.set_yticks(range(len(a_vals)))
    ax.set_yticklabels([f"{v:.2f}" for v in a_vals], fontsize=8)
    ax.set_xlabel("肿瘤生长率 r", fontsize=9)
    ax.set_ylabel("诱导效率 α", fontsize=9)
    for i in range(len(a_vals)):
        for j in range(len(r_vals)):
            v = data[i, j]
            ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=7.5,
                    color="white" if v > 55 else DARK)
    ax.set_title(title, fontsize=9.5, pad=6)
    return im


def fig1():
    d1 = pd.read_csv(os.path.join(RES, "E1_phase_constant.csv"))
    d2 = pd.read_csv(os.path.join(RES, "E2_phase_proportional.csv"))
    m1, m2 = clear_matrix(d1), clear_matrix(d2)

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.5), sharey=True)
    draw_heat(axes[0], m1, "A  恒定供给（每日 0.05 个）")
    im = draw_heat(axes[1], m2, "B  比例供给（随瘤负荷增长）")
    cax = fig.add_axes([0.92, 0.18, 0.015, 0.64])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("清除率（%）", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.subplots_adjust(left=0.09, right=0.90, top=0.90, bottom=0.16, wspace=0.12)
    out = os.path.join(IMG, "Fig1_phase_boundary.png")
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("已保存", out)
    return m1, m2


def fig2():
    d = pd.read_csv(os.path.join(RES, "E3_ablation.csv"))
    arms = [("alpha", "诱导效率 α", "A", None),
            ("kappa", "杀伤效率 κ", "B", None),
            ("inflow", "恒定供给量（个/日）", "C", ["0.02", "0.05", "0.2", "0.5", "1", "2"]),
            ("lam_scale", "比例供给系数（$\\times10^{-3}$）", "D", ["0.5", "1", "2", "4", "8"])]
    fig, axes = plt.subplots(1, 4, figsize=(8.6, 2.9), sharey=True)
    for ax, (arm, label, letter, ticklabels) in zip(axes, arms):
        sub = d[d["_arm"] == arm]
        g = sub.groupby("_level")["cleared"].mean() * 100
        ax.plot(g.index, g.values, "-o", color=GREEN, ms=3.6, lw=1.3)
        for x, y in g.items():
            ax.annotate(f"{y:.0f}", (x, y), textcoords="offset points",
                        xytext=(0, 6), ha="center", fontsize=6.2, color=DARK)
        ax.set_xlabel(label, fontsize=8)
        ax.set_title(letter, fontsize=9, loc="left")
        ax.set_ylim(-6, 112)
        ax.grid(color="#E6E6E6", lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        vals = sorted(sub["_level"].unique())
        if ticklabels is not None:
            ax.set_xscale("log")
        ax.set_xticks(vals)
        if ticklabels is not None:
            ax.set_xticklabels(ticklabels, fontsize=6.8)
        else:
            ax.set_xticklabels([f"{v:g}" for v in vals], fontsize=7)
    axes[0].set_ylabel("清除率（%）", fontsize=9)
    fig.text(0.5, -0.02, "工作点：r = 0.20，α = 0.10，κ = 0.80，恒定供给 0.05 个/日；每个水平 30 次重复",
             ha="center", fontsize=6.8, color="#555555")
    fig.subplots_adjust(left=0.08, right=0.99, top=0.88, bottom=0.24, wspace=0.16)
    out = os.path.join(IMG, "Fig2_ablation.png")
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("已保存", out)


def fig3():
    d = pd.read_csv(os.path.join(RES, "E4_safety.csv"))
    d["_level"] = d["kappa_off"]
    g = d.groupby("_level").agg(clear=("cleared", "mean"),
                                off=("off_target", "mean"),
                                off_med=("off_target", "median")).reset_index()
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    ax.plot(g["off"], g["clear"] * 100, "-o", color=GREEN, lw=1.4, ms=4.2)
    off_map = {0.10: ((0, -17), "center"),
               0.20: ((-3, 11), "right"),
               0.30: ((3, 11), "left")}
    for _, row in g.iterrows():
        (dx, dy), ha = off_map.get(round(row["_level"], 2), ((6, 5), "left"))
        ax.annotate(f"κ_off={row['_level']:.2f}", (row["off"], row["clear"] * 100),
                    textcoords="offset points", xytext=(dx, dy),
                    ha=ha, fontsize=6.5, color="#444444")
    ax.set_xlabel("每次仿真中正常细胞被误杀的平均数（个）", fontsize=8.5)
    ax.set_ylabel("肿瘤清除率（%）", fontsize=8.5)
    ax.set_ylim(55, 100)
    ax.set_xlim(-8, 150)
    ax.grid(color="#E6E6E6", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.subplots_adjust(left=0.16, right=0.98, top=0.95, bottom=0.16)
    out = os.path.join(IMG, "Fig3_safety_tradeoff.png")
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("已保存", out)


def fig4():
    d = pd.read_csv(os.path.join(RES, "E5_montecarlo.csv"))
    days = d["cleared_day"].dropna().astype(float)

    sat = SimConfig(r=0.25, alpha=0.50, kappa=0.95, lam_scale=0.02,
                    supply_mode="proportional", chemo_days=(30, 60, 90),
                    chemo_kill=0.60, seed=7)
    ctrl = SimConfig(r=0.25, alpha=0.10, kappa=0.80, inflow=0.05,
                     supply_mode="constant", seed=7)
    sim1 = Sim(sat)
    sim1.run()
    sim2 = Sim(ctrl)
    sim2.run()
    h1 = sim1.history
    h2 = sim2.history

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1))
    ax = axes[0]
    ax.plot(range(1, len(h2) + 1), h2, color=GREY, lw=1.3, label="单药（α=0.10, κ=0.80）")
    ax.plot(range(1, len(h1) + 1), h1, color=GREEN, lw=1.5, label="饱和方案（α=0.50, κ=0.95）")
    for dday in sat.chemo_days:
        ax.axvline(dday, color="#C9C9C9", ls=":", lw=0.8)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(0, 2500)
    ax.set_yticks([0, 1, 10, 100, 1000])
    ax.set_yticklabels(["0", "1", "10", r"$10^2$", r"$10^3$"])
    ax.set_xlabel("时间（天）", fontsize=8.5)
    ax.set_ylabel("肿瘤细胞数", fontsize=8.5)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    ax.set_title("A  典型轨迹", fontsize=9, loc="left")
    ax.grid(color="#EFEFEF", lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    ax = axes[1]
    ax.hist(days, bins=16, color=LGREEN, edgecolor="white", linewidth=0.8)
    ax.axvline(days.median(), color=GREEN, lw=1.4)
    ax.annotate(f"中位 {days.median():.1f} 天", (days.median(), ax.get_ylim()[1] * 0.9),
                textcoords="offset points", xytext=(6, 0), fontsize=7, color=GREEN)
    ax.set_xlabel("清除时间（天）", fontsize=8.5)
    ax.set_ylabel("次数", fontsize=8.5)
    ax.set_title("B  清除时间分布（n = %d）" % len(days), fontsize=9, loc="left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.18, wspace=0.30)
    out = os.path.join(IMG, "Fig4_monte_carlo.png")
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("已保存", out)


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
