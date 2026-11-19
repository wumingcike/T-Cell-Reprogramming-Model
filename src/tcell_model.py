# -*- coding: utf-8 -*-
"""自体诱导型 T 细胞 - 肿瘤相互作用模型（v2）

相对 v1 的主要修改
------------------
1. 修复更新不一致：v1 的癌细胞读取旧网格、T 细胞写入新网格，同一格可能被重复占用。
   本版统一采用「随机顺序异步更新」（random sequential update），所有个体在当步只能看到已发生的更新。
2. 补上论文声称但 v1 未实现的机制：诱导型 T 细胞可在抗原刺激下原位增殖。
3. 补上正常细胞群体，使「只杀癌细胞」这一安全性主张可以被检验；同时引入显式脱靶杀伤概率 κ_off，
   并明确这是对「模型不具备分子选择性」这一事实的补偿参数，而非真实选择性机制。
4. 免疫供给提供两种标度律：
   - constant     ：每日流入恒定（v1 的设定，最不利情形）
   - proportional ：每日流入正比于肿瘤负荷（更接近「浸润随瘤负荷增长」的情形）
5. 化疗改为按日程批量减瘤，不再是随机单细胞杀伤。

所有参数集中在 SimConfig 中，运行结果可复现（给定 seed）。
"""

from dataclasses import dataclass, asdict
import numpy as np

EMPTY, NORMAL, CANCER, NAIVE, INDUCED = 0, 1, 2, 3, 4

OFFSETS = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))


@dataclass
class SimConfig:
    L: int = 60                     # 网格边长，L*L 个位点
    r: float = 0.05                 # 癌细胞每步分裂概率
    alpha: float = 0.30             # 初始 T 细胞被诱导的概率
    kappa: float = 0.80             # 诱导型 T 细胞对相邻癌细胞的杀伤概率
    kappa_off: float = 0.0          # 对相邻正常细胞的脱靶杀伤概率
    rho: float = 0.30               # 诱导型 T 细胞在抗原刺激下的原位增殖概率
    max_gen: int = 4                # 克隆扩增代数上限（效应 T 细胞的有限扩增）
    t_life: int = 25                # 效应 T 细胞寿命（步）
    inflow: float = 0.05            # 每日流入的初始 T 细胞数（constant 模式）
    supply_mode: str = "constant"   # constant | proportional
    lam_scale: float = 0.002        # proportional 模式下，每日流入 = lam_scale * 瘤负荷
    lam_cap: float = 20.0           # proportional 模式下的流入上限
    init_tumor: int = 6             # 初始肿瘤边长（init_tumor^2 个细胞）
    init_naive: int = 20            # 初始 T 细胞数
    normal_frac: float = 0.0        # 正常细胞占非肿瘤位点的比例
    chemo_days: tuple = ()          # 化疗日
    chemo_kill: float = 0.60        # 每次化疗清除的瘤负荷比例
    max_days: int = 200
    seed: int = 0


class Sim:
    def __init__(self, cfg: SimConfig):
        self.cfg = cfg
        self.rng = np.random.default_rng(cfg.seed)
        L = cfg.L
        self.grid = np.zeros((L, L), dtype=np.int8)
        self.age = np.zeros((L, L), dtype=np.int16)
        self.gen = np.zeros((L, L), dtype=np.int8)

        c = L // 2
        h = cfg.init_tumor // 2
        lo, hi = c - h, c - h + cfg.init_tumor
        self.grid[lo:hi, lo:hi] = CANCER

        if cfg.normal_frac > 0:
            free = np.argwhere(self.grid == EMPTY)
            n = int(len(free) * cfg.normal_frac)
            if n:
                pick = self.rng.choice(len(free), size=n, replace=False)
                for i in pick:
                    self.grid[free[i, 0], free[i, 1]] = NORMAL

        self._spawn_naive(cfg.init_naive)

        self.peak = int((self.grid == CANCER).sum())
        self.cleared_day = None
        self.off_target = 0
        self.history = []

    # ---------- 辅助 ----------
    def _spawn_naive(self, n, near=True):
        """流入的初始 T 细胞归巢到肿瘤边缘（趋化），而非均匀撒在全网格。

        v1 采用全网格随机放置，会引入「肿瘤越大越容易被随机游走撞上」的伪效应；
        改为边缘归巢后，供给量才与局部免疫压力直接对应。
        """
        if n <= 0:
            return
        rng, grid, L = self.rng, self.grid, self.cfg.L
        cancers = None
        for _ in range(n):
            placed = False
            if near:
                if cancers is None or len(cancers) == 0:
                    cancers = np.argwhere(grid == CANCER)
                if len(cancers):
                    for _ in range(8):
                        cx, cy = cancers[rng.integers(len(cancers))]
                        nx = cx + int(rng.integers(-3, 4))
                        ny = cy + int(rng.integers(-3, 4))
                        if 0 <= nx < L and 0 <= ny < L and grid[nx, ny] == EMPTY:
                            grid[nx, ny] = NAIVE
                            placed = True
                            break
            if not placed:
                empty = np.argwhere(grid == EMPTY)
                if len(empty) == 0:
                    return
                x, y = empty[rng.integers(len(empty))]
                grid[x, y] = NAIVE

    def _adjacent(self, x, y, kind):
        L = self.cfg.L
        for dx, dy in OFFSETS:
            nx, ny = x + dx, y + dy
            if 0 <= nx < L and 0 <= ny < L and self.grid[nx, ny] == kind:
                return nx, ny
        return None

    def _adjacent_list(self, x, y, kind):
        L = self.cfg.L
        out = []
        for dx, dy in OFFSETS:
            nx, ny = x + dx, y + dy
            if 0 <= nx < L and 0 <= ny < L and self.grid[nx, ny] == kind:
                out.append((nx, ny))
        return out

    def _empty_neighbor(self, x, y):
        L = self.cfg.L
        cand = []
        for dx, dy in OFFSETS:
            nx, ny = x + dx, y + dy
            if 0 <= nx < L and 0 <= ny < L and self.grid[nx, ny] == EMPTY:
                cand.append((nx, ny))
        if not cand:
            return None
        return cand[self.rng.integers(len(cand))]

    # ---------- 单步 ----------
    def step(self, day):
        cfg, rng, grid, age, gen = self.cfg, self.rng, self.grid, self.age, self.gen
        L = cfg.L
        K = L * L

        # 化疗：按日程批量减瘤
        if day in cfg.chemo_days:
            cancer = np.argwhere(grid == CANCER)
            if len(cancer):
                n_kill = int(len(cancer) * cfg.chemo_kill)
                if n_kill:
                    pick = rng.choice(len(cancer), size=n_kill, replace=False)
                    for i in pick:
                        x, y = cancer[i]
                        grid[x, y] = EMPTY

        # 免疫细胞流入
        nc = int((grid == CANCER).sum())
        if cfg.supply_mode == "proportional":
            lam = min(cfg.lam_cap, cfg.lam_scale * max(nc, 1))
        else:
            lam = cfg.inflow
        n_in = int(rng.poisson(lam)) if lam > 0 else 0
        if n_in:
            self._spawn_naive(n_in)

        # 随机顺序异步更新
        occupied = np.argwhere(grid > 0)
        order = rng.permutation(len(occupied))
        growth_p = cfg.r * max(0.0, 1.0 - nc / K)

        for oi in order:
            x, y = occupied[oi]
            t = grid[x, y]
            if t == CANCER:
                if growth_p > 0 and rng.random() < growth_p:
                    tgt = self._empty_neighbor(x, y)
                    if tgt is not None:
                        grid[tgt[0], tgt[1]] = CANCER
            elif t == NAIVE:
                cancers = self._adjacent_list(x, y, CANCER)
                if cancers and rng.random() < cfg.alpha:
                    grid[x, y] = INDUCED
                    age[x, y] = 0
                    gen[x, y] = 0
                else:
                    tgt = self._empty_neighbor(x, y)
                    if tgt is not None:
                        grid[x, y] = EMPTY
                        grid[tgt[0], tgt[1]] = NAIVE
            elif t == INDUCED:
                age[x, y] += 1
                if age[x, y] > cfg.t_life:
                    grid[x, y] = EMPTY
                    age[x, y] = 0
                    continue
                killed = False
                cancers = self._adjacent_list(x, y, CANCER)
                if cancers:
                    tx, ty = cancers[rng.integers(len(cancers))]
                    if rng.random() < cfg.kappa:
                        grid[tx, ty] = EMPTY
                        killed = True
                    # 原位增殖，受局部密度与克隆代数上限约束
                    if (rng.random() < cfg.rho
                            and gen[x, y] < cfg.max_gen
                            and len(self._adjacent_list(x, y, INDUCED)) < 6):
                        tgt = self._empty_neighbor(x, y)
                        if tgt is not None:
                            grid[tgt[0], tgt[1]] = INDUCED
                            age[tgt[0], tgt[1]] = 0
                            gen[tgt[0], tgt[1]] = gen[x, y] + 1
                if not killed and cfg.kappa_off > 0:
                    normals = self._adjacent_list(x, y, NORMAL)
                    if normals and rng.random() < cfg.kappa_off:
                        tx, ty = normals[rng.integers(len(normals))]
                        grid[tx, ty] = EMPTY
                        self.off_target += 1
                if not killed and not cancers and rng.random() < 0.7:
                    tgt = self._empty_neighbor(x, y)
                    if tgt is not None:
                        grid[x, y] = EMPTY
                        grid[tgt[0], tgt[1]] = INDUCED
                        age[tgt[0], tgt[1]] = age[x, y]
                        gen[tgt[0], tgt[1]] = gen[x, y]
                        age[x, y] = 0

        # 记录
        nc = int((grid == CANCER).sum())
        if nc > self.peak:
            self.peak = nc
        self.history.append(nc)
        if nc == 0 and self.cleared_day is None:
            self.cleared_day = day
        return nc

    def run(self):
        for day in range(1, self.cfg.max_days + 1):
            if self.step(day) == 0 and self.cleared_day is not None:
                # 清零后继续观察 20 天，用于检验是否复发
                for extra in range(day + 1, min(day + 21, self.cfg.max_days + 1)):
                    self.step(extra)
                break
        return self.result()

    def result(self):
        final = int(self.history[-1]) if self.history else int((self.grid == CANCER).sum())
        return {
            "cleared_day": self.cleared_day,
            "cleared": int(self.cleared_day is not None),
            "final_cancer": final,
            "peak_cancer": self.peak,
            "off_target": self.off_target,
            "relapsed": int(self.cleared_day is not None and final > 0),
        }


def run_one(**kwargs):
    cfg = SimConfig(**kwargs)
    out = Sim(cfg).run()
    out.update({k: v for k, v in asdict(cfg).items() if k in
                ("r", "alpha", "kappa", "kappa_off", "rho", "inflow",
                 "supply_mode", "lam_scale", "normal_frac")})
    out["seed"] = cfg.seed
    return out
