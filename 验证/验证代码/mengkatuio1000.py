import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import os

# ==========================================
# 配置
# ==========================================
OUTPUT_DIR = "."
MC_RUNS = 1000  # 蒙特卡洛测试次数1000

# 基础参数
HIGH_MALIGNANCY_GROWTH = 0.08
EMPTY = 0
CANCER = 2
NAIVE_T = 3
INDUCED_T = 4

print(f"启动复杂度与稳健性验证系统... (目标: {MC_RUNS} 轮测试)")

# ==========================================
# 1. 极速版仿真引擎 (针对大规模测试优化)
# ==========================================
class FastSim:
    def __init__(self, size=50, mode='Saturation'):
        self.size = size
        self.grid = np.zeros((self.size, self.size), dtype=np.int8) # 优化内存
        self.mode = mode
        
        # 核心参数
        self.prob_growth = HIGH_MALIGNANCY_GROWTH
        self.prob_induction = 0.9
        self.kill_efficiency = 0.95
        self.signal_threshold = 5
        
        # 初始化
        center = self.size // 2
        # 放置肿瘤
        self.grid[center-2:center+2, center-2:center+2] = CANCER
        
        # 初始兵力
        init_t = 50 
        for _ in range(init_t):
            self.spawn_cell(NAIVE_T)

    def spawn_cell(self, cell_type):
        for _ in range(5):
            x, y = np.random.randint(0, self.size, 2)
            if self.grid[x, y] == EMPTY:
                self.grid[x, y] = cell_type
                return

    def get_neighbors_indices(self, x, y):
        # 向量化邻居查找的简化版，为了Python速度，这里用显式列表
        indices = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0: continue
                nx, ny = (x + dx) % self.size, (y + dy) % self.size
                indices.append((nx, ny))
        return indices

    def run_until_end(self, max_steps=150):
        """跑完一局，返回 (是否成功, 耗时天数, 最终癌细胞数)"""
        chemo_days = [30, 60, 90]
        chemo_end = -1
        
        for day in range(max_steps):
            # 1. 判断胜负
            cancer_count = np.sum(self.grid == CANCER)
            if cancer_count == 0:
                return True, day, 0 # 胜利
            if cancer_count > (self.size**2 * 0.95):
                return False, day, cancer_count # 失败(占满)

            # 2. 化疗逻辑
            curr_growth = self.prob_growth
            if day in chemo_days:
                # 随机清除 60%
                mask = (self.grid == CANCER)
                # 快速随机掩码
                random_mask = np.random.random(mask.shape) < 0.6
                self.grid[mask & random_mask] = EMPTY
                chemo_end = day + 15
                curr_growth *= 0.4
            elif day < chemo_end:
                curr_growth *= 0.4

            signal = cancer_count > self.signal_threshold
            new_grid = self.grid.copy()

            # 3. 饱和攻击 (每天空降5个)
            if self.mode == 'Saturation':
                for _ in range(5):
                    rx, ry = np.random.randint(0, self.size, 2)
                    if new_grid[rx, ry] == EMPTY:
                        new_grid[rx, ry] = NAIVE_T
            
            # --- 极速细胞循环 ---
            # 为了速度，我们只遍历非空格子
            # 获取所有非空坐标
            xs, ys = np.where(self.grid != EMPTY)
            # 打乱顺序模拟随机性
            indices = np.arange(len(xs))
            np.random.shuffle(indices)
            
            for i in indices:
                x, y = xs[i], ys[i]
                ctype = self.grid[x, y]
                
                # 如果这个细胞刚才被别的细胞干掉了/覆盖了，跳过
                if ctype == EMPTY: continue 
                
                # 获取邻居
                nbs = self.get_neighbors_indices(x, y)
                
                if ctype == CANCER:
                    if np.random.random() < curr_growth:
                        # 找空地繁殖
                        empties = [n for n in nbs if self.grid[n[0], n[1]] == EMPTY]
                        if empties:
                            nx, ny = empties[np.random.randint(len(empties))]
                            new_grid[nx, ny] = CANCER
                            
                elif ctype == NAIVE_T:
                    # 移动
                    empties = [n for n in nbs if new_grid[n[0], n[1]] == EMPTY]
                    if empties:
                        nx, ny = empties[np.random.randint(len(empties))]
                        new_grid[nx, ny] = NAIVE_T
                        new_grid[x, y] = EMPTY
                        # 诱导
                        if signal and np.random.random() < self.prob_induction:
                            new_grid[nx, ny] = INDUCED_T

                elif ctype == INDUCED_T:
                    # 杀伤
                    cancers = [n for n in nbs if self.grid[n[0], n[1]] == CANCER]
                    if cancers:
                        if np.random.random() < self.kill_efficiency:
                            # --- 修正点在这里 ---
                            tx, ty = cancers[np.random.randint(len(cancers))]
                            new_grid[tx, ty] = EMPTY # 直接使用 tx, ty，不要加 [0][1]
                    else:
                        # 移动
                        empties = [n for n in nbs if new_grid[n[0], n[1]] == EMPTY]
                        if empties:
                            nx, ny = empties[np.random.randint(len(empties))]
                            new_grid[nx, ny] = INDUCED_T
                            new_grid[x, y] = EMPTY
                        # 衰减
                        if not signal and np.random.random() < 0.2:
                            new_grid[x, y] = EMPTY
            
            self.grid = new_grid
            
        return False, max_steps, np.sum(self.grid == CANCER) # 超时未解决

# ==========================================
# 任务 1: 算法时间复杂度验证 (Big-O)
# ==========================================
print("\n[任务 1/2] 正在进行算法复杂度压力测试 (Grid Size Scaling)...")
sizes = [10, 20, 30, 40, 50, 60, 80]
times = []

for s in sizes:
    start_t = time.time()
    # 跑 5 局取平均，减少误差
    for _ in range(5):
        sim = FastSim(size=s, mode='Saturation')
        sim.run_until_end(max_steps=50) # 只跑50步测速度
    avg_time = (time.time() - start_t) / 5
    times.append(avg_time)
    print(f" -> 网格 {s}x{s}: 平均耗时 {avg_time:.4f} 秒/局")

# 绘图: 时间复杂度
plt.figure(figsize=(10, 6))
x = np.array(sizes)
y = np.array(times)
# 拟合曲线 (O(N^2) 因为网格点数是 N^2)
p = np.polyfit(x, y, 2)
y_fit = np.polyval(p, x)

plt.plot(x, y, 'ro', label='Measured Time')
plt.plot(x, y_fit, 'b--', label=f'Quadratic Fit ($O(L^2)$)\nAlgorithm scales linearly with cell count')
plt.title('Computational Complexity Verification', fontsize=14)
plt.xlabel('Grid Side Length (L)')
plt.ylabel('Simulation Time per 50 Steps (Seconds)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig(os.path.join(OUTPUT_DIR, "Time_Complexity.png"))
plt.close()
print(" -> 复杂度验证图已生成: Time_Complexity.png")


# ==========================================
# 任务 2: 1000 轮蒙特卡洛稳健性测试
# ==========================================
print(f"\n[任务 2/2] 正在进行 {MC_RUNS} 轮蒙特卡洛饱和攻击测试...")
print("这将模拟 1000 名患者的治疗过程，请稍候...")

outcomes = [] # 记录每一局的结果
cure_days = [] # 记录治愈的天数

start_mc = time.time()
for i in range(MC_RUNS):
    sim = FastSim(size=50, mode='Saturation')
    is_success, day, final_cancer = sim.run_until_end(max_steps=200)
    
    outcomes.append(1 if is_success else 0)
    if is_success:
        cure_days.append(day)
    
    if (i+1) % 100 == 0:
        print(f" -> 已完成 {i+1}/{MC_RUNS} 轮... 当前成功率: {sum(outcomes)/(i+1):.2%}")

total_time = time.time() - start_mc
success_rate = sum(outcomes) / MC_RUNS
avg_cure_day = np.mean(cure_days) if cure_days else 0

print(f"\n=== 测试报告 ===")
print(f"总耗时: {total_time:.2f} 秒")
print(f"总样本数 (N): {MC_RUNS}")
print(f"治疗成功率: {success_rate:.2%}")
print(f"平均治愈时间: {avg_cure_day:.1f} 天")

# 绘图: 治愈时间分布直方图
plt.figure(figsize=(10, 6))
plt.hist(cure_days, bins=30, color='green', alpha=0.7, edgecolor='black')
plt.axvline(avg_cure_day, color='red', linestyle='dashed', linewidth=2, label=f'Mean Cure Time: {avg_cure_day:.1f} Days')
plt.title(f'Monte Carlo Robustness Verification (N={MC_RUNS})\nSuccess Rate: {success_rate:.1%}', fontsize=14)
plt.xlabel('Days to Complete Remission')
plt.ylabel('Number of Patients (Simulations)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(os.path.join(OUTPUT_DIR, "Monte_Carlo_Distribution.png"))
plt.close()
print(" -> 稳健性验证图已生成: Monte_Carlo_Distribution.png")

print("\n所有验证完成！")
