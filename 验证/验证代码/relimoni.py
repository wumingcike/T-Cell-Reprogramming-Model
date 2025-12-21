import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd
import os
import time

# ==========================================
# 1. 配置中心 (你只需调整这里)
# ==========================================
# 模拟的总次数 = 生长率测试数 * 诱导率测试数
# 这里是从 0.01 到 0.1 测 10 次，从 0.1 到 0.9 测 9 次，总共跑 90 局
GROWTH_RATES = np.linspace(0.01, 0.1, 10)  # 癌细胞恶性程度范围
INDUCTION_RATES = np.linspace(0.1, 0.9, 9) # 治疗有效率范围

GRID_SIZE = 50
MAX_STEPS = 200 # 每一局最多跑多少天
OUTPUT_FOLDER = f"Sim_Result_{time.strftime('%Y%m%d_%H%M%S')}" # 自动生成带时间的文件夹名

# 确保文件夹存在
if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)
    os.makedirs(f"{OUTPUT_FOLDER}/images") # 专门存图的子文件夹

print(f"任务启动！结果将保存在文件夹: {OUTPUT_FOLDER}")
print(f"预计执行仿真次数: {len(GROWTH_RATES) * len(INDUCTION_RATES)} 次")

# ==========================================
# 2. 核心仿真引擎 (去掉了动画，纯计算版)
# ==========================================
# 状态定义
EMPTY = 0
NORMAL = 1
CANCER = 2
NAIVE_T = 3
INDUCED_T = 4

class HeadlessSim:
    def __init__(self, growth_rate, induction_rate):
        self.size = GRID_SIZE
        self.grid = np.zeros((self.size, self.size), dtype=int)
        self.t_cell_ages = np.zeros((self.size, self.size), dtype=int)
        
        # 参数
        self.prob_growth = growth_rate
        self.prob_induction = induction_rate
        self.prob_naive_inflow = 0.05 # 固定骨髓造血能力
        self.kill_efficiency = 0.8    # 固定杀伤效率
        self.signal_threshold = 5
        
        # 初始化中心癌细胞
        center = self.size // 2
        self.grid[center-2:center+2, center-2:center+2] = CANCER
        # 初始T细胞
        for _ in range(20):
            x, y = np.random.randint(0, self.size, 2)
            if self.grid[x, y] == EMPTY:
                self.grid[x, y] = NAIVE_T

    def get_neighbors(self, x, y):
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0: continue
                nx, ny = (x + dx) % self.size, (y + dy) % self.size
                neighbors.append((nx, ny))
        return neighbors

    def run(self):
        """跑完一整局，返回最终结果"""
        history = []
        
        for day in range(MAX_STEPS):
            new_grid = self.grid.copy()
            new_ages = self.t_cell_ages.copy()
            cancer_count = np.sum(self.grid == CANCER)
            signal_active = cancer_count > self.signal_threshold
            
            # 如果癌细胞清零，提前胜利
            if cancer_count == 0:
                return "Success", day, 0, self.grid

            # 如果癌细胞占满90%，提前失败
            if cancer_count > (self.size * self.size * 0.9):
                return "Fail", day, cancer_count, self.grid

            # --- 极速逻辑循环 ---
            # 为了速度，这里不再遍历每个格子的坐标，而是使用矩阵操作或简化循环
            # 但为了保持逻辑一致性，我们还是用循环，Python处理50x50还算快
            for x in range(self.size):
                for y in range(self.size):
                    cell_type = self.grid[x, y]
                    neighbors = self.get_neighbors(x, y)
                    
                    if cell_type == CANCER:
                        if np.random.random() < self.prob_growth:
                            empties = [(nx, ny) for nx, ny in neighbors if self.grid[nx, ny] == EMPTY]
                            if empties:
                                nx, ny = empties[np.random.randint(len(empties))]
                                new_grid[nx, ny] = CANCER
                                
                    elif cell_type == NAIVE_T:
                        empties = [(nx, ny) for nx, ny in neighbors if new_grid[nx, ny] == EMPTY]
                        if empties:
                            nx, ny = empties[np.random.randint(len(empties))]
                            new_grid[nx, ny] = NAIVE_T
                            new_grid[x, y] = EMPTY
                            if signal_active and np.random.random() < self.prob_induction:
                                new_grid[nx, ny] = INDUCED_T
                                new_ages[nx, ny] = 0
                                
                    elif cell_type == INDUCED_T:
                        cancers = [(nx, ny) for nx, ny in neighbors if self.grid[nx, ny] == CANCER]
                        if cancers:
                            if np.random.random() < self.kill_efficiency:
                                tx, ty = cancers[np.random.randint(len(cancers))]
                                new_grid[tx, ty] = EMPTY
                        else:
                            empties = [(nx, ny) for nx, ny in neighbors if new_grid[nx, ny] == EMPTY]
                            if empties:
                                nx, ny = empties[np.random.randint(len(empties))]
                                new_grid[nx, ny] = INDUCED_T
                                new_grid[x, y] = EMPTY
                                new_ages[nx, ny] = new_ages[x, y] + 1
                                new_ages[x, y] = 0
                            else:
                                new_ages[x, y] += 1
                        
                        if (not signal_active and np.random.random() < 0.2) or new_ages[x, y] > 20:
                            new_grid[x, y] = EMPTY
                            new_ages[x, y] = 0

            # 骨髓造血
            if np.random.random() < self.prob_naive_inflow:
                rx, ry = np.random.randint(0, self.size, 2)
                if new_grid[rx, ry] == EMPTY:
                    new_grid[rx, ry] = NAIVE_T
            
            self.grid = new_grid
            self.t_cell_ages = new_ages
            
        return "Stalemate", MAX_STEPS, np.sum(self.grid == CANCER), self.grid

# ==========================================
# 3. 批量执行控制器
# ==========================================
results = []
cmap = ListedColormap(['white', '#eeeeee', '#ff3333', '#33cc33', '#0000ff'])

print("--- 开始批量测试 ---")

count = 0
total_runs = len(GROWTH_RATES) * len(INDUCTION_RATES)

for g_rate in GROWTH_RATES:
    for i_rate in INDUCTION_RATES:
        count += 1
        print(f"[{count}/{total_runs}] 测试参数: 生长率={g_rate:.3f}, 诱导率={i_rate:.3f} ...", end="\r")
        
        # 运行仿真
        sim = HeadlessSim(g_rate, i_rate)
        outcome, day, final_cancer, final_grid = sim.run()
        
        # 记录数据
        results.append({
            "Growth_Rate": g_rate,
            "Induction_Rate": i_rate,
            "Outcome": outcome,
            "Day_Ended": day,
            "Final_Cancer_Count": final_cancer
        })
        
        # 保存图片 (只有当结果有代表性时保存，或者全部保存)
        # 这里为了方便，我们全部保存
        plt.figure(figsize=(5, 5))
        plt.imshow(final_grid, cmap=cmap, vmin=0, vmax=4)
        plt.title(f"G={g_rate:.3f} | I={i_rate:.3f} | Result: {outcome}")
        plt.axis('off')
        plt.savefig(f"{OUTPUT_FOLDER}/images/Run_{count}_G{g_rate:.2f}_I{i_rate:.2f}_{outcome}.png")
        plt.close()

print(f"\n全部测试完成！")

# ==========================================
# 4. 生成分析报告
# ==========================================

# 1. 保存 CSV 数据表
df = pd.DataFrame(results)
df.to_csv(f"{OUTPUT_FOLDER}/Summary_Report.csv", index=False)
print(f"数据报告已保存至: {OUTPUT_FOLDER}/Summary_Report.csv")

# 2. 生成热力图 (最关键的图表)
# 我们需要把数据转换成矩阵格式：X轴是生长率，Y轴是诱导率，颜色是成功/失败
pivot_table = df.pivot(index="Induction_Rate", columns="Growth_Rate", values="Final_Cancer_Count")

plt.figure(figsize=(10, 8))
# 使用热力图：蓝色代表癌细胞少(成功)，红色代表癌细胞多(失败)
plt.imshow(pivot_table, cmap='coolwarm', origin='lower', aspect='auto')
plt.colorbar(label='Final Cancer Cell Count')

# 设置坐标轴
plt.xlabel('Cancer Growth Rate (Malignancy)')
plt.ylabel('Immune Induction Rate (Therapy Efficacy)')
plt.title('Phase Space Analysis: Therapy Success Zone vs. Failure Zone')

# 调整刻度标签
plt.xticks(range(len(GROWTH_RATES)), [f"{x:.2f}" for x in GROWTH_RATES], rotation=45)
plt.yticks(range(len(INDUCTION_RATES)), [f"{x:.2f}" for x in INDUCTION_RATES])

plt.tight_layout()
plt.savefig(f"{OUTPUT_FOLDER}/Final_Analysis_Heatmap.png")
print(f"最终分析热力图已保存至: {OUTPUT_FOLDER}/Final_Analysis_Heatmap.png")

print("--- 所有工作已完成，请查看文件夹 ---")
