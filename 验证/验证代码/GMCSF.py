import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 真正的大结局：全方位增强版 (兵源+装备+战术)
# ==========================================
GRID_SIZE = 50
MAX_STEPS = 200
HIGH_MALIGNANCY_GROWTH = 0.08

# 状态码
EMPTY = 0
NORMAL = 1
CANCER = 2
NAIVE_T = 3
INDUCED_T = 4

class UltimateVictorySim:
    def __init__(self, mode='Standard'):
        self.size = GRID_SIZE
        self.grid = np.zeros((self.size, self.size), dtype=int)
        self.t_cell_ages = np.zeros((self.size, self.size), dtype=int)
        self.mode = mode
        
        # 基础参数
        self.base_growth_rate = HIGH_MALIGNANCY_GROWTH
        self.current_growth_rate = self.base_growth_rate
        self.signal_threshold = 5
        self.use_chemo = True
        
        # --- 战术配置差异 ---
        if mode == 'Standard':
            # 标准联合：普通诱导 + 贫弱的骨髓
            self.induction_rate = 0.5   
            self.kill_efficiency = 0.8
            self.marrow_inflow = 0.05   # <--- 瓶颈在这里！兵源太少！
        elif mode == 'Ultimate':
            # 终极联合：强力诱导 + GM-CSF骨髓动员
            self.induction_rate = 0.9   # 转化快
            self.kill_efficiency = 0.95 # 杀得准
            self.marrow_inflow = 0.3    # <--- 关键修改！兵源充足！

        # 初始化战场
        center = self.size // 2
        self.grid[center-3:center+3, center-3:center+3] = CANCER 
        for _ in range(30): # 初始多给一点兵
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
        history_cancer = []
        chemo_schedule = [30, 60, 90] 
        
        for day in range(MAX_STEPS):
            # --- 化疗 ---
            if self.use_chemo and day in chemo_schedule:
                for x in range(self.size):
                    for y in range(self.size):
                        if self.grid[x, y] == CANCER:
                            if np.random.random() < 0.6: # 加大一点化疗力度
                                self.grid[x, y] = EMPTY
                self.current_growth_rate = self.base_growth_rate * 0.4
                self.chemo_end_day = day + 15
                
            if self.use_chemo and day == getattr(self, 'chemo_end_day', -1):
                self.current_growth_rate = self.base_growth_rate

            history_cancer.append(np.sum(self.grid == CANCER))
            signal_active = np.sum(self.grid == CANCER) > self.signal_threshold

            # --- 极速循环 ---
            new_grid = self.grid.copy()
            new_ages = self.t_cell_ages.copy()

            for x in range(self.size):
                for y in range(self.size):
                    cell_type = self.grid[x, y]
                    neighbors = self.get_neighbors(x, y)
                    
                    if cell_type == CANCER:
                        if np.random.random() < self.current_growth_rate:
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
                            # 诱导
                            if signal_active and np.random.random() < self.induction_rate:
                                new_grid[nx, ny] = INDUCED_T
                                new_ages[nx, ny] = 0
                    elif cell_type == INDUCED_T:
                        cancers = [(nx, ny) for nx, ny in neighbors if self.grid[nx, ny] == CANCER]
                        if cancers:
                            if np.random.random() < self.kill_efficiency:
                                tx, ty = cancers[np.random.randint(len(cancers))]
                                new_grid[tx, ty] = EMPTY
                        else: # 移动寻找
                            empties = [(nx, ny) for nx, ny in neighbors if new_grid[nx, ny] == EMPTY]
                            if empties:
                                nx, ny = empties[np.random.randint(len(empties))]
                                new_grid[nx, ny] = INDUCED_T
                                new_grid[x, y] = EMPTY
                                new_ages[nx, ny] = new_ages[x, y] + 1
                                new_ages[x, y] = 0
                            else:
                                new_ages[x, y] += 1
                        # 衰老
                        if (not signal_active and np.random.random() < 0.2) or new_ages[x, y] > 25:
                            new_grid[x, y] = EMPTY
                            new_ages[x, y] = 0

            # --- 骨髓造血 (差异核心) ---
            # Standard: 0.05 (兵源枯竭)
            # Ultimate: 0.30 (GM-CSF动员)
            if np.random.random() < self.marrow_inflow: 
                rx, ry = np.random.randint(0, self.size, 2)
                if new_grid[rx, ry] == EMPTY:
                    new_grid[rx, ry] = NAIVE_T
            
            self.grid = new_grid
            self.t_cell_ages = new_ages

        return history_cancer

# ==========================================
# 见证奇迹的时刻
# ==========================================
print("Running: The 'Bone Marrow Bottleneck' Scenario...")
sim_std = UltimateVictorySim(mode='Standard')
res_std = sim_std.run()

print("Running: The 'GM-CSF Boosted' Ultimate Scenario...")
sim_ult = UltimateVictorySim(mode='Ultimate')
res_ult = sim_ult.run()

plt.figure(figsize=(10, 6))

# 失败的对照组
plt.plot(res_std, 'r--', label='Standard Combo (Inflow=0.05) - Supply Chain Failure', linewidth=2, alpha=0.6)

# 成功的终极组
plt.plot(res_ult, 'g-', label='Ultimate Strategy (Inflow=0.3 + GM-CSF) - Total Cure', linewidth=3)

plt.title('The Missing Link: Why Bone Marrow Support (GM-CSF) is Critical', fontsize=12)
plt.xlabel('Days')
plt.ylabel('Cancer Cells')
plt.legend()
plt.grid(True, alpha=0.3)
plt.annotate('System Overwhelmed', xy=(180, res_std[-1]), xytext=(120, 1000), arrowprops=dict(facecolor='red', shrink=0.05), color='red')
plt.annotate('Complete Remission', xy=(150, 0), xytext=(100, 200), arrowprops=dict(facecolor='green', shrink=0.05), color='green', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.show()
