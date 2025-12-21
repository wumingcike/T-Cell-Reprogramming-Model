import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 终极仿真：优化协议 vs 高恶性肿瘤
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

class OptimizedWarSim:
    def __init__(self, mode='Standard'):
        self.size = GRID_SIZE
        self.grid = np.zeros((self.size, self.size), dtype=int)
        self.t_cell_ages = np.zeros((self.size, self.size), dtype=int)
        self.mode = mode
        
        # 初始参数设定
        self.base_growth_rate = HIGH_MALIGNANCY_GROWTH
        self.current_growth_rate = self.base_growth_rate
        self.signal_threshold = 5
        
        # --- 关键差异 ---
        if mode == 'Standard':
            # 之前的配置：普通免疫 + 化疗
            self.induction_rate = 0.5   # 普通诱导
            self.kill_efficiency = 0.8  # 普通杀伤
            self.use_chemo = True
        elif mode == 'Optimized':
            # 终极配置：增强免疫 + 化疗
            self.induction_rate = 0.9   # 强力诱导！(模拟高剂量IL-12)
            self.kill_efficiency = 0.95 # 精准杀伤！(模拟TRAIL增强)
            self.use_chemo = True

        # 初始化肿瘤
        center = self.size // 2
        self.grid[center-3:center+3, center-3:center+3] = CANCER 
        # 初始化T细胞
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
        history_cancer = []
        # 三疗程化疗
        chemo_schedule = [30, 60, 90] 
        
        for day in range(MAX_STEPS):
            # --- 化疗逻辑 ---
            if self.use_chemo and day in chemo_schedule:
                # 每次清除 50%
                for x in range(self.size):
                    for y in range(self.size):
                        if self.grid[x, y] == CANCER:
                            if np.random.random() < 0.5: 
                                self.grid[x, y] = EMPTY
                # 抑制生长
                self.current_growth_rate = self.base_growth_rate * 0.4
                self.chemo_end_day = day + 15
                
            if self.use_chemo and day == getattr(self, 'chemo_end_day', -1):
                self.current_growth_rate = self.base_growth_rate

            # --- 记录数据 ---
            cancer_count = np.sum(self.grid == CANCER)
            signal_active = cancer_count > self.signal_threshold
            history_cancer.append(cancer_count)

            # --- 细胞行为循环 ---
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
                            # 诱导逻辑
                            if signal_active and np.random.random() < self.induction_rate:
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
                        # 衰老逻辑
                        if (not signal_active and np.random.random() < 0.2) or new_ages[x, y] > 20:
                            new_grid[x, y] = EMPTY
                            new_ages[x, y] = 0

            if np.random.random() < 0.05: # 骨髓造血
                rx, ry = np.random.randint(0, self.size, 2)
                if new_grid[rx, ry] == EMPTY:
                    new_grid[rx, ry] = NAIVE_T
            
            self.grid = new_grid
            self.t_cell_ages = new_ages

        return history_cancer

# ==========================================
# 运行并绘图
# ==========================================
print("Running Standard Protocol (Resistance Model)...")
sim_std = OptimizedWarSim(mode='Standard')
res_std = sim_std.run()

print("Running Optimized Protocol (High Induction)...")
sim_opt = OptimizedWarSim(mode='Optimized')
res_opt = sim_opt.run()

plt.figure(figsize=(10, 6))
# 绘制之前那个耐药的曲线（作为反面教材）
plt.plot(res_std, 'r--', label='Standard Combo (Induction=0.5) - Resistance', linewidth=2, alpha=0.6)
# 绘制我们期待的完美曲线
plt.plot(res_opt, 'g-', label='Optimized Combo (Induction=0.9) - Eradication', linewidth=3)

plt.title('Final Victory: Optimized Autologous Induction vs. High-Grade Tumor', fontsize=12)
plt.xlabel('Days')
plt.ylabel('Cancer Cells')
plt.legend()
plt.grid(True, alpha=0.3)
plt.axhline(y=0, color='k', linestyle='-', alpha=0.2)
plt.show()
