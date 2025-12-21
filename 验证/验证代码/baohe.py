import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 终局之战：饱和攻击模式 (Saturation Attack)
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

class SaturationAttackSim:
    def __init__(self):
        self.size = GRID_SIZE
        self.grid = np.zeros((self.size, self.size), dtype=int)
        self.t_cell_ages = np.zeros((self.size, self.size), dtype=int)
        
        # 终极参数：饱和攻击
        self.base_growth_rate = HIGH_MALIGNANCY_GROWTH
        self.current_growth_rate = self.base_growth_rate
        self.signal_threshold = 5
        self.use_chemo = True
        self.induction_rate = 0.9   
        self.kill_efficiency = 0.95 
        
        # 初始化
        center = self.size // 2
        self.grid[center-3:center+3, center-3:center+3] = CANCER 
        for _ in range(50): # 初始兵力加强
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
            # 化疗逻辑
            if self.use_chemo and day in chemo_schedule:
                for x in range(self.size):
                    for y in range(self.size):
                        if self.grid[x, y] == CANCER:
                            if np.random.random() < 0.6: 
                                self.grid[x, y] = EMPTY
                self.current_growth_rate = self.base_growth_rate * 0.4
                self.chemo_end_day = day + 15
            
            if self.use_chemo and day == getattr(self, 'chemo_end_day', -1):
                self.current_growth_rate = self.base_growth_rate

            history_cancer.append(np.sum(self.grid == CANCER))
            signal_active = np.sum(self.grid == CANCER) > self.signal_threshold

            new_grid = self.grid.copy()
            new_ages = self.t_cell_ages.copy()

            # --- 核心修改：饱和式造血 (模拟超高剂量动员/外源输注) ---
            # 每一天，强制支援 5 个单位的兵力！(这就是线性 vs 指数的破局关键)
            for _ in range(5): 
                rx, ry = np.random.randint(0, self.size, 2)
                if new_grid[rx, ry] == EMPTY:
                    new_grid[rx, ry] = NAIVE_T

            # 细胞行为循环
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
                        if (not signal_active and np.random.random() < 0.2) or new_ages[x, y] > 25:
                            new_grid[x, y] = EMPTY
                            new_ages[x, y] = 0
            
            self.grid = new_grid
            self.t_cell_ages = new_ages

        return history_cancer

# ==========================================
# 见证真正的归零
# ==========================================
print("Running Saturation Attack Mode...")
sim = SaturationAttackSim()
result = sim.run()

plt.figure(figsize=(10, 6))
plt.plot(result, 'g-', label='Saturation Strategy (Inflow=5/day) - True Cure', linewidth=3)
plt.title('Final Validation: Overcoming Exponential Growth with Saturation Inflow', fontsize=12)
plt.xlabel('Days')
plt.ylabel('Cancer Cells')
plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
