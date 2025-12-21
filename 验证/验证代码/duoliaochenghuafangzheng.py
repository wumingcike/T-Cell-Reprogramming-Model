import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 升级版仿真：支持多疗程化疗
# ==========================================
GRID_SIZE = 50
MAX_STEPS = 200
HIGH_MALIGNANCY_GROWTH = 0.08
NORMAL_INDUCTION = 0.5 

EMPTY = 0
NORMAL = 1
CANCER = 2
NAIVE_T = 3
INDUCED_T = 4

class MultiCycleChemoSim:
    def __init__(self, use_chemo=False):
        self.size = GRID_SIZE
        self.grid = np.zeros((self.size, self.size), dtype=int)
        self.t_cell_ages = np.zeros((self.size, self.size), dtype=int)
        self.use_chemo = use_chemo
        
        # 初始参数
        self.base_growth_rate = HIGH_MALIGNANCY_GROWTH
        self.current_growth_rate = self.base_growth_rate
        self.induction_rate = NORMAL_INDUCTION
        self.kill_efficiency = 0.8
        self.signal_threshold = 5
        
        # 初始化中心癌细胞
        center = self.size // 2
        self.grid[center-3:center+3, center-3:center+3] = CANCER 
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
        
        # 设定化疗日程表：第30天, 60天, 90天 各打一次
        chemo_schedule = [30, 60, 90] 
        
        for day in range(MAX_STEPS):
            # --- 多疗程化疗逻辑 ---
            if self.use_chemo and day in chemo_schedule:
                # 每次清除 50%
                for x in range(self.size):
                    for y in range(self.size):
                        if self.grid[x, y] == CANCER:
                            if np.random.random() < 0.5: 
                                self.grid[x, y] = EMPTY
                # 抑制生长15天
                self.current_growth_rate = self.base_growth_rate * 0.4
                self.chemo_end_day = day + 15
                
            # 恢复生长
            if self.use_chemo and day == getattr(self, 'chemo_end_day', -1):
                self.current_growth_rate = self.base_growth_rate

            # --- 常规循环 ---
            cancer_count = np.sum(self.grid == CANCER)
            signal_active = cancer_count > self.signal_threshold
            history_cancer.append(cancer_count)

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
                        
                        if (not signal_active and np.random.random() < 0.2) or new_ages[x, y] > 20:
                            new_grid[x, y] = EMPTY
                            new_ages[x, y] = 0

            # 骨髓造血
            if np.random.random() < 0.05:
                rx, ry = np.random.randint(0, self.size, 2)
                if new_grid[rx, ry] == EMPTY:
                    new_grid[rx, ry] = NAIVE_T
            
            self.grid = new_grid
            self.t_cell_ages = new_ages

        return history_cancer

# ==========================================
# 运行对比
# ==========================================
print("Running simulation...")
sim_control = MultiCycleChemoSim(use_chemo=False)
cancer_ctrl = sim_control.run()
sim_combo = MultiCycleChemoSim(use_chemo=True)
cancer_combo = sim_combo.run()

plt.figure(figsize=(10, 6))
plt.plot(cancer_ctrl, 'r--', label='Control (No Chemo)', linewidth=2)
plt.plot(cancer_combo, 'g-', label='Combo (3-Cycle Chemo)', linewidth=3)
plt.title('Multi-Cycle Chemo Strategy: Overcoming High Malignancy (Growth Rate 0.08)')
plt.xlabel('Days')
plt.ylabel('Cancer Cells')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
