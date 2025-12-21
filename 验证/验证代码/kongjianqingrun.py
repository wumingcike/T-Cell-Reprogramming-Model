import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. 模拟参数设定
# ==========================================
GRID_SIZE = 50
MAX_STEPS = 200

# 设定一个“地狱模式”参数 (之前热力图里的红色死亡区)
HIGH_MALIGNANCY_GROWTH = 0.08  # 生长极快
NORMAL_INDUCTION = 0.5         # 免疫诱导正常

# 化疗参数
CHEMO_DAY = 30                 # 第30天开始化疗
CHEMO_KILL_RATE = 0.6          # 瞬间清除 60% 的癌细胞
CHEMO_GROWTH_SUPPRESSION = 0.5 # 化疗后，癌细胞生长速度减半 (持续一段时间)
CHEMO_DURATION = 20            # 抑制效果持续20天

# ==========================================
# 2. 仿真引擎 (升级版：支持化疗)
# ==========================================
EMPTY = 0
NORMAL = 1
CANCER = 2
NAIVE_T = 3
INDUCED_T = 4

class ChemoSim:
    def __init__(self, use_chemo=False):
        self.size = GRID_SIZE
        self.grid = np.zeros((self.size, self.size), dtype=int)
        self.t_cell_ages = np.zeros((self.size, self.size), dtype=int)
        self.use_chemo = use_chemo
        self.chemo_triggered = False
        
        # 初始参数
        self.base_growth_rate = HIGH_MALIGNANCY_GROWTH
        self.current_growth_rate = self.base_growth_rate
        self.induction_rate = NORMAL_INDUCTION
        self.kill_efficiency = 0.8
        self.signal_threshold = 5
        
        # 初始化
        center = self.size // 2
        self.grid[center-3:center+3, center-3:center+3] = CANCER # 初始肿瘤大一点
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
        history_t_cells = []
        
        for day in range(MAX_STEPS):
            # --- 化疗事件 ---
            if self.use_chemo and day == CHEMO_DAY:
                print(f"Day {day}: 化疗介入！炸开空间...")
                # 遍历所有格子，如果是癌细胞，按概率清除
                for x in range(self.size):
                    for y in range(self.size):
                        if self.grid[x, y] == CANCER:
                            if np.random.random() < CHEMO_KILL_RATE:
                                self.grid[x, y] = EMPTY # 变回空地！T细胞可以进来了！
                
                # 生长抑制开始
                self.current_growth_rate = self.base_growth_rate * CHEMO_GROWTH_SUPPRESSION
                self.chemo_end_day = day + CHEMO_DURATION
                self.chemo_triggered = True

            # --- 化疗效果结束 ---
            if self.use_chemo and self.chemo_triggered and day == getattr(self, 'chemo_end_day', -1):
                print(f"Day {day}: 化疗抑制结束，癌细胞恢复生长速度。")
                self.current_growth_rate = self.base_growth_rate

            # --- 常规细胞行为 ---
            new_grid = self.grid.copy()
            new_ages = self.t_cell_ages.copy()
            cancer_count = np.sum(self.grid == CANCER)
            t_cell_count = np.sum(self.grid == INDUCED_T)
            signal_active = cancer_count > self.signal_threshold
            
            # 记录数据
            history_cancer.append(cancer_count)
            history_t_cells.append(t_cell_count)

            # 极速逻辑循环
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

        return history_cancer, history_t_cells

# ==========================================
# 3. 运行对比并绘图
# ==========================================
print("正在运行对照组 (无化疗)...")
sim_control = ChemoSim(use_chemo=False)
cancer_ctrl, t_ctrl = sim_control.run()

print("正在运行实验组 (联合化疗)...")
sim_combo = ChemoSim(use_chemo=True)
cancer_combo, t_combo = sim_combo.run()

# 绘图
plt.figure(figsize=(10, 6))

# 画癌细胞曲线
plt.plot(cancer_ctrl, 'r--', label='Immune Therapy Only (Control)', linewidth=2, alpha=0.7)
plt.plot(cancer_combo, 'g-', label='Immune + Chemo Combination', linewidth=3)

# 标注化疗时间点
plt.axvline(x=CHEMO_DAY, color='k', linestyle=':', label='Chemo Intervention (Day 30)')

# 装饰图表
plt.title('Rescue Effect: Chemotherapy Creating "Space" for Immune Infiltration', fontsize=14)
plt.xlabel('Days')
plt.ylabel('Cancer Cell Count')
plt.legend()
plt.grid(True, alpha=0.3)

# 添加注释
plt.annotate('Tumor Escape (Failure)', xy=(180, cancer_ctrl[-1]), xytext=(120, 2000),
             arrowprops=dict(facecolor='red', shrink=0.05), fontsize=12, color='red')

plt.annotate('Tumor Eradicated', xy=(100, 0), xytext=(120, 500),
             arrowprops=dict(facecolor='green', shrink=0.05), fontsize=12, color='green')

plt.tight_layout()
plt.show()
