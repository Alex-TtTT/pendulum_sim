本项目将同一套代码拆分为 **核心仿真模块**（Python 脚本）与 **可视化/示例**（Jupyter Notebook）两部分：

- **`pendulum_sim.py`**
  - 提供参数数据类 `Params` 与主函数 `simulate(params)`；
  - 只做数值积分和数据产出，不进行绘图，便于你在任何环境/工程中复用。
- **`demo.ipynb`**
  - 演示如何设置参数、调用 `simulate` 并用 `matplotlib` 绘制速度/加速度/角位移曲线；
  - 作为实验记录与可视化入口。



---

## 1. 模型简介



系统是长度为 L、质量为 m 的摆，考虑重力、气动（或阻力）项与 **受控水平推力**（在小角窗内，随角速度方向施加），以及 **底部反向刹车**（在靠近竖直附近按角度与角速度自适应地给反向力）。

动力学采用非线性方程：

θ¨ = -(g/L)·sinθ + (Fₕ/(mL))·cosθ + (κ₀g/L)·sinθ·cosθ·sgn(θ˙) - (cₜ/m)·θ˙

其中：
- g 为重力加速度；
- Fₕ 为水平推力；
- κ₀ 为扰动幅度系数；
- sgn(θ˙) 表示角速度的符号函数；
- cₜ 为切向阻尼系数。

该方程与代码中 `alpha(theta, omega, Fh)` 的计算逻辑完全一致。

---

**命中-减推控制**：当达到目标极角 \\(\\theta_{\\text{target}}\\)（允许容差）计一次“命中”；每次命中后将推力幅值按 `hit_decay` 衰减（不低于下限 `F0_min_frac*F0`）。当累计命中次数达到 `n_hits` 时：**永久停推**，并进入**底部反向刹车**阶段以迅速稳定。

---

## 2. 文件结构

```
.
├── pendulum_sim.py              # 仿真核心：Params, simulate()
├── demo.ipynb                   # 可视化示例：参数设置 + 调用 + 画图
└── README.md                    # 你正在看的文档
```

---

## 3. 依赖与环境

- Python 3.9+（建议）
- 仅仿真：标准库 `math`, `dataclasses`
- 可视化（在 Notebook 中）：`matplotlib`

安装：

```bash
pip install matplotlib
```

---

## 4. 快速开始

### 4.1 在 Jupyter Notebook 中

1. 将 `pendulum_sim.py` 与 `demo.ipynb` 放在同一目录；
2. 打开 `demo.ipynb`，运行所有单元即可得到 4 张常用曲线：
   - 切向线速度 `v_t(t)`
   - 切向线加速度 `a_t(t)`
   - 角位移 `θ(t)`（含目标角 ±target 参考线）
   - 径向加速度 `a_r(t) = -L·ω²`

### 4.2 在 Python 脚本中调用

```python
from pendulum_sim import Params, simulate

p = Params(
    m=250.0, L=24.0, g=9.81,
    kappa0=0.03, c_t=10.0,
    F0=6000.0,
    theta_on_deg=5.0,
    theta_target_deg=80.0,
    dt=0.002, t_max=200.0,
    n_hits=10, hit_decay=0.95, F0_min_frac=0.80,
    peak_eps_deg=0.5,
    theta_tol_deg=0.2, omega_tol_deg=0.2,
    quiet_time=5.0,
    brake_after_hits=None,          # 默认与 n_hits 相同，或可手动指定
    F_brake0=1000.0,
    theta_brake_on_deg=5.0,
    omega_eps=0.05,
    theta_eps_deg=0.3,
    stall_timeout=60.0,
    stall_cycles=4
)

out = simulate(p)
T = out["t"]
theta_deg = out["theta_deg"]
# 其余键见下文“输出说明”
```

---

## 5. 主要参数（`Params`）

> 以下字段名可能与你最初脚本中的中文注释一一对应；此处统一为可读的英文字段并提供“*_deg”便捷角度单位。

| 参数名 | 含义 | 典型值/建议 |
|---|---|---|
| `m` | 质量 (kg) | 250.0 |
| `L` | 摆长 (m) | 24.0 |
| `g` | 重力加速度 (m/s²) | 9.81 |
| `kappa0` | 动态支持力切向分量系数基准 | 0.03 |
| `c_t` | 切向阻尼系数 | 10.0 |
| `F0` | 初始推力幅值 (N) | 6000.0 |
| `theta_on_deg` | 小角窗阈值（仅在 |θ|≤该角时给力） | 5.0° |
| `theta_target_deg` | 目标极角（命中阈值） | 80.0° |
| `dt` | 积分步长 (s) | 0.002 |
| `t_max` | 最大仿真时长 (s) | 200.0 |
| `n_hits` | 需要命中的次数（达到目标极角） | 10 |
| `hit_decay` | 每次命中后推力衰减因子 | 0.95 |
| `F0_min_frac` | 推力下限（相对初始 `F0`） | 0.80 |
| `peak_eps_deg` | 命中角容差 | 0.5° |
| `theta_tol_deg` | 停机判据的角度阈值 | 0.2° |
| `omega_tol_deg` | 停机判据的角速度阈值 | 0.2°/s |
| `quiet_time` | 满足小角小速需持续的时间 (s) | 5.0 |
| `brake_after_hits` | 触发底部刹车的命中次数（默认与 `n_hits` 相同） | - |
| `F_brake0` | 底部反向刹车峰值 (N) | 1000.0 |
| `theta_brake_on_deg` | 刹车只在该小角窗内生效 | 5.0° |
| `omega_eps` | 近零速门限 (rad/s) | 0.05 |
| `theta_eps_deg` | 近零角门限 (deg) | 0.3° |
| `stall_timeout` | 停推防卡死：自上次命起超过该秒数则停推 | 60 s |
| `stall_cycles` | 停推防卡死：连续错过峰值的次数阈值 | 4 |

---

## 6. 输出说明（`simulate(params)` 返回）

`simulate` 返回一个 `dict`，包含时间序列与元信息：

- **时间序列**
  - `t`：时间 (s)
  - `theta`：角位移 (rad)
  - `theta_deg`：角位移 (deg)
  - `omega`：角速度 (rad/s)
  - `Fh`：实际施加的水平推力 (N)，含命中-减推与底部刹车叠加
  - `v_t`：切向线速度 (m/s)，\\(v_t = L\\,\\omega\\)
  - `a_t`：切向线加速度 (m/s²)，来自模型右端项
  - `a_r`：径向加速度 (m/s²)，\\(a_r = -L\\,\\omega^2\\)
- **元信息**
  - `meta.hit_count`：最终命中次数
  - `meta.brake_on`：是否进入过底部刹车阶段
  - `meta.sim_end_time`：实际仿真停止时刻（可能早于 `t_max`）

> 注：Notebook 中演示了如何从 `out` 中取 `T/VT/AT/AR/θ(t)` 等并绘图。

---

## 7. 可视化（Notebook 中示例）

- `Tangential linear velocity vs time`
- `Tangential linear acceleration vs time`
- `Angular displacement vs time`（含 ±target 参考线）
- `Radial acceleration vs time (a_r = -L·ω²)`

你可按需新增：
- 推力、高度时序 `Fh(t)`；
- 命中事件标记（在峰值处散点/竖线）；

---



