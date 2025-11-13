# detect_phases_fixed.py
import numpy as np
import math

def detect_phases(res, brake_threshold=1e-3):
    """
    根据仿真结果分出三个阶段：
      1) 0 -> 第一次命中目标角
      2) 第一次命中 -> 第 n_hits 次命中
      3) 刹车开始 -> 结束
    """
    T = np.array(res["T"])
    theta = np.array(res["theta"])
    omega = np.array(res["omega"])           # 用来判定峰值：ω 过零
    Fh_brake = np.array(res["Fh_brake"])
    params = res["params"]
    meta = res["meta"]

    theta_target = math.radians(params["theta_target_deg"])
    peak_eps = math.radians(params["peak_eps_deg"])
    n_hits = params["n_hits"]

    # ---- 找“命中”的峰值：ω 过零的时刻作为峰（或谷），再看 |θ| 是否达到目标 ----
    # 过零点索引（从正到负或负到正）
    zero_cross_idx = np.where(np.signbit(omega[1:]) != np.signbit(omega[:-1]))[0] + 1
    hit_indices = [i for i in zero_cross_idx if abs(theta[i]) >= (theta_target - peak_eps)]
    hit_times = T[hit_indices] if len(hit_indices) else np.array([])

    # 阶段1结束 = 第一次命中时间（若没命中就用第一次峰值或末尾兜底）
    if len(hit_indices) > 0:
        t1 = hit_times[0]
    elif len(zero_cross_idx) > 0:
        t1 = T[zero_cross_idx[0]]
    else:
        t1 = T[-1]

    # 阶段2结束 = 第 n_hits 次命中（不足 n_hits 用最后一次命中；若没有命中，用 t1）
    if len(hit_indices) >= n_hits:
        t2 = T[hit_indices[n_hits - 1]]
    elif len(hit_indices) > 0:
        t2 = T[hit_indices[-1]]
    else:
        t2 = t1

    # # 阶段3开始 = 刹车力出现（考虑平滑后用阈值）
    # if np.any(np.abs(Fh_brake) > brake_threshold):
    #     t3 = T[np.argmax(np.abs(Fh_brake) > brake_threshold)]
    # else:
    #     # 若没有刹车就让阶段3为零长度
    #     t3 = T[-1]

    t_end = float(meta.get("sim_end_time", T[-1]))

    phases = {
        "phase1_start": 0.0,
        "phase1_end": float(t1),
        "phase1_duration": float(t1 - 0.0),

        "phase2_start": float(t1),
        "phase2_end": float(t2),
        "phase2_duration": float(t2 - t1),

        "phase3_start": float(t2),
        "phase3_end": float(t_end),
        "phase3_duration": float(t_end - t2)
    }

    # 简洁打印
    print("=== Phase timing (fixed by first-hit logic) ===")
    print(f"Phase 1: 0.00  → {phases['phase1_end']:.2f} s  (dur {phases['phase1_duration']:.2f}s)")
    print(f"Phase 2: {phases['phase2_start']:.2f} → {phases['phase2_end']:.2f} s  (dur {phases['phase2_duration']:.2f}s)")
    print(f"Phase 3: {phases['phase3_start']:.2f} → {phases['phase3_end']:.2f} s  (dur {phases['phase3_duration']:.2f}s)")

    return phases
