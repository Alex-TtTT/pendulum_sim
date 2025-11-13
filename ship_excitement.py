"""
ship_excitement.py

计算三点乘客（左端A、中间B、右端C）的运动学参数与兴奋度。
输入来自 pendulum_sim.py 的仿真结果（res）。
输出包含三点的高度、速度、加速度、jerk、兴奋度与累计兴奋度。
"""

import math
import numpy as np
from scipy.integrate import cumtrapz


def compute_ship_kinematics(res, L_d=20.0,
                            k_a=0.5/12.22,
                            k_j=0.3/11.38,
                            k_h=0.2/5.89):
    """
    计算船上三点（A: 左端, B: 中心, C: 右端）的运动学与兴奋度。

    参数
    ----------
    res : dict
        来自 pendulum_sim.simulate() 的仿真结果。
        至少应包含 "T", "theta", "v_t"。
    L_d : float, optional
        甲板长度（米）。
    k_a, k_j, k_h : float, optional
        兴奋度权重系数（对应 a_total, jerk, height）。

    返回
    -------
    ship_kinematics : dict
        包含各点的高度、速度、加速度、jerk、兴奋度和累计兴奋度。
    """

    # === 从仿真结果读取基础数据 ===
    L = res["params"]["L"]
    T = np.array(res["T"])
    theta = np.array(res["theta"])
    v_t = np.array(res["v_t"])  # 中心点线速度

    # === 计算几何参数 ===
    L_prime = math.sqrt(L**2 + (L_d**2) / 4)
    theta_prime = math.atan((L_d / 2) / L)

    # === 三个点的角度 ===
    theta_A = theta - theta_prime
    theta_B = theta
    theta_C = theta + theta_prime

    # === 三点的高度 ===
    h_A = L - L_prime * np.cos(theta_A)
    h_B = L - L * np.cos(theta)
    h_C = L - L_prime * np.cos(theta_C)

    # === 径向与切向分量 ===
    v_r_A = -v_t * np.sin(theta_prime) * L_prime / L
    v_r_B = np.zeros_like(v_t)
    v_r_C = v_t * np.sin(theta_prime) * L_prime / L

    v_t_A = v_t * np.cos(theta_prime) * L_prime / L
    v_t_B = v_t.copy()
    v_t_C = -v_t * np.cos(theta_prime) * L_prime / L

    # === 加速度 ===
    a_r_A = np.gradient(v_r_A, T)
    a_r_B = np.gradient(v_r_B, T)
    a_r_C = np.gradient(v_r_C, T)

    a_t_A = np.gradient(v_t_A, T)
    a_t_B = np.gradient(v_t_B, T)
    a_t_C = np.gradient(v_t_C, T)

    a_tot_A = np.sqrt(a_r_A**2 + a_t_A**2)
    a_tot_B = np.sqrt(a_r_B**2 + a_t_B**2)
    a_tot_C = np.sqrt(a_r_C**2 + a_t_C**2)

    # === jerk ===
    j_r_A = np.gradient(a_r_A, T)
    j_r_B = np.gradient(a_r_B, T)
    j_r_C = np.gradient(a_r_C, T)

    j_t_A = np.gradient(a_t_A, T)
    j_t_B = np.gradient(a_t_B, T)
    j_t_C = np.gradient(a_t_C, T)

    j_tot_A = np.sqrt(j_r_A**2 + j_t_A**2)
    j_tot_B = np.sqrt(j_r_B**2 + j_t_B**2)
    j_tot_C = np.sqrt(j_r_C**2 + j_t_C**2)

    # === 兴奋度 ===
    E_A = k_a * a_tot_A + k_j * j_tot_A + k_h * h_A
    E_B = k_a * a_tot_B + k_j * j_tot_B + k_h * h_B
    E_C = k_a * a_tot_C + k_j * j_tot_C + k_h * h_C

    # === 累积兴奋度 ===
    E_cum_A = cumtrapz(E_A, T, initial=0)
    E_cum_B = cumtrapz(E_B, T, initial=0)
    E_cum_C = cumtrapz(E_C, T, initial=0)

    # === 打包输出 ===
    ship_kinematics = {
        "T": T,
        "theta_A": theta_A,
        "theta_B": theta_B,
        "theta_C": theta_C,
        "h_A": h_A, "h_B": h_B, "h_C": h_C,
        "v_r_A": v_r_A, "v_r_B": v_r_B, "v_r_C": v_r_C,
        "v_t_A": v_t_A, "v_t_B": v_t_B, "v_t_C": v_t_C,
        "a_r_A": a_r_A, "a_r_B": a_r_B, "a_r_C": a_r_C,
        "a_t_A": a_t_A, "a_t_B": a_t_B, "a_t_C": a_t_C,
        "a_tot_A": a_tot_A, "a_tot_B": a_tot_B, "a_tot_C": a_tot_C,
        "j_r_A": j_r_A, "j_r_B": j_r_B, "j_r_C": j_r_C,
        "j_t_A": j_t_A, "j_t_B": j_t_B, "j_t_C": j_t_C,
        "j_tot_A": j_tot_A, "j_tot_B": j_tot_B, "j_tot_C": j_tot_C,
        "E_A": E_A, "E_B": E_B, "E_C": E_C,
        "E_cum_A": E_cum_A, "E_cum_B": E_cum_B, "E_cum_C": E_cum_C,
        "L_prime": L_prime,
        "theta_prime": theta_prime,
    }

    print(f"L' = {L_prime:.3f} m, θ' = {math.degrees(theta_prime):.3f}°")

    return ship_kinematics
