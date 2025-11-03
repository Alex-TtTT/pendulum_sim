# pendulum_sim.py
from dataclasses import dataclass, asdict
import math

@dataclass
class Params:
    # --- 物理参数 ---
    m: float = 250.0
    L: float = 24.0
    g: float = 9.81

    # --- 扰动/阻尼 ---
    kappa0: float = 0.03
    c_t: float = 10.0

    # --- 驱动相关 ---
    F0: float = 6000.0
    theta_on_deg: float = 5.0           # 给力角窗（度）
    theta_target_deg: float = 80.0      # 目标峰值角（度）
    dt: float = 0.002
    t_max: float = 200.0

    # --- 命中-减推控制 ---
    n_hits: int = 10
    hit_decay: float = 0.95
    F0_min_frac: float = 0.80
    peak_eps_deg: float = 0.5

    theta_tol_deg: float = 0.2
    omega_tol_deg: float = 0.2
    quiet_time: float = 5.0

    # --- 底部反向刹车设置 ---
    brake_after_hits: int = None  # 若为 None 则等于 n_hits
    theta_brake_on_deg: float = 5.0
    F_brake0: float = 1000.0

    # --- 刹车死区/门限 ---
    omega_eps: float = 0.05
    theta_eps_deg: float = 0.3

    # --- 防卡死 ---
    stall_timeout: float = 60.0
    stall_cycles: int = 4

    # --- 初始条件 ---
    theta0_deg: float = 0.0
    omega0_deg: float = 0.0


def simulate(p: Params):
    # 角度统一转弧度
    theta_on     = math.radians(p.theta_on_deg)
    theta_target = math.radians(p.theta_target_deg)
    peak_eps     = math.radians(p.peak_eps_deg)
    theta_tol    = math.radians(p.theta_tol_deg)
    omega_tol    = math.radians(p.omega_tol_deg)
    theta_brake_on = math.radians(p.theta_brake_on_deg)
    theta_eps      = math.radians(p.theta_eps_deg)

    # 控制参数
    brake_after_hits = p.brake_after_hits if p.brake_after_hits is not None else p.n_hits

    # 初值
    t = 0.0
    theta = math.radians(p.theta0_deg)
    omega = math.radians(p.omega0_deg)

    max_abs_theta = 0.0
    drive_enabled = True
    quiet_left = p.quiet_time

    F0_curr = p.F0
    F0_min  = p.F0 * p.F0_min_frac
    hit_count = 0
    prev_omega = 0.0
    last_hit_time = -1e9
    min_hit_interval = 0.5
    miss_cycles = 0

    brake_on = False

    # 轨迹存储
    T, TH, OM, FH = [], [], [], []
    VT, AT, AR = [], [], []   # 切向速度、切向加速度、径向加速度
    H = [] 

    # 动力学（内部函数，闭包捕获 brake_on 等状态）
    def alpha(theta, omega, Fh):
        # 动态 kappa
        kappa_eff = p.kappa0 * math.sin(theta) * (1.0 if omega >= 0.0 else -1.0)

        # 底部反向刹车
        Fh_brake = 0.0
        if brake_on and abs(theta) <= theta_brake_on:
            if abs(omega) > p.omega_eps or abs(theta) > theta_eps:
                amp_theta = max(0.0, 1.0 - abs(theta)/theta_brake_on)    # 越靠近 0 刹得越强
                amp_speed = math.tanh(abs(omega)/p.omega_eps)            # 速度越大越强
                Fh_brake = - math.copysign(p.F_brake0 * amp_theta * amp_speed, omega)

        Fh_total = Fh + Fh_brake

        return ( - (p.g/p.L) * math.sin(theta)
                 + (Fh_total/(p.m*p.L)) * math.cos(theta)
                 + (kappa_eff * p.g / p.L) * math.cos(theta)
                 - (p.c_t/p.m) * omega )

    def drive(theta, omega, drive_enabled, F0_curr):
        # 达到目标极角 n 次后永久停推
        if drive_enabled and hit_count >= p.n_hits:
            drive_enabled = False

        in_window = (abs(theta) <= theta_on)
        if drive_enabled and in_window:
            amp = max(0.0, 1.0 - abs(theta)/theta_on)  # 三角窗
            Fh = math.copysign(F0_curr * amp, omega)   # 顺着当前角速度方向推
        else:
            Fh = 0.0
        return Fh, drive_enabled

    def rk4_step(theta, omega, Fh, dt):
        def f(th, om, Fh_in):
            return om, alpha(th, om, Fh_in)
        k1_th, k1_om = f(theta,                omega,                Fh)
        k2_th, k2_om = f(theta+0.5*p.dt*k1_th, omega+0.5*p.dt*k1_om, Fh)
        k3_th, k3_om = f(theta+0.5*p.dt*k2_th, omega+0.5*p.dt*k2_om, Fh)
        k4_th, k4_om = f(theta+p.dt*k3_th,     omega+p.dt*k3_om,     Fh)
        theta_new = theta + (p.dt/6.0)*(k1_th + 2*k2_th + 2*k3_th + k4_th)
        omega_new = omega + (p.dt/6.0)*(k1_om + 2*k2_om + 2*k3_om + k4_om)
        return theta_new, omega_new

    # 主循环
    while t <= p.t_max:
        Fh, drive_enabled = drive(theta, omega, drive_enabled, F0_curr)

        # 线速度/加速度
        v_t = p.L * omega
        a_t = p.L * alpha(theta, omega, Fh)
        a_r = -p.L * (omega**2)
        h = p.L * (1.0 - math.cos(theta))   


        # 记录
        T.append(t); TH.append(theta); OM.append(omega); FH.append(Fh)
        VT.append(v_t); AT.append(a_t); AR.append(a_r); H.append(h) 

        # 峰值检测：角速度过零
        if (prev_omega > 0 and omega <= 0) or (prev_omega < 0 and omega >= 0):
            peak_angle = abs(theta)
            if peak_angle >= (theta_target - peak_eps) and (t - last_hit_time) >= min_hit_interval:
                hit_count += 1
                last_hit_time = t
                # 命中后衰减推力（下限保护）
                F0_curr = max(F0_min, F0_curr * p.hit_decay)
                max_abs_theta = 0.0
                # 触发刹车阶段
                if (not brake_on) and (hit_count >= brake_after_hits):
                    drive_enabled = False
                    brake_on = True
        prev_omega = omega

        # 反卡死停机
        if drive_enabled:
            if (hit_count > 0 and (t - last_hit_time) > p.stall_timeout) or (miss_cycles >= p.stall_cycles):
                drive_enabled = False

        # RK4 推进
        theta, omega = rk4_step(theta, omega, Fh, p.dt)
        t += p.dt
        max_abs_theta = max(max_abs_theta, abs(theta))

        # 静止判据
        if not drive_enabled and (abs(theta) < theta_tol and abs(omega) < omega_tol):
            quiet_left -= p.dt
            if quiet_left <= 0.0:
                break
        else:
            quiet_left = p.quiet_time

    # 结果（角度以度返回一份，便于画图）
    theta_deg = [math.degrees(th) for th in TH]
    out = {
        "params": asdict(p),
        "T": T,
        "theta": TH,
        "theta_deg": theta_deg,
        "omega": OM,
        "Fh": FH,
        "v_t": VT,
        "a_t": AT,
        "a_r": AR,
        "h": H,   
        "meta": {
            "hit_count": hit_count,
            "brake_on": brake_on,
            "sim_end_time": t
        }
    }
    return out
