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
    theta_on_deg: float = 5.0
    theta_target_deg: float = 80.0
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
    brake_after_hits: int = None
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
    # 参数换算
    theta_on        = math.radians(p.theta_on_deg)
    theta_target    = math.radians(p.theta_target_deg)
    peak_eps        = math.radians(p.peak_eps_deg)
    theta_tol       = math.radians(p.theta_tol_deg)
    omega_tol       = math.radians(p.omega_tol_deg)
    theta_brake_on  = math.radians(p.theta_brake_on_deg)
    theta_eps       = math.radians(p.theta_eps_deg)
    brake_after_hits = p.brake_after_hits if p.brake_after_hits is not None else p.n_hits

    # 初值
    t = 0.0
    theta = math.radians(p.theta0_deg)
    omega = math.radians(p.omega0_deg)

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
    # === 平滑刹车控制状态 ===
    Fh_brake_state = 0.0      # 当前平滑刹车力
    tau_brake = 0.15          # 刹车力上升时间常数 (s)


    # 数据记录
    T, TH, OM = [], [], []
    FH, FHB, FHT = [], [], []  # ✅ Fh, Fh_brake, Fh_total
    VT, AT, AR, H = [], [], [], []

    # ---------- 内部函数 ----------
    def compute_brake_force(theta, omega):
        """计算刹车力"""
        if brake_on and abs(theta) <= theta_brake_on:
            if abs(omega) > p.omega_eps or abs(theta) > theta_eps:
                amp_theta = max(0.0, 1.0 - abs(theta)/theta_brake_on)
                amp_speed = math.tanh(abs(omega)/p.omega_eps)
                return - math.copysign(p.F_brake0 * amp_theta * amp_speed, omega)
        return 0.0

    def alpha(theta, omega, Fh):
    # 动态 kappa
        kappa_eff = p.kappa0 * math.sin(theta) * (1.0 if omega >= 0.0 else -1.0)

        # === 改这里：平滑刹车力 ===
        nonlocal Fh_brake_state
        Fh_brake_target = 0.0
        if brake_on and abs(theta) <= theta_brake_on:
            if abs(omega) > p.omega_eps or abs(theta) > theta_eps:
                amp_theta = max(0.0, 1.0 - abs(theta)/theta_brake_on)
                amp_speed = math.tanh(abs(omega)/p.omega_eps)
                Fh_brake_target = -math.copysign(p.F_brake0 * amp_theta * amp_speed, omega)

        # 平滑过渡（指数逼近）
        Fh_brake_state += (p.dt / tau_brake) * (Fh_brake_target - Fh_brake_state)
        Fh_brake = Fh_brake_state

        Fh_total = Fh + Fh_brake

        return ( - (p.g/p.L) * math.sin(theta)
                + (Fh_total/(p.m*p.L)) * math.cos(theta)
                + (kappa_eff * p.g / p.L) * math.cos(theta)
                - (p.c_t/p.m) * omega )

    def drive(theta, omega, drive_enabled, F0_curr):
        """驱动力控制"""
        if drive_enabled and hit_count >= p.n_hits:
            drive_enabled = False
        in_window = (abs(theta) <= theta_on)
        if drive_enabled and in_window:
            amp = max(0.0, 1.0 - abs(theta)/theta_on)
            Fh = math.copysign(F0_curr * amp, omega)
        else:
            Fh = 0.0
        return Fh, drive_enabled

    def rk4_step(theta, omega, Fh_total, dt):
        """RK4 积分"""
        def f(th, om, Ft):
            return om, alpha(th, om, Ft)
        k1_th, k1_om = f(theta, omega, Fh_total)
        k2_th, k2_om = f(theta+0.5*dt*k1_th, omega+0.5*dt*k1_om, Fh_total)
        k3_th, k3_om = f(theta+0.5*dt*k2_th, omega+0.5*dt*k2_om, Fh_total)
        k4_th, k4_om = f(theta+dt*k3_th, omega+dt*k3_om, Fh_total)
        theta_new = theta + (dt/6.0)*(k1_th + 2*k2_th + 2*k3_th + k4_th)
        omega_new = omega + (dt/6.0)*(k1_om + 2*k2_om + 2*k3_om + k4_om)
        return theta_new, omega_new

    # ---------- 主循环 ----------
    while t <= p.t_max:
        Fh, drive_enabled = drive(theta, omega, drive_enabled, F0_curr)
        Fh_brake = compute_brake_force(theta, omega)
        Fh_total = Fh + Fh_brake

        # 动力学
        v_t = p.L * omega
        a_t = p.L * alpha(theta, omega, Fh_total)
        a_r = -p.L * (omega**2)
        h = p.L * (1.0 - math.cos(theta))

        # 记录
        T.append(t)
        TH.append(theta)
        OM.append(omega)
        FH.append(Fh)
        FHB.append(Fh_brake)
        FHT.append(Fh_total)  
        VT.append(v_t)
        AT.append(a_t)
        AR.append(a_r)
        H.append(h)

        # 峰值检测
        if (prev_omega > 0 and omega <= 0) or (prev_omega < 0 and omega >= 0):
            peak_angle = abs(theta)
            if peak_angle >= (theta_target - peak_eps) and (t - last_hit_time) >= min_hit_interval:
                hit_count += 1
                last_hit_time = t
                F0_curr = max(F0_min, F0_curr * p.hit_decay)
                if (not brake_on) and (hit_count >= brake_after_hits):
                    drive_enabled = False
                    brake_on = True
        prev_omega = omega

        if drive_enabled:
            if (hit_count > 0 and (t - last_hit_time) > p.stall_timeout) or (miss_cycles >= p.stall_cycles):
                drive_enabled = False

        # 数值积分
        theta, omega = rk4_step(theta, omega, Fh_total, p.dt)
        t += p.dt

        # 静止判据
        if not drive_enabled and (abs(theta) < theta_tol and abs(omega) < omega_tol):
            quiet_left -= p.dt
            if quiet_left <= 0.0:
                break
        else:
            quiet_left = p.quiet_time

    # ---------- 输出 ----------
    theta_deg = [math.degrees(th) for th in TH]
    out = {
        "params": asdict(p),
        "T": T,
        "theta": TH,
        "theta_deg": theta_deg,
        "omega": OM,
        "Fh": FH,
        "Fh_brake": FHB,
        "Fh_total": FHT,   # ✅ 新增合力输出
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
