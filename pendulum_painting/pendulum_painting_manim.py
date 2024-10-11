from manim import *
import numpy as np
from scipy.integrate import solve_ivp

# 物理常数
g = 9.81
l = 4  # 绳长，同时也是球面半径

# 初始条件
theta_0 = np.pi / 6
d_theta_0 = 0
phi_0 = 0
d_phi_0 = np.pi / 4.2  # 初始角速度
Air_friction = 8 * 10**(-2)

# 方程
def equations(t, y):
    theta, d_theta, phi, d_phi = y
    dd_theta = np.sin(theta) * np.cos(theta) * d_phi**2 - (g / l) * np.sin(theta) - Air_friction * d_theta
    dd_phi = -2 * np.cos(theta) / np.sin(theta) * d_theta * d_phi - Air_friction * d_phi
    return [d_theta, dd_theta, d_phi, dd_phi]

# 时间跨度和初始条件
t_span = (0, 20)
y0 = [theta_0, d_theta_0, phi_0, d_phi_0]
t_eval = np.linspace(t_span[0], t_span[1], 1000)  # 增加点的密度，提升平滑度

# 求解微分方程
sol = solve_ivp(equations, t_span, y0, t_eval=t_eval, method='RK45')

# 提取解
theta = sol.y[0]
phi = sol.y[2]

# 转换为笛卡尔坐标
x = l * np.sin(theta) * np.cos(phi)
y = l * np.sin(theta) * np.sin(phi)
z = 1.5 - l * np.cos(theta)

class SphericalPendulum(ThreeDScene):
    def construct(self):
        # 绳子的顶端位置
        top_point = np.array([0, 0, 1.5])

        # 创建摆球
        pendulum_ball = Sphere(radius=0.1, color=RED)
        pendulum_ball.move_to([x[0], y[0], z[0]])

        # 绳子
        rope = Line(top_point, pendulum_ball.get_center(), color=WHITE)

        # 创建地面轨迹
        trace = TracedPath(pendulum_ball.get_center, stroke_color=GREEN, stroke_width=2)

        # 将摆球、绳子和轨迹添加到场景中
        self.add(pendulum_ball, rope, trace)

        # 更新函数
        def update_pendulum(mob, dt):
            current_time = self.renderer.time  # 获取当前时间
            current_frame = int(current_time * 25) % len(t_eval)
            pendulum_ball.move_to([x[current_frame], y[current_frame], z[current_frame]])
            rope.put_start_and_end_on(top_point, pendulum_ball.get_center())

        pendulum_ball.add_updater(update_pendulum)
        rope.add_updater(update_pendulum)

        # 显示俯视图、侧视图和3D视图
        self.set_camera_orientation(phi=75 * DEGREES, theta=45 * DEGREES)  # 初始 3D 视图
        self.play(Create(rope), Create(pendulum_ball))
        self.wait(5)

        # 切换到俯视图
        self.move_camera(phi=0 * DEGREES, theta=90 * DEGREES, run_time=3)
        self.wait(5)

        # 切换到侧视图
        self.move_camera(phi=90 * DEGREES, theta=0 * DEGREES, run_time=3)
        self.wait(5)

        # 切换回 3D 视图
        self.move_camera(phi=75 * DEGREES, theta=45 * DEGREES, run_time=3)
        self.wait(10)
