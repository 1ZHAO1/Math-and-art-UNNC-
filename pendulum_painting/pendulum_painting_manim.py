from manim import *
import numpy as np
from scipy.integrate import solve_ivp
import random

# 物理常数
g = 9.81
l = 4  # 绳长，同时也是球面半径

# 初始条件
theta_0 = np.pi / 6
d_theta_0 = 0
phi_0 = 0
d_phi_0 = np.pi / 4 # 初始角速度
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
t_eval = np.linspace(t_span[0], t_span[1], 1000)

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
        # 初始化粒子群体和拖尾路径
        particles = VGroup()
        trails = VGroup()

        # 发射源为初始位置
        emission_source = np.array([x[0], y[0], z[0]])

        # 粒子发射时间间隔控制
        last_emission_time = 0
        emission_interval = 0.3  # 粒子发射的时间间隔，单位秒

        # 动态更新粒子和拖尾的函数
        def update_particles(mob, dt):
            nonlocal last_emission_time

            current_time = self.renderer.time
            current_frame = int(current_time * 25) % len(t_eval)

            # 只有当间隔达到指定时间才发射新粒子
            if current_time - last_emission_time > emission_interval:
                # 创建新粒子作为小球，并从发射源开始
                new_particle = Sphere(radius=0.1, color=BLUE)
                new_particle.move_to(emission_source)

                # 为粒子创建自己的拖尾
                new_trail = VMobject(stroke_color=BLUE, stroke_width=1)
                # 初始化拖尾的第一个点为发射源位置，确保有点可以更新
                new_trail.set_points_as_corners([emission_source])
                particles.add(new_particle)
                trails.add(new_trail)

                # 更新发射时间
                last_emission_time = current_time

            # 更新粒子位置和拖尾路径
            for i, particle in enumerate(particles):
                frame = min(current_frame, len(t_eval) - 1)

                # 移动粒子到轨迹的相应位置
                particle.move_to([x[frame], y[frame], z[frame]])

                # 更新拖尾的路径
                new_point = particle.get_center()
                
                # 检查拖尾中是否已经有点，确保拖尾可以更新
                if len(trails[i].points) > 0:
                    trails[i].add_points_as_corners([new_point])
                else:
                    # 如果没有点则初始化
                    trails[i].set_points_as_corners([new_point])

        # 添加更新器以动态更新粒子和拖尾
        particles.add_updater(update_particles)
        trails.add_updater(update_particles)

        # 添加到场景中
        self.add(particles, trails)

        # 设置相机角度
        self.set_camera_orientation(phi=75 * DEGREES, theta=45 * DEGREES)

        # 播放动画
        self.wait(20)

