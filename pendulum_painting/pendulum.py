from manim import *
import numpy as np
from scipy.integrate import solve_ivp

class MovingSphere(ThreeDScene):
    def construct(self):
        # 物理常数
        g = 9.81  # 重力加速度
        l = 4     # 摆长
        Air_friction = 0.03  # 空气阻力系数

        # 方程定义
        def equations(t, y):
            theta, d_theta, phi, d_phi = y
            dd_theta = np.sin(theta) * np.cos(theta) * d_phi**2 - (g / l) * np.sin(theta) - Air_friction * d_theta
            dd_phi = -2 * np.cos(theta) / np.sin(theta) * d_theta * d_phi - Air_friction * d_phi
            return [d_theta, dd_theta, d_phi, dd_phi]

        # 初始条件和时间间隔
        y0 = [np.pi / 6, 0, 0, np.pi / 4]
        t_span = (0, 30)
        t_eval = np.linspace(t_span[0], t_span[1], 300)

        # 解方程
        sol = solve_ivp(equations, t_span, y0, t_eval=t_eval, method='RK45')



        # 设置更远的摄像机距离
        self.set_camera_orientation(phi=45 * DEGREES, theta=45 * DEGREES, distance=40)

        # 计算初始位置
        initial_position = np.array([
            np.sin(y0[0]) * np.cos(y0[2]),
            np.sin(y0[0]) * np.sin(y0[2]),
            -np.cos(y0[0]) + 1
        ]) * l

        # 创建小球
        sphere = Sphere(radius=0.05, color=WHITE)
        sphere.move_to(initial_position)

        # 创建轨迹对象
        trace = VMobject()
        trace.set_points_as_corners([sphere.get_center()])
        trace.start_new_path(sphere.get_center())
        trace.set_stroke(WHITE, 2)

        # 时间索引
        time_index = ValueTracker(0)

        # 获取小球的位置
        def get_sphere_position():
            idx = int(time_index.get_value())
            theta, phi = sol.y[0, idx], sol.y[2, idx]
            return np.array([
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi),
                -np.cos(theta) + 1
            ]) * l

        # 更新小球位置的函数
        def update_sphere(mob):
            new_pos = get_sphere_position()
            mob.move_to(new_pos)
            trace.add_line_to(new_pos)

        # tilie
        title = Text("Spherical Pendulum", font="Consolas").scale(0.9)


        # 方程文本
        equations_text = VGroup(
            MathTex(r"\left\{"),
            MathTex(
                r"\ddot{\theta} &= \sin(\theta) \cos(\theta) \dot{\phi}^2 - \frac{g}{l} \sin(\theta) - k \dot{\theta},\\",
                r"\ddot{\phi} &= -2 \frac{\cos(\theta)}{\sin(\theta)} \dot{\theta} \dot{\phi} - k \dot{\phi}",
            ).scale(0.5)
        ).arrange(RIGHT)



        # 播放标题
        self.add_fixed_in_frame_mobjects(title)
        self.play(Write(title),run_time=1)
        self.play(title.animate.move_to(UP * 3).scale(0.8), run_time=1)


        self.add_fixed_in_frame_mobjects(equations_text)
        self.play(FadeIn(equations_text), run_time=1)
        self.wait(1)
        self.play(equations_text.animate.to_corner(DL).scale(0.8), run_time=1)
        # 设置场景
        axes = ThreeDAxes()
        self.add(axes)
        # 渐入坐标轴和方程文本
        self.play(FadeIn(axes))



        # 在坐标轴和方程文本渐入之后添加小球和轨迹
        self.add(sphere)
        self.add(trace)
        sphere.add_updater(update_sphere)

        # 开始动画
        self.play(time_index.animate.set_value(299), run_time=10, rate_func=linear)

        # 改变摄像机视角
        self.add_fixed_in_frame_mobjects(equations_text)
        self.move_camera(phi=0 * DEGREES, theta=0 * DEGREES, distance=40, run_time=3)

        # 结束动画
        self.wait(2)
        sphere.remove_updater(update_sphere)
