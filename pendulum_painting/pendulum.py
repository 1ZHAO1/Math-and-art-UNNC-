from manim import *
import numpy as np
from scipy.integrate import solve_ivp

class MovingSphere(ThreeDScene):
    def construct(self):
        # Physical constants
        g = 9.81  # Gravity
        l = 4     # Length of pendulum
        Air_friction = 3 * 10**(-2) # Air friction coefficient

        # Equations definition
        def equations(t, y):
            theta, d_theta, phi, d_phi = y
            dd_theta = np.sin(theta) * np.cos(theta) * d_phi**2 - (g / l) * np.sin(theta) - Air_friction * d_theta
            dd_phi = -2 * np.cos(theta) / np.sin(theta) * d_theta * d_phi - Air_friction * d_phi
            return [d_theta, dd_theta, d_phi, dd_phi]

        # Initial conditions and time interval
        y0 = [np.pi / 6, 0, 0, np.pi / 4]
        t_span = (0, 300)
        t_eval = np.linspace(t_span[0], t_span[1], 36000)

        # Solve the equations
        sol = solve_ivp(equations, t_span, y0, t_eval=t_eval, method='RK45')

        # Calculate initial position
        initial_position = np.array([
            np.sin(y0[0]) * np.cos(y0[2]),
            np.sin(y0[0]) * np.sin(y0[2]),
            -np.cos(y0[0]) + 1
        ]) * l
        # Set further camera distance
        self.set_camera_orientation(phi=45 * DEGREES, theta=45 * DEGREES, distance=40)
        # Create sphere
        sphere = Sphere(radius=0.05, color=WHITE)
        sphere.move_to(initial_position)

        # Create trace object
        trace = VMobject()
        trace.set_points_as_corners([sphere.get_center()])
        trace.start_new_path(sphere.get_center())
        trace.set_stroke(WHITE, 2)

        # Time index
        time_index = ValueTracker(0)

        # Function to get sphere's position
        def get_sphere_position():
            idx = int(time_index.get_value())
            theta, phi = sol.y[0, idx], sol.y[2, idx]
            return np.array([
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi),
                -np.cos(theta) + 1
            ]) * l

        # Function to update sphere's position
        def update_sphere(mob):
            new_pos = get_sphere_position()
            mob.move_to(new_pos)
            trace.add_line_to(new_pos)

        # Title
        title = Text("Spherical Pendulum", font="Consolas").scale(0.9)

        # Equations text
        equations_text = VGroup(
            MathTex(r"\left\{"),
            MathTex(
                r"\ddot{\theta} &= \sin(\theta) \cos(\theta) \dot{\phi}^2 - \frac{g}{l} \sin(\theta) - k \dot{\theta},\\",
                r"\ddot{\phi} &= -2 \frac{\cos(\theta)}{\sin(\theta)} \dot{\theta} \dot{\phi} - k \dot{\phi}",
            ).scale(0.5)
        ).arrange(RIGHT)

        # Play title
        self.add_fixed_in_frame_mobjects(title)
        self.play(Write(title), run_time=1)
        self.play(title.animate.move_to(UP * 3).scale(0.8), run_time=1)

        self.add_fixed_in_frame_mobjects(equations_text)
        self.play(FadeIn(equations_text), run_time=1)
        self.wait(1)
        self.play(equations_text.animate.to_corner(DL).scale(0.8), run_time=1)
        # Set the scene
        axes = ThreeDAxes()
        self.add(axes)
        # Fade in axes and equations text
        self.play(FadeIn(axes))

        # Add sphere and trace after fading in axes and equations text
        self.add(sphere)
        self.add(trace)
        sphere.add_updater(update_sphere)

        # Start animation
        self.play(time_index.animate.set_value(7999), run_time=20, rate_func=linear)

        # Change camera view
        self.add_fixed_in_frame_mobjects(equations_text)
        self.move_camera(phi=0 * DEGREES, theta=0 * DEGREES, distance=40, run_time=3)

        # End animation
        self.wait(2)
        sphere.remove_updater(update_sphere)
