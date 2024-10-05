import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# 物理常数
g = 9.81
l = 4  # 绳长，同时也是球面半径

# 初始条件
theta_0 = np.pi / 6
d_theta_0 = 0
phi_0 = 0
d_phi_0 = np.pi / 4.2  # 初始角速度
Air_friction = 3.4 * 10**(-2)

# 方程
def equations(t, y):
    theta, d_theta, phi, d_phi = y
    dd_theta = np.sin(theta) * np.cos(theta) * d_phi**2 - (g / l) * np.sin(theta) - Air_friction * d_theta
    dd_phi = -2 * np.cos(theta) / np.sin(theta) * d_theta * d_phi - Air_friction * d_phi
    return [d_theta, dd_theta, d_phi, dd_phi]

# 时间跨度和初始条件
t_span = (0, 100)
y0 = [theta_0, d_theta_0, phi_0, d_phi_0]
t_eval = np.linspace(t_span[0], t_span[1], 5000)

# 求解微分方程
sol = solve_ivp(equations, t_span, y0, t_eval=t_eval, method='RK45')

# 提取解
theta = sol.y[0]
phi = sol.y[2]

# 转换为笛卡尔坐标
x = l * np.sin(theta) * np.cos(phi)
y = l * np.sin(theta) * np.sin(phi)
z = 1.5 - l * np.cos(theta)

# 创建图和轴
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 绘制球面
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(np.pi, np.pi*1.4, 100)
x_sphere = l * np.outer(np.sin(v), np.cos(u))
y_sphere = l * np.outer(np.sin(v), np.sin(u))
z_sphere = 1.5 + l * np.outer(np.cos(v), np.ones_like(u))
ax.plot_surface(x_sphere, y_sphere, z_sphere, color='c', alpha=0.1, linewidth=0)

# 绘制地面
ax.plot_surface(np.array([[-3.5, 3.5], [-3.5, 3.5]]), np.array([[-3.5, -3.5], [3.5, 3.5]]), np.zeros((2, 2)), color='gray', alpha=0.5)

# 设置图的边界
ax.set_xlim([-3.5, 3.5])
ax.set_ylim([-3.5, 3.5])
ax.set_zlim([-2, 5])

line, = ax.plot([], [], [], 'r-', linewidth=2)
point, = ax.plot([], [], [], 'bo')
rope, = ax.plot([], [], [], 'k-', linewidth=2)  # 绘制绳子
trace, = ax.plot([], [], [], 'g-', linewidth=2)  # 地面轨迹

trace_x, trace_y = [], []

def init():
    line.set_data([], [])
    line.set_3d_properties([])
    point.set_data([], [])
    point.set_3d_properties([])
    rope.set_data([], [])
    rope.set_3d_properties([])
    trace.set_data([], [])
    trace.set_3d_properties([])
    return line, point, rope, trace

def update(frame):
    line.set_data(x[:frame], y[:frame])
    line.set_3d_properties(z[:frame])
    point.set_data([x[frame]], [y[frame]])
    point.set_3d_properties([z[frame]])
    rope.set_data([0, x[frame]], [0, y[frame]])
    rope.set_3d_properties([1.5, z[frame]])

    # 计算交点
    if z[frame] != 1.5:
        factor = (0 - 1.5) / (z[frame] - 1.5)
        inter_x = factor * (x[frame] - 0) + 0
        inter_y = factor * (y[frame] - 0) + 0
        trace_x.append(inter_x)
        trace_y.append(inter_y)
        trace.set_data(trace_x, trace_y)
        trace.set_3d_properties([0] * len(trace_x))

    return line, point, rope, trace

# 创建动画
ani = FuncAnimation(fig, update, frames=len(t_eval), init_func=init, blit=True, interval=50)

# 显示图形
plt.show()
