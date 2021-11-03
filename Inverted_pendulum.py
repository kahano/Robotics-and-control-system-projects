import math
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import scipy.signal as signal
import scipy.linalg as linalg
import matplotlib.lines as lines
from scipy.integrate import solve_ivp
import matplotlib.patches as patches
from time import time


"""
I am going to implement feedback linearization control algorithms in the inverted pendulum

"""




# model parameters

L = 2 #[m] length of the bar
M =  1 #[kg] cart mass
m = 0.2 #[kg] pendulum mass
l = L/2             # CoM of uniform rod
I = (1./12) * m* (L**2)
b = 0.5           # dampending coefficient
g = -9.81
dt = 0.001
frame = 250
time_span = [0,frame*dt]
radius = 0.03
s = np.array([0,0,0,0])

# input_noise =
#Matrix weights for the cost function "LQR" they must be diagonal
Q = np.array([[1,0,0,0],
              [0,1,0,0],
              [0,0,1,0],
              [0,0,0,1]]) # weights for the outputs
Q = Q*dt

C = np.array([[1,0,0,0],
              [0,1,0,0]])

# C = C*dt

D = np.array([[0],[0]])


R = [[0.001]] # weights for the inputs
# animation = True

pendulum_Arm = lines.Line2D(s, s, color='r')
cart = patches.Rectangle(s, 0.5, 0.15, color='b')

fig = plt.figure()
ax = fig.add_subplot(111, aspect = 'equal', xlim = (-5, 5), ylim = (-1, 1), title = "Inverted Pendulum Simulation")
ax.grid()


def init():
    ax.add_patch(cart)
    ax.add_line(pendulum_Arm)
    return pendulum_Arm, cart



def model_matrix():
    A_22 = m*g/M
    A_33 = g*(M+m)/(L*M)
    A = np.array([[0,1,0,0],
                  [0,0,A_22,0],
                  [0,0,0,1],
                  [0,0,0,A_33]])
    A = np.eye(4) + dt*A # discritize

    B = np.array([[0],
                  [1/M],
                  [1],
                  [1/L*M]])

    B = dt*B #discritize

    return A, B


def LQR_optimization(): # u = -K*x
    A,B = model_matrix()
    K_r = np.matrix(linalg.solve_discrete_are(A,B,Q,R))
    K = np.matrix(linalg.inv(R)*(B.T*K_r))
    return K

def get_numpy_array_from_matrix(x):
    """
    get build-in list from matrix
    """
    return np.array(x).flatten()



def state_space(t,j):
    #model state space evaluated parameters in initial condtions
    K = LQR_optimization()
    S_x = np.sin(j[1])
    C_x = np.cos(j[1])
    dm = (m*l*C_x)**2 - (M+m) * (m*(l**2) + I)
    # nonlinear system
    s[0] = j[2]
    s[1] = j[3]
    s[2] = (1/dm)*(-(m**2)*g*(l**2)*C_x*S_x - m*l*S_x*(m*(l**2)+I)*(j[3]**2) + b*(m*(l**2) + I)*j[2]) - ((m*(l**2)+I)/dm)*u(j)
    s[3] = (1/dm)*((M+m)*m*g*l*S_x + (m**2)*(l**2)*S_x*C_x*(j[3]**2) - b*m*l*C_x*j[2]) + ((m*l*C_x)/dm)*u(j)
    return s

K = LQR_optimization()
init = np.array([-1,np.pi-0.5,-5,0])
end = np.array([1,np.pi,0,0])
u = lambda x : np.matmul(-K,(x-end))
ts = np.linspace(time_span[0],time_span[1],frame)
pendulum_state = solve_ivp(state_space, time_span, init, t_eval = ts, rtol=1e-8)



def flatten(a):
    return np.array(a).flatten()


def plot_cart(xt, theta):
    cart_w = 1.0
    cart_h = 0.5
    radius = 0.01

    c_x = np.array([-cart_w / 2.0, cart_w / 2.0, cart_w /
                   2.0, -cart_w / 2.0, -cart_w / 2.0])
    c_y = np.array([0.0, 0.0, cart_h, cart_h, 0.0])
    c_y += radius * 2.0

    c_x = c_x + xt

    b_x = np.array([0.0, L * math.sin(-theta)])
    b_x += xt
    b_y = np.array([cart_h, L * math.cos(-theta) + cart_h])
    b_y += radius * 2.0

    angles = np.arange(0.0, math.pi * 2.0, math.radians(3.0))
    o_x = np.array([radius * math.cos(a) for a in angles])
    o_y = np.array([radius * math.sin(a) for a in angles])

    rwx = np.copy(o_x) + cart_w / 4.0 + xt
    rwy = np.copy(o_y) + radius
    lwx = np.copy(o_x) - cart_w / 4.0 + xt
    lwy = np.copy(o_y) + radius

    wx = np.copy(o_x) + b_x[-1]
    wy = np.copy(o_y) + b_y[-1]

    plt.plot(flatten(c_x), flatten(c_y), "-b")
    plt.plot(flatten(b_x), flatten(b_y), "-k")
    plt.plot(flatten(rwx), flatten(rwy), "-k")
    plt.plot(flatten(lwx), flatten(lwy), "-k")
    plt.plot(flatten(wx), flatten(wy), "-k")
    plt.title("x:" + str(round(xt, 2)) + ",theta:" +
              str(round(math.degrees(theta), 2)))

    plt.axis("equal")


def animate(j):

    plt.clf()
    px = pendulum_state.y[0][j]
    theta = pendulum_state.y[1][j]
    x = [s[0] + px, s[0] + px + L * np.sin(theta)]
    y = [s[1], s[1] - L * np.cos(theta)]
    plot_cart(px,theta)
    pendulum_Arm.set_xdata(x)
    pendulum_Arm.set_ydata(y)
    cart_Pos = [s[0] + px - cart.get_width()/2, s[1] - cart.get_height()]
    cart.set_xy(cart_Pos)

    return pendulum_Arm, cart





def main():
    A,B = model_matrix()
    K = LQR_optimization()
    eigVals,eigVecs = linalg.eig(A-B*K)
    print("eigVecs")
    print(eigVecs)
    print("eigVals")
    print(eigVals)
    print("K")
    print(K)


    if animation:
        # plot_cart(px, theta)
        animate(0)
        plt.xlim([-5.0, 2.0])
        plt.pause(0.1)




if __name__ == '__main__':
    main()
