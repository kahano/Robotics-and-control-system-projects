import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import math
L1 = 100.9
L2 =  222.1
L3 = 136.2
# Oppgave 3 a & b
def jacobian_kinematics(joint_angles,joint_velocities):
    t1 = joint_angles[0]
    t2 = joint_angles[1]
    t3 = joint_angles[2]
    t1 = (t1/180.0)*np.pi
    t2 = (t2/180.0)*np.pi
    t3 = (t3/180.0)*np.pi

    q_1 = joint_velocities[0]
    q_2 = joint_velocities[1]
    q_3 = joint_velocities[2]


    z_0 = np.array([0, 0, 1])
    z_1 = np.array([np.sin(t1), -np.cos(t1), 0])
    z_2 = np.array([np.sin(t1), -np.cos(t1), 0])



    O3_0 = np.array([np.cos(t1)*(L3*np.cos(t2+t3)+L2*np.cos(t2)), np.sin(t1)*(L3*np.cos(t2+t3)+L2*np.cos(t2)), L3*np.sin(t2+t3)+L2*np.sin(t2)+L1])
    O3_1 = np.array([np.cos(t1)*(L3*np.cos(t2+t3)+L2*np.cos(t2)), np.sin(t1)*(L3*np.cos(t2+t3)+L2*np.cos(t2)), L3*np.sin(t2+t3)+L2*np.sin(t2)])
    O3_2 = np.array([np.cos(t1)*(L3*np.cos(t2+t3)), np.sin(t1)*(L3*np.cos(t2+t3)), L3*np.sin(t2+t3)])

    jv_1 = np.cross(z_0,O3_0)
    jv_2 = np.cross(z_1,O3_1)
    jv_3 = np.cross(z_2,O3_2)


    jv = np.array([jv_1,jv_2,jv_3])
    jv_transpose = jv.transpose() # her tok jeg transpose så kan det stemme med konfigurasjon jeg skal ha for videre regning av tip p hastighet.
    v = np.matmul(jv_transpose, joint_velocities) # her finner jeg cartesian hastighet av tip p
    return v



cart_velocities = jacobian_kinematics([-90,30,-45],[0.1,0.05,0.05])
print("carteisan velocities are: ",cart_velocities)
