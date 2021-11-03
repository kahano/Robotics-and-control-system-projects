import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import math

# oppgave 2 a,b,c,d
L1 = 100.9
L2 = 222.1
L3 = 136.2
"""
Her skal jeg løse alle delsoppgaver av Task 1. Jeg starter først med å sette eulervinklene i en liste "joint_angles"
men først omgjøre grader til radiner på linje 21-23 og bruke dem videre på for uttrykker av siste kolonnen av forward_kinematics
fra oblig1 der finner jeg posisjoner x,y,z for tip "p". Jeg finner også invers kinematics ved å ta posisjoner jeg fant
i deloppgave "a" og bruke dem videre for å finne eulervinkler ved geometrisk approach metode. De fire løsningene
av "elbow up og elbow down på linje 64 "

"""
def forward(joint_angles):
    teta1 = joint_angles[0]
    teta2 = joint_angles[1]
    teta3 = joint_angles[2]
    teta1 = (teta1/180.0)*np.pi
    teta2 = (teta2/180.0)*np.pi
    teta3 = (teta3/180.0)*np.pi
    a1 = ((L3*np.cos(teta1)*np.cos(teta2)*np.cos(teta3)) - (L3*np.cos(teta1)*np.sin(teta2)*np.sin(teta3)) + (L1*np.cos(teta1)*np.cos(teta2)))
    a2 = ((L3*np.sin(teta1)*np.cos(teta2)*np.cos(teta3)) - (L3*np.sin(teta1)*np.sin(teta2)*np.sin(teta3)) + (L2*np.sin(teta1)*np.cos(teta2)))
    a3 = ((L3*np.sin(teta2)*np.cos(teta3)) + (L3*np.cos(teta2)*np.sin(teta3)) + (L2*np.sin(teta2)+L1))

    cart_cord = np.array([[round(a1,4)],
                  [round(a2, 4)],
                  [round(a3, 4)]]
                  )

    return cart_cord

print("posistions x,y,z of tip p are: ",forward([270,30,-45]))
print("\n")

cart_cord = forward([-90,30,-45])

def inverse(cart_cord):
    x_03 = cart_cord[0]
    y_03 = cart_cord[1]
    z_03 = cart_cord[2]
    teta1 = math.atan2(y_03,x_03)

    teta1_1 = round((teta1/np.pi)*180,1)
    teta1_2 = 180 + teta1_1
    r2 = z_03-L1
    r1 = np.sqrt((x_03)**2 + (y_03)**2)
    r3 = np.sqrt((r1)**2+(r2)**2)
    phi2 = math.atan2(r2,r1)
    verdi1 = ((L3**2)-(L2**2)-(r3**2))/(-2*L2*r3)
    phi1 = math.atan2(np.sqrt(1-verdi1**2),verdi1)
    phi11 = math.acos(verdi1)
    teta21 = phi2+phi1
    teta22 = phi2-phi11
    teta2_1 = round((teta21/np.pi)*180,1)
    teta2_2 = round((teta22/np.pi)*180,1)
    verdi2 = (r3**2-L2**2-L3**2)/(-2*L2*L3)
    phi3 = math.atan2(np.sqrt(1-verdi2**2),verdi2)
    phi3 = round((phi3/np.pi)*180,1)
    teta3_1 = phi3 - 180
    teta3_2 = 180-phi3
    joint_angles = [[teta1_1,teta2_1,teta3_1], #four solutions
                    [teta1_1,teta2_2,teta3_2],
                    [teta1_2,180-teta2_1,-teta3_1],
                    [teta1_2,180-teta2_2,-teta3_2]]


    return joint_angles

print("Euler vinkles are: ",inverse(cart_cord))
