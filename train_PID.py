import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import math
import random
import matplotlib.animation as animation

# from sim_plot_train.py import plot


############################################
# import inputs

forsok = 3
skraa_plan_vinkel = np.pi/6 * 0
plan_weight = 100 #[kg]
g = 10 # gravitasjon [m/s^2]

#Tunning constantene
K_P = 300
K_D = 300
K_I = 10


global_forsok = forsok

# generisk x_referanse til fallende cube

def set_x_referanse(skraa_plan_vinkel):
    vannrettetPos_random = random.uniform(0,120) # equal probability for distribuition
    loddrettetPos_random = random.uniform(20+120*np.tan(skraa_plan_vinkel)+6.5,40+120*np.tan(skraa_plan_vinkel)+6.5)

    return vannrettetPos_random,loddrettetPos_random

dt  = 0.02
t0 = 0
t_end = 5 # sekunder

total_tid = np.arange(t0,t_end+dt,dt) # array for  tid
F_g = -plan_weight*g # gravitasjonskraft
display_jernbane = np.zeros((forsok,len(total_tid)))
v_gjernbane = np.zeros((forsok,len(total_tid)))
h_gjernbane = np.zeros((forsok,len(total_tid)))
x_pos_tog = np.zeros((forsok,len(total_tid)))
y_pos_tog = np.zeros((forsok,len(total_tid)))

tang_kraft = F_g * math.sin(skraa_plan_vinkel)
error = np.zeros((forsok,len(total_tid)))
error_dot = np.zeros((forsok,len(total_tid)))
error_Integral = np.zeros((forsok,len(total_tid)))

x_pos_cube = np.zeros((forsok,len(total_tid)))
y_pos_cube = np.zeros((forsok,len(total_tid)))

init_posx = 120
init_posy = 120*math.tan(skraa_plan_vinkel) + 6.5
init_disp = math.sqrt(init_posx**2+init_posy**2)**(0.5)
start_fart = 0
start_aks = 0

init_globalt_posx = init_posx
backup = np.ones(forsok)
trials = forsok

while(forsok > 0):  # bestemmer hvor mange ganger cube faller ned
    x_pos_cube_ref = set_x_referanse(skraa_plan_vinkel)[0]
    y_pos_cube_ref = set_x_referanse(skraa_plan_vinkel)[1]
    antallganger = trials - forsok
    x_pos_cube[antallganger] = x_pos_cube_ref
    y_pos_cube[antallganger] = y_pos_cube_ref - g/2*total_tid**2
    ferdig = False
    delta = 1

    for i in range(1,len(total_tid)): ## Implementere PID Kontroller

        if(i == 1):
            display_jernbane[antallganger][0] = init_disp
            x_pos_tog[antallganger][0] = init_posx
            y_pos_tog[antallganger][0] = init_posy
            v_gjernbane[antallganger][0] = start_fart
            h_gjernbane[antallganger][0] = start_aks

        error[antallganger][i-1] = x_pos_cube_ref - x_pos_tog[antallganger][i-1]

        if(i > 1):
            error_dot[antallganger][i-1] = (error[antallganger][i-1] -  error[antallganger][i-2])/dt
            error_Integral[antallganger][i-1] = error_Integral[antallganger][i-2] + (error[antallganger][i-2] + error[antallganger][i-1])/2*dt

        if(i == len(total_tid)-1):
            error[antallganger][-1] = error_dot[antallganger][-2]
            error_dot[antallganger][-1] = error_dot[antallganger][-2]
            error_Integral[antallganger][-1] = error_Integral[antallganger][-2]

        F_a = K_P*error[antallganger][i-1] + K_D*error_dot[antallganger][i-1] + K_I*error_Integral[antallganger][i-1] # applied force
        F_net = F_a + F_g
        h_gjernbane[antallganger][i] = F_net / plan_weight
        v_gjernbane[antallganger][i] = v_gjernbane[antallganger][i-1] + (h_gjernbane[antallganger][i-1]+ v_gjernbane[antallganger][i])/2*dt
        display_jernbane[antallganger][i] = display_jernbane[antallganger][i-1] + (v_gjernbane[antallganger][i-1] + v_gjernbane[antallganger][i])/2*dt
        x_pos_tog[antallganger][i] = display_jernbane[antallganger][i]*math.cos(skraa_plan_vinkel)
        y_pos_tog[antallganger][i] = display_jernbane[antallganger][i] *math.cos(skraa_plan_vinkel) + 6.5

        ## prøv å fange cubet

        if(x_pos_tog[antallganger][i] - 5 < x_pos_cube[antallganger][0] + 3 and x_pos_tog[antallganger][i] + 5 > x_pos_cube[antallganger][1]-3)or ferdig == True:
            if(y_pos_tog[antallganger][i] + 3 < y_pos_cube[antallganger][i] - 2 and y_pos_tog[antallganger][i] + 8 > y_pos_cube[antallganger][i]+ 2 )or ferdig == True :
                ferdig = True

                if delta == 1 :
                    forskjell = x_pos_tog[antallganger][i] - x_pos_cube[antallganger][i]
                    delta = 0
                x_pos_cube[antallganger][i] = x_pos_tog[antallganger][i] - forskjell
                y_pos_cube[antallganger][i] = y_pos_tog[antallganger][i] + 5

    init_disp = display_jernbane[antallganger][-1]
    init_posx = x_pos_tog[antallganger][-1] + v_gjernbane[antallganger][-1]*math.cos(skraa_plan_vinkel)*dt
    init_posy = y_pos_tog[antallganger][-1] +  v_gjernbane[antallganger][-1]*math.sin(skraa_plan_vinkel)*dt
    start_fart = v_gjernbane[antallganger][-1]
    start_aks = h_gjernbane[antallganger][-1]
    backup[antallganger] = delta
    forsok = forsok - 1
################################################################### slutt av simulering##########################################################





##################################### animation ########################################3


len_t=len(total_tid)
frame_amount=len(total_tid)*trials
def update_plot(num):

    platform.set_data([x_pos_tog[int(num/len_t)][num-int(num/len_t)*len_t]-3.1,\
    x_pos_tog[int(num/len_t)][num-int(num/len_t)*len_t]+3.1],\
    [y_pos_tog[int(num/len_t)][num-int(num/len_t)*len_t],\
    y_pos_tog[int(num/len_t)][num-int(num/len_t)*len_t]])

    cube.set_data([x_pos_cube[int(num/len_t)][num-int(num/len_t)*len_t]-1,\
    x_pos_cube[int(num/len_t)][num-int(num/len_t)*len_t]+1],\
    [y_pos_cube[int(num/len_t)][num-int(num/len_t)*len_t],\
    y_pos_cube[int(num/len_t)][num-int(num/len_t)*len_t]])

    if global_forsok*len_t==num+1 and num>0: # All attempts must be successful
        if sum(backup)==0:
            success.set_text('congratulation', 'YOU DID IT !' )
        else:
            again.set_text('DONT GIVE UP! YOU CAN DO IT')

    display_jernbane_f.set_data(total_tid[0:(num-int(num/len_t)*len_t)],
        display_jernbane[int(num/len_t)][0:(num-int(num/len_t)*len_t)])

    v_gjernbane_f.set_data(total_tid[0:(num-int(num/len_t)*len_t)],
        v_gjernbane[int(num/len_t)][0:(num-int(num/len_t)*len_t)])

    h_gjernbane_f.set_data(total_tid[0:(num-int(num/len_t)*len_t)],
        h_gjernbane[int(num/len_t)][0:(num-int(num/len_t)*len_t)])

    error_f.set_data(total_tid[0:(num-int(num/len_t)*len_t)],
        error[int(num/len_t)][0:(num-int(num/len_t)*len_t)])

    error_dot_f.set_data(total_tid[0:(num-int(num/len_t)*len_t)],
        error_dot[int(num/len_t)][0:(num-int(num/len_t)*len_t)])

    error_Integral_f.set_data(total_tid[0:(num-int(num/len_t)*len_t)],
        error_Integral[int(num/len_t)][0:(num-int(num/len_t)*len_t)])

    return display_jernbane_f,v_gjernbane_f,h_gjernbane_f,error_f,error_dot_f,error_Integral_f,platform,cube,success,again

fig=plt.figure(figsize=(16,9),dpi=120,facecolor=(0.8,0.8,0.8))
gs=gridspec.GridSpec(4,3)

# Create main window
ax_main=fig.add_subplot(gs[0:3,0:2],facecolor=(0.9,0.9,0.9))
plt.xlim(0,init_globalt_posx)
plt.ylim(0,init_globalt_posx)
plt.xticks(np.arange(0,init_globalt_posx+1,10))
plt.yticks(np.arange(0,init_globalt_posx+1,10))
plt.grid(True)

copyright=ax_main.text(0,122,'© Mohamed Siraj A Ali',size=12)

rail=ax_main.plot([0,init_globalt_posx],[5,init_globalt_posx*np.tan(skraa_plan_vinkel)+5],'k',linewidth=6)
platform,=ax_main.plot([],[],'b',linewidth=18)
cube,=ax_main.plot([],[],'k',linewidth=14)

bbox_props_success=dict(boxstyle='square',fc=(0.9,0.9,0.9),ec='g',lw='1')
success=ax_main.text(40,60,'',size='20',color='g',bbox=bbox_props_success)

bbox_props_again=dict(boxstyle='square',fc=(0.9,0.9,0.9),ec='r',lw='1')
again=ax_main.text(30,60,'',size='20',color='r',bbox=bbox_props_again)

# Plot windows
ax1v=fig.add_subplot(gs[0,2],facecolor=(0.9,0.9,0.9))
display_jernbane_f,=ax1v.plot([],[],'-b',linewidth=2,label='displ. on rails [m]')
plt.xlim(t0,t_end)
plt.ylim(np.min(display_jernbane)-abs(np.min(display_jernbane))*0.1,np.max(display_jernbane)+abs(np.max(display_jernbane))*0.1)
plt.grid(True)
plt.legend(loc='lower left',fontsize='small')

ax2v=fig.add_subplot(gs[1,2],facecolor=(0.9,0.9,0.9))
v_gjernbane_f,=ax2v.plot([],[],'-b',linewidth=2,label='velocity on rails [m/s]')
plt.xlim(t0,t_end)
plt.ylim(np.min(v_gjernbane)-abs(np.min(v_gjernbane))*0.1,np.max(v_gjernbane)+abs(np.max(v_gjernbane))*0.1)
plt.grid(True)
plt.legend(loc='lower left',fontsize='small')

ax3v=fig.add_subplot(gs[2,2],facecolor=(0.9,0.9,0.9))
h_gjernbane_f,=ax3v.plot([],[],'-b',linewidth=2,label='accel. on rails [m/s^2] = F_net/m_platf.')
plt.xlim(t0,t_end)
plt.ylim(np.min(h_gjernbane)-abs(np.min(h_gjernbane))*0.1,np.max(h_gjernbane)+abs(np.max(h_gjernbane))*0.1)
plt.grid(True)
plt.legend(loc='lower left',fontsize='small')

ax1h=fig.add_subplot(gs[3,0],facecolor=(0.9,0.9,0.9))
error_f,=ax1h.plot([],[],'-b',linewidth=2,label='horizontal error [m]')
plt.xlim(t0,t_end)
plt.ylim(np.min(error)-abs(np.min(error))*0.1,np.max(error)+abs(np.max(error))*0.1)
plt.grid(True)
plt.legend(loc='lower left',fontsize='small')

ax2h=fig.add_subplot(gs[3,1],facecolor=(0.9,0.9,0.9))
error_dot_f,=ax2h.plot([],[],'-b',linewidth=2,label='change of horiz. error [m/s]')
plt.xlim(t0,t_end)
plt.ylim(np.min(error_dot)-abs(np.min(error_dot))*0.1,np.max(error_dot)+abs(np.max(error_dot))*0.1)
plt.grid(True)
plt.legend(loc='lower left',fontsize='small')

ax3h=fig.add_subplot(gs[3,2],facecolor=(0.9,0.9,0.9))
error_Integral_f,=ax3h.plot([],[],'-b',linewidth=2,label='sum of horiz. error [m*s]')
plt.xlim(t0,t_end)
plt.ylim(np.min(error_Integral)-abs(np.min(error_Integral))*0.1,np.max(error_Integral)+abs(np.max(error_Integral))*0.1)
plt.grid(True)
plt.legend(loc='lower left',fontsize='small')

pid_ani=animation.FuncAnimation(fig,update_plot,
    frames=frame_amount,interval=20,repeat=False,blit=True)
plt.show()
