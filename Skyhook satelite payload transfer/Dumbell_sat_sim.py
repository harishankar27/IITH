import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import scienceplots

# plt.style.use('dark_background')
plt.style.use(["science","notebook"])

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",             # specify font family here
    "font.serif": ["Times New Roman"],  # specify font here
    "font.size":12})                    # specify font size here

#======== Parameters ========#
Run_ID      = "070525_01_spinning";
plot_flag   = "N";                      # "Y" string for saving
savenpz_flag= "N";                      # "Y" string for saving

G   = 6.67430e-11                       # Gravitational constant [m^3/kg^1 s^2]
M_e = 5.972e24                          # Mass of Earth [kg]
Re  = 6.37e6                            # Radius of Earth [m]
m0  = 500            
M_t = 500                               # Mass of central part of satellite (Not COM)
L1  = 1e1                               # Distance of each mass from M_t [m]
L2  = L1

#========= Function definitions =================================================================================#
def get_r1(r, phi, theta, l):
    x       = r*np.cos(theta) 
    y       = r*np.sin(theta)
    x_bob   = l*np.cos(phi)
    y_bob   = l*np.sin(phi)
    x_new   = x - x_bob
    y_new   = y - y_bob
    return np.sqrt((x_new**2) + (y_new**2))

def get_r2(r, phi, theta, l):
    # phi     = phi + np.pi
    x       = r*np.cos(theta) 
    y       = r*np.sin(theta)
    x_bob   = l*np.cos(phi)
    y_bob   = l*np.sin(phi)
    x_new   = x + x_bob
    y_new   = y + y_bob
    return np.sqrt((x_new**2) + (y_new**2))

# Define the system of ODEs
def system(t, y, t_points, m1_arr, m2_arr, M_t, L1,L2, G, M_e):
    r, r_dot, theta, theta_dot, phi, phi_dot = y
    
    r1_val      = get_r1(r, phi, theta, L1)
    r2_val      = get_r2(r, phi, theta, L2)

    m1          = np.interp(t, t_points, m1_arr)
    m2          = np.interp(t, t_points, m2_arr)
    
    I           = m1*(L1**2) + m2*(L2**2)
    M_tot       = m1 + m2 + M_t 
    mlr_val1    = m1*L1/r1_val**3
    mlr_val2    = m2*L2/r2_val**3
    #----- Equation for r_ddot ------#  
    r_ddot      = r*theta_dot**2 - (G*M_e)/(M_tot)*( (m1/r1_val**3 + m2/r2_val**3)*r + (mlr_val1 + mlr_val2)*np.cos(phi - theta) + M_t/r**2 ) 
    #----- Equation for theta_ddot ------#  
    theta_ddot  = ( 1/r )*( -2*r_dot*theta_dot + (G*M_e)/(M_tot)*(mlr_val1 - mlr_val2)*np.sin(phi-theta) )
    #----- Equation for phi_ddot ------#  
    phi_ddot    = ( G*M_e*(mlr_val2 - mlr_val1)*r*np.sin(phi-theta) )/I

    return [r_dot, r_ddot, theta_dot, theta_ddot, phi_dot, phi_ddot]


#================ Initial conditions =================================================================================================================================#
#[r, r_dot, theta, theta_dot, phi, phi_dot]

r0          = 0.5*Re + Re
theta0      = 0.0
phi0        = np.pi*0.5

# Assuming Circular orbit: Vc = r*theta_dot
Vc          = np.sqrt(G*M_e/(r0))

r_dot0      = 0.0                           # No radial acceleration
theta_dot0  = Vc/r0                         # For circular orbit
phi_dot0    = 16*theta_dot0                  

y0 = [r0,r_dot0,theta0,theta_dot0,phi0,phi_dot0]

#=========== Simulation time ====================================================================================#
ang_vel = np.sqrt(G*M_e/(r0**3))           # Same as theta_dot0
time_period = (2*np.pi/(ang_vel))

Tperiod_num = 20                             # Number of time period to be simulated
t0      = 0
tf      = Tperiod_num*time_period
t_steps = 1e6
t_eval  = np.linspace(t0, tf, int(t_steps))

m_net           = 0.001*m0;
add_sub_flag    = 1         # +1 for adding mass, -1 for removing the mass 
# Mass changes:
m1_arr                  = np.ones_like(t_eval)*m0
m1_arr[int(t_steps/4):] = np.ones_like(m1_arr[int(t_steps/4):])* (m0 + add_sub_flag* m_net);
# m1_arr[int(t_steps/4):int(t_steps/4)+int(t_steps/4)] = np.ones_like(m1_arr[int(t_steps/4):int(t_steps/4)+int(t_steps/4)])* 1.001*m0
# m1_arr[int(t_steps/4)+int(t_steps/4):] = np.ones_like(m1_arr[int(t_steps/4)+int(t_steps/4):])*m0
m2_arr                  = m1_arr
# m2_arr                    = np.ones_like(t_eval)*m0

#================ Solve the ODE ==================================================================================#
sol = solve_ivp(system, [t0, tf], y0, t_eval=t_eval, method='RK45', args=[t_eval, m1_arr, m2_arr, M_t, L1, L2, G, M_e])

#================ Post-Processing ================================================================================#
r_t     = sol.y[0,:]; r_dott        = sol.y[1,:]
theta_t = sol.y[2,:]; theta_dott    = sol.y[3,:]
phi_t   = sol.y[4,:]; phi_dott      = sol.y[5,:]

time_sol = sol.t


I = m1_arr*L1**2 + m2_arr*L2**2
# Angular momentum: H = I*omega
H_t = I * phi_dott

# Rotaional kinetic energy: RE = 1/2 I omega**2
REnergy = 0.5 * I * phi_dott**2

r_net = get_r1(r_t,r_dott,phi_t,L2)
r_net[:int(t_steps/4)] = r_net[0]
# Gravitational potential energy: GPE = SumOver_i(-G*Me*mi/Ri)
GPEnergy_Mt     = -(G*M_t*M_e / r_t)
GPEnergy_m1     = -(G*m1_arr*M_e / get_r1(r_t,r_dott,phi_t,L1))
GPEnergy_m2     = -(G*m2_arr*M_e / get_r2(r_t,r_dott,phi_t,L2))
GPEnergy_mNet   = -(G*(m_net)*M_e / r_net)
GPEnergy        = GPEnergy_Mt + GPEnergy_m1 + GPEnergy_m2

# Translational kineatic energy: TKE = 1/2 M_total * V_cm***2 
M_total = M_t+m1_arr+m2_arr
TKEnergy = 0.5 * M_total * (r_t*theta_dott)**2

# Total Energy: SumOver_allEnergy
TotalEngery = REnergy + GPEnergy + TKEnergy
#================ Plotting =====================================================================================
row,col     = 3,7;
grid        = plt.GridSpec(row,col, wspace=1, hspace=0.31);
fig_handle  = plt.figure(figsize=(12,9));


font_size       = 21
line_thickness  = 3

ax1  = fig_handle.add_subplot(grid[0,:]);
ax1.set_xlim((t_eval[0],t_eval[-1]));
ax1.set_xticks([])
ax1.set_ylabel(r"$r_{(t)}$ [m]", fontsize=font_size);

ax2  = fig_handle.add_subplot(grid[1,:]);
ax2.set_xlim((t_eval[0],t_eval[-1]));
ax2.set_xlabel("Time [s]", fontsize=font_size)
ax2.set_ylabel(r"$\theta_{(t)}^{o}$", fontsize=font_size);

ax3  = fig_handle.add_subplot(grid[2,:]);
ax3.set_xlim(13000,15500);
# ax3.set_xlim((t_eval[0],t_eval[-1]));
ax3.set_xlabel("Time [s]", fontsize=font_size)
ax3.set_ylabel(r"$\phi_{(t)}^{o}$", fontsize=font_size);

#Plotting results
ax1.plot(time_sol, r_t, linewidth=line_thickness, color="#8c6db7")
ax2.plot(time_sol, theta_t*(180/np.pi)%360, linewidth=line_thickness, color="green")
ax3.plot(time_sol, phi_t*(180/np.pi)%360, linewidth=line_thickness, color="#be491a")

if plot_flag == "Y":
    plt.savefig(f"{Run_ID}_set1.pdf",bbox_inches="tight")

#####################################################################
#####################################################################

row,col     = 3,7;
grid        = plt.GridSpec(row,col, wspace=1, hspace=0.31);
fig_handle  = plt.figure(figsize=(12,9));


ax1  = fig_handle.add_subplot(grid[0,:]);
ax1.set_xlim((t_eval[0],t_eval[-1]));
# ax1.set_xlim((2e5,3e6));
ax1.set_xticks([])
ax1.set_ylabel(r"$H_{(t)}$", fontsize=font_size);

ax2  = fig_handle.add_subplot(grid[1,:]);
ax2.set_xlim((t_eval[0],t_eval[-1]));
# ax1.set_xlim((2e5,3e6));/
ax2.set_xlabel("Time [s]", fontsize=font_size)
# ax2.set_ylabel("R.K.E. [J]", fontsize=font_size);
ax2.set_ylabel(r"$\dot{\phi}_{(t)}$ [rad/s]", fontsize=font_size);

ax3  = fig_handle.add_subplot(grid[2,:]);
# ax3.set_xlim(13000,15500);
# ax3.set_xlim(13500,30000);
ax3.set_xlim((t_eval[0],t_eval[-1]));
ax3.set_xlabel("Time [s]", fontsize=font_size)
ax3.set_ylabel("G.P.E. [J]", fontsize=font_size);


#Plotting results
ax1.plot(time_sol, H_t, linewidth=line_thickness, color="#8c6db7")
# ax2.plot(time_sol, REnergy, linewidth=line_thickness,color="green")
ax2.plot(time_sol, phi_dott, linewidth=line_thickness,color="green")
ax3.plot(time_sol, GPEnergy_m1, linewidth=line_thickness*0.5, label="m1")
ax3.plot(time_sol, GPEnergy_m2, linewidth=line_thickness*0.5, label="m2")
ax3.plot(time_sol, GPEnergy_mNet, linewidth=line_thickness*0.5, label="m_net")


plt.legend()



if plot_flag == "Y":
    plt.savefig(f"{Run_ID}_set2.pdf",bbox_inches="tight")

else: plt.show()

file_name   = f"{Run_ID}_Dumbell_Sat_data.npz"
if savenpz_flag == "Y":
    np.savez(file=file_name, L1=L1, L2=L2, M1=m1_arr, M2=m2_arr, Mt=M_t, t_eval=time_sol,
                                                        R=r_t, R_dot=r_dott,
                                                        theta=theta_t, theta_dot=theta_dott,
                                                        phi=phi_t, phi_dot=phi_dott)