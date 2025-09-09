from vpython import canvas, vector, sphere, cylinder, color, rate, textures, arrow
import numpy as np


# -----------------------
# Load data from file
# -----------------------
file_name = "Data_set/050525_02_Stability_Dumbell_Sat_data.npz"
Data = np.load(file=file_name)

t       = Data["t_eval"]
r       = Data["R"]
theta   = Data["theta"]
phi     = Data["phi"]         # Body rotation (spinning)

# -----------------------
# Satellite parameters
# -----------------------
Mt = Data["Mt"]      # central mass (visual only)
m1 = Data["M1"]      # end mass 1
m2 = Data["M2"]      # end mass 2
rod_length = 3e6#Data["L"]*1 # meters (length of dumbbell rod)

# -----------------------
# Create dual canvas
# -----------------------
# scene = canvas(title="Spinning Dumbbell Satellite Orbiting Earth", width=1000, height=700, background=color.black)
scene1 = canvas(title="Earth Orbit Overview", width=1250, height=1000, background=color.black, align='left')
scene2 = canvas(title="Zoomed COM-Tracking View", width=700, height=400, background=color.black, align='right')

# -----------------------
# Objects for scene1 (Overview)
# -----------------------
rad_Mt  = 1e5
rad_M1  = 1e5
rad_M2  = 1e5
rad_rod = 1e4

earth1 = sphere(canvas=scene1, pos=vector(0, 0, 0), radius=6.4e6, texture=textures.earth)
Mt1 = sphere(canvas=scene1, radius=rad_Mt, color=color.white, make_trail=True)
m1_sphere1 = sphere(canvas=scene1, radius=rad_M1, color=color.red)
m2_sphere1 = sphere(canvas=scene1, radius=rad_M2, color=color.cyan)
rod1 = cylinder(canvas=scene1, radius=rad_rod, color=color.yellow)

# Orbit frame arrows for scene1
r_hat_arrow1 = arrow(canvas=scene1, shaftwidth=1e3, color=color.green)
t_hat_arrow1 = arrow(canvas=scene1, shaftwidth=1e3, color=color.orange)
z_hat_arrow1 = arrow(canvas=scene1, shaftwidth=1e3, color=color.blue)

# -----------------------
# Objects for scene2 (COM-Following)
# -----------------------
earth2 = sphere(canvas=scene2, pos=vector(0, 0, 0), radius=6.4e6, texture=textures.earth, opacity=0.2)
Mt2 = sphere(canvas=scene2, radius=rad_Mt, color=color.white, make_trail=True)
m1_sphere2 = sphere(canvas=scene2, radius=rad_M1, color=color.red)
m2_sphere2 = sphere(canvas=scene2, radius=rad_M2, color=color.cyan)
rod2 = cylinder(canvas=scene2, radius=rad_rod, color=color.yellow)

# Orbit frame arrows for scene2
r_hat_arrow2 = arrow(canvas=scene2, shaftwidth=5e4, color=color.green)
t_hat_arrow2 = arrow(canvas=scene2, shaftwidth=5e4, color=color.orange)
z_hat_arrow2 = arrow(canvas=scene2, shaftwidth=5e4, color=color.blue)



# -----------------------
# Animation loop
# -----------------------
for i in range(len(t)):
    rate(10000)

    # COM of satellite
    com_pos = vector(r[i]*np.cos(theta[i]), r[i]*np.sin(theta[i]), 0)

    # Relative orientation of body
    dx = (rod_length / 2) * np.cos(phi[i])
    dy = (rod_length / 2) * np.sin(phi[i])
    rel_vec = vector(dx, dy, 0)

    # -------------------------
    # Scene 1 (Overview)
    # -------------------------
    Mt1.pos = com_pos
    m1_sphere1.pos = com_pos + rel_vec
    m2_sphere1.pos = com_pos - rel_vec
    rod1.pos = m1_sphere1.pos
    rod1.axis = m2_sphere1.pos - m1_sphere1.pos


    # -------------------------
    # Scene 2 (Zoomed-in COM-tracking)
    # -------------------------
    Mt2.pos = com_pos
    m1_sphere2.pos = com_pos + rel_vec
    m2_sphere2.pos = com_pos - rel_vec
    rod2.pos = m1_sphere2.pos
    rod2.axis = m2_sphere2.pos - m1_sphere2.pos

    # Zoomed-in camera: always above COM, looking straight down
    scene2.camera.pos = com_pos + vector(0, 0, 3e6)
    scene2.camera.axis = com_pos - scene2.camera.pos

    # -------------------------
    # Orbit Frame (unit vectors)
    # -------------------------
    r_hat = -com_pos.norm()  # Inverted radial
    v_orbit = vector(-r[i]*np.sin(theta[i])*1.1, r[i]*np.cos(theta[i])*1.1, 0)
    t_hat = v_orbit.norm()
    z_hat = vector(0, 0, 1)

    arrow_scale = 9e5

    # Scene 1 Arrows
    r_hat_arrow1.pos = com_pos
    r_hat_arrow1.axis = arrow_scale * r_hat

    t_hat_arrow1.pos = com_pos
    t_hat_arrow1.axis = arrow_scale * t_hat

    z_hat_arrow1.pos = com_pos
    z_hat_arrow1.axis = arrow_scale * z_hat

    # Scene 2 Arrows
    r_hat_arrow2.pos = com_pos
    r_hat_arrow2.axis = arrow_scale * r_hat

    t_hat_arrow2.pos = com_pos
    t_hat_arrow2.axis = arrow_scale * t_hat

    z_hat_arrow2.pos = com_pos
    z_hat_arrow2.axis = arrow_scale * z_hat
