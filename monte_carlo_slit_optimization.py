### Monte Carlo Collimation Slit Optimization ###
### Purpose: To determine the appropriate combination of detector and collimation tip slit sizes
### Author: Tyler Eastmond (teastmond@anl.gov)
### Last revision: 11/20/2023

# import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random

# user inputs
tth_deg = 20  # two theta in degrees
tth = np.radians(tth_deg) # convert two theta to radians
x1 = 50 # tip distance (mm) (from sample to tip)
x2 = 500 # detector distance (mm) (from sample to detector)
T_0 = 0.050 # tip horizontal slit size (mm)
tipV = 0.150 # tip vertical slit size (mm)
T_min = 0.025 # minimum horizontal tip slit size (mm)
H_0 = 0.1 # INITIAL GUESS: detector horizontal slit size (mm)
detV = 0.05 # detector vertical slit size (mm)
B = 0.1 # beam width (mm)
d_sample = 1.5 # sample diameter (mm)
Amax = 1.0949E-4*tth_deg**2 - 1.1883E-5*tth_deg + 3.7146E-2 # max allowable compound area (detA x tipA)
phi = np.arange(0,2*np.pi + 0.1,0.1)
spread_max = np.radians(0.160) # max angluar spread (degrees)



# number of iterations
itr = 5000
H_tweak_max = 0.0005 # maximum amount that detector slit size can be adjusted per loop (enter a value an order of magnitude less than desired)
T_tweak_max = 0.00005 # max amount that tip slit can be adjusted per loop (enter a value an order of magnitude less than desired)
H_opt_array = []
H_new_array = []
T_opt_array = []
T_new_array = []
iteration = []
# begin optimization
H_opt = H_0
T_opt = T_0
for i in range(itr):
    H_tweak = H_tweak_max*random.randint(-10,10)
    H_new = H_opt + H_tweak
    if H_new < 0:
        H_new = H_opt + abs(H_tweak)
    T_tweak = T_tweak_max*random.randint(-10,10)
    T_new = T_opt + T_tweak
    if T_new < T_min:
        T_new = T_opt + abs(T_tweak)
    # derived values
    slitA = tipV*T_new*detV*H_new
    l1 = H_new*(x2 - x1)/(T_new + H_new) # distance from focus point to detector
    l2 = x2 - x1 - l1 # distance from tip to focus point
    gamma = np.arctan(H_new/2/l1)*2 # angle of divergence of diffracted, collimated X-rays into the detector
    alpha = np.pi - tth - gamma/2 # angle of downstream 
    beta = np.pi - alpha - gamma
    a1 = np.sin(gamma/2)/np.sin(alpha)*(x1 + l2)
    a2 = np.sin(gamma/2)/np.sin(beta)*(x1 + l2)
    D0 = a1 + a2 # collimation depth with zero beam width 
    l3 = x1 + l2 - B/(2*np.sin(tth)) # distance from inboard edge of beam to focus point along the line at 2theta from centerline
    a3 = np.sin(gamma/2)*l3/np.sin(alpha) 
    a4 = np.sin(gamma/2)*l3/np.sin(beta)
    b1 = a3 + a4 # length of edge that is the lower bound of the scattering volume
    b2 = B/np.tan(beta) # extra portion of length on the upper bound of the scattering volume
    D1 = b1 + b2 # collimation depth accounting for beam width
    l4 = B/2/np.sin(tth) # distance along centerline between rotation center and inboard edge of beam
    # accept/reject proposed changes depending on whether the following criteria are met:
    if D1 <= d_sample and H_new > H_opt and slitA <= Amax and gamma <= spread_max and T_new > T_opt: # if gauge volume is within the sample volume and slit size increases
        H_opt = H_new # accept the change, update H_opt
        T_opt = T_new
    # elif D1 > d_sample and H_new < H_opt: # provides a way to reduce the slit size if the initial guess results in D1 > gauge length
    #     H_opt = H_new # accept the change, update H_opt
    #     T_opt = T_new 
    else:
        H_opt = H_opt # reject the change, retain previous H_opt
        T_opt = T_opt
    H_opt_array.append(H_opt)
    H_new_array.append(H_new)
    T_opt_array.append(T_opt)
    T_new_array.append(T_new)
    iteration.append(i)

# recalculate the derived values with the optimized detector slit size
slitA = tipV*T_opt*detV*H_opt
l1 = H_opt*(x2 - x1)/(T_opt + H_opt) # distance from focus point to detector
l2 = x2 - x1 - l1 # distance from tip to focus point
gamma = np.arctan(H_opt/2/l1)*2 # angle of divergence of diffracted, collimated X-rays into the detector
alpha = np.pi - tth - gamma/2 # angle of downstream 
beta = np.pi - alpha - gamma
a1 = np.sin(gamma/2)/np.sin(alpha)*(x1 + l2)
a2 = np.sin(gamma/2)/np.sin(beta)*(x1 + l2)
D0 = a1 + a2 # collimation depth with zero beam width 
l3 = x1 + l2 - B # distance from inboard edge of beam to focus point along the line at 2theta from centerline
a3 = np.sin(gamma/2)*l3/np.sin(alpha) 
a4 = np.sin(gamma/2)*l3/np.sin(beta)
b1 = a3 + a4 # length of edge that is the lower bound of the scattering volume
b2 = B/np.tan(beta) # extra portion of length on the upper bound of the scattering volume
D1 = b1 + b2 # collimation depth accounting for beam width
l4 = B/2/np.sin(tth) # distance along centerline between rotation center and inboard edge of beam


print("2theta = " + str(np.round(np.degrees(tth),2)))
print("D0 = " + str(np.round(D0,2)))
print("D1 = " + str(np.round(D1,2)))
print("Slit area = " + str(np.round(slitA,8)))
print("Slit max allowable area = " + str(np.round(Amax,4)))
print("Angular spread = " + str(np.round(np.degrees(gamma),4)) + " degrees")
print("Max allowable angular spread = " + str(np.round(np.degrees(spread_max),4)) + " degrees")
print("Distance between slits = " + str(np.round(l1 + l2,2)))
print("Optimized detector H slit size = " + str(round(H_opt,2)) + " mm")
print("Optimized tip H slit size = " + str(round(T_opt,4)) + " mm")

# plot the variation in detector H slit size with iteration number
plt.figure(0)
plt.plot(iteration,H_opt_array,'or:',markersize = 3)
plt.xlabel('Iteration #')
plt.ylabel('Detector horizontal slit size (mm)')
plt.show()

# plot the variation in tip H slit size with iteration number
plt.figure(0)
plt.plot(iteration,T_opt_array,'ob:',markersize = 3)
plt.xlabel('Iteration #')
plt.ylabel('Tip horizontal slit size (mm)')
plt.show()


# now plot a geometric schematic for visualization
c_rot = [0,0] # center of rotation
pt_det = [-x2*np.cos(tth),-x2*np.sin(tth)] # detector position
pt_det_out = [pt_det[0] - H_opt/2*np.sin(tth),pt_det[1] + H_opt/2*np.cos(tth)] # detector outboard slit position
pt_det_in = [pt_det[0] + H_opt/2*np.sin(tth),pt_det[1] - H_opt/2*np.cos(tth)] # detector inboard slit position
pt_focus  = [-(x1+l2)*np.cos(tth),-(x1+l2)*np.sin(tth)]
pt_tip = [pt_focus[0] + l2*np.cos(tth),pt_focus[1] + l2*np.sin(tth)] # tip center position
pt_center_inboard = [c_rot[0] - l4*np.cos(tth),c_rot[1] - B/2] # point on centerline at the inboard edge of the beam
pt_right_inboard = [pt_center_inboard[0] + a4,pt_center_inboard[1]] # rightmost point on the inboard side of the gauge volume
pt_right_outboard = [pt_right_inboard[0] + b2,pt_right_inboard[1] + B] # rightmost point on the outboard side of the gauge volume
pt_left_inboard = [pt_center_inboard[0] - a3, pt_center_inboard[1]] # leftmost point on the inboard side of the gauge volume
pt_left_outboard = [pt_left_inboard[0] + B*np.tan(alpha - np.pi/2),pt_left_inboard[1] + B] # leftmost point on the outboard side of the gauge volume


# plot the whole collimation system
plt.figure(1)
plt.axhline(y = 0,linestyle = ':',color = 'k') # center of incident beam
# plt.axhline(y = -B/2) # inboard edge of incident beam
# plt.axhline(y = B/2) # outboard edge of incident beam
plt.axhspan(B/2,-B/2,color = 'lightcoral', alpha = 0.75, lw = 0)
# top_bnd_gauge = [pt_left_inboard[0],pt_left_outboard[0],pt_right_outboard[0]]
# bot_bnd_gauge = [pt_left_inboard[1],pt_left_outboard[1],pt_right_outboard[1]]
# plt.fill_between(top_bnd_gauge,bot_bnd_gauge,-2,facecolor = 'green')
plt.plot(c_rot[0],c_rot[1],'ro') # center of rotation
plt.plot(pt_det[0],pt_det[1],'ro') # center of detector plane
plt.plot([c_rot[0],pt_det[0]],[c_rot[1],pt_det[1]],':k') # center line from sample to detector
plt.plot(pt_focus[0],pt_focus[1],'ro') # focusing point between tip and detector slits 
plt.plot(pt_det_in[0],pt_det_in[1],'bx') # detector inboard slit position
plt.plot(pt_det_out[0],pt_det_out[1],'bx') # detector outboard slit position
plt.plot(pt_tip[0],pt_tip[1],'bo') # tip center position
# plt.plot(pt_center_inboard[0],pt_center_inboard[1],'go')
plt.plot(pt_right_inboard[0],pt_right_inboard[1],'go')
plt.plot(pt_right_outboard[0],pt_right_outboard[1],'go')
plt.plot(pt_left_inboard[0],pt_left_inboard[1],'go')
plt.plot(pt_left_outboard[0],pt_left_outboard[1],'go')
plt.plot([pt_right_outboard[0],pt_det_out[0]],[pt_right_outboard[1],pt_det_out[1]],':g')
plt.plot([pt_left_outboard[0],pt_det_in[0]],[pt_left_outboard[1],pt_det_in[1]],':g')
plt.plot(d_sample/2*np.cos(phi),d_sample/2*np.sin(phi),'r')
plt.axis('equal')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.show()


# plot a zoomed-in view of the gauge volume
plt.figure(2)
plt.axhline(y = 0,linestyle = ':',color = 'k') # center of incident beam
plt.axhspan(B/2,-B/2,color = 'lightcoral', alpha = 0.75, lw = 0)
gauge_x_left = [pt_left_inboard[0], pt_left_outboard[0]]
gauge_x_right = [pt_right_inboard[0],pt_right_outboard[0]]
gauge_y = [pt_left_inboard[1],pt_left_outboard[1]]
top_bnd_gauge_x = [pt_left_inboard[0],pt_left_outboard[0],pt_right_outboard[0]]
top_bnd_gauge_y = [pt_left_inboard[1],pt_left_outboard[1],pt_right_outboard[1]]
bot_bnd_gauge_x = [pt_left_inboard[0],pt_right_inboard[0],pt_right_outboard[0]]
bot_bnd_gauge_y = [pt_left_inboard[1],pt_right_inboard[1],pt_right_outboard[1]]
plt.fill_betweenx(gauge_y,gauge_x_left,gauge_x_right, facecolor = 'green')
plt.plot(c_rot[0],c_rot[1],'ro') # center of rotation
plt.plot(pt_det[0],pt_det[1],'ro') # center of detector plane
plt.plot([c_rot[0],pt_det[0]],[c_rot[1],pt_det[1]],':k') # center line from sample to detector
plt.plot(pt_focus[0],pt_focus[1],'ro') # focusing point between tip and detector slits 
plt.plot(pt_det_in[0],pt_det_in[1],'bx') # detector inboard slit position
plt.plot(pt_det_out[0],pt_det_out[1],'bx') # detector outboard slit position
plt.plot(pt_tip[0],pt_tip[1],'bo') # tip center position
plt.plot([pt_right_outboard[0],pt_det_out[0]],[pt_right_outboard[1],pt_det_out[1]],':g')
plt.plot([pt_left_outboard[0],pt_det_in[0]],[pt_left_outboard[1],pt_det_in[1]],':g')
plt.plot(d_sample/2*np.cos(phi),d_sample/2*np.sin(phi),'r')
plt.axis('equal')
plt.xlim(1.2*np.min(d_sample/2*np.cos(phi)),1.2*np.max(d_sample/2*np.cos(phi)))
plt.ylim(1.2*np.min(d_sample/2*np.sin(phi)),1.2*np.max(d_sample/2*np.sin(phi)))
plt.title('D1 = ' + str(np.round(D1,2)) + ' mm')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')

plt.show()

