import numpy as np
import JulianDay

#Dati
R_Earth = 6308                                                   #[km]
ALPHA_G0 = 280.46061837                                          #Pos_ang di Greenwich in dat 01/01/2000 alle 12:00 UTC
t = [15,10,2025]                                                 #Osservazione satellite in LEO
D = t(0)
M = t(1)
Y = t(2)

LONG_GS = 41.9028*np.pi/180                                      #Latitudine Roma [rad]
LAT_GS = 12.4964*np.pi/180                                       #Longitudine Roma [rad]
z = 0                                                            #quota (km)

OMEGA_Earth = ((360 + 360/365.2422)*np.pi)(24*3600)*[0, 0, 1]

rho_sez = np.array([-520.0, 2120.0, 1180.0])                     # raggio vettore[km]
dot_rho_sez = np.array([-0.45, 1.82, 0.28])                      # velocità vettore[km/s]



#Calcoliamo i JD passati dal 1/01/2000 alle 12 UTC
J2000 = JulianDay.JD(1,1,2000)
JD_0 = JulianDay.JD(D,M,Y)

d = JD_0 - J2000
T = d/36525

#Definizione di alphaG = alphaG_0 + 360.98564736629*d + 0.0003875*T^2
alpha_G = ALPHA_G0 + 360.98564736629*d + 0.0003875*T**2                                #[deg]

#Calcolo di ascensione retta (alpha) e declinazione (delta) di GS

alpha_GS = (alpha_G + LONG_GS)*np.pi/180
delta_GS = (alpha_G + LAT_GS)*np.pi/180

#Definizione di raggio vettore e velocità, rispetto al sistema topico, nel sistema di riferimento inerziale

T_matrix_SEZ_In = np.array( [[np.sin(delta_GS)*np.cos(alpha_GS), -np.sin(alpha_GS), np.cos(delta_GS)*np.cos(alpha_GS)], 
                     [np.sin(delta_GS)*np.sin(alpha_GS), np.cos(alpha_GS), np.cos(delta_GS)*np.sin(alpha_GS)],
                     [- np.cos(delta_GS), 0, np.sin(delta_GS)]])


rho_In = T_matrix_SEZ_In* rho_sez

#Calcolo raggi vettore e velocità relative al sistema di riferimento inerziale
#Raggio vettore
r_GS_In = R_Earth*np.array([np.cos(alpha_GS)*np.cos(delta_GS), np.cos(delta_GS)*np.sin(alpha_GS), np.sin(delta_GS)])

r_vect_In = r_GS_In + rho_In

# Vettore velocità
r_GS_sez = np.transpose(T_matrix_SEZ_In)*r_GS_In
omega_Earth_sez = np.transpose(T_matrix_SEZ_In)*OMEGA_Earth

v_vect_sez = dot_rho_sez + np.cross(omega_Earth_sez,r_GS_sez)
v_vect_In = T_matrix_SEZ_In*v_vect_sez





