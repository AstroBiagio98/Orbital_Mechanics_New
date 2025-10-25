import numpy as np
import JulianDay

#Dati
MU_EARTH = 398600.4418                                                                      #[km^3 * s^-2]
R_Earth = 6378                                                                              #[km]
ALPHA_G0 = 280.46061837                                                                     #Pos_ang di Greenwich in dat 01/01/2000 alle 12:00 UTC
t = [15,10,2025]                                                                            #Osservazione satellite in LEO
D = t[0]
M = t[1]
Y = t[2]

LAT_GS = 41.9028                                                                            #Latitudine Roma [rad]
LONG_GS = 12.4964                                                                           #Longitudine Roma [rad]
z = 0                                                                                       #quota (km)

OMEGA_Earth = ((((2*np.pi))/(86164.09)))*np.array([0, 0, 1])
print(OMEGA_Earth)

rho_sez = np.array([-520.0, 2120.0, 1180.0])                                                # raggio vettore[km]
dot_rho_sez = np.array([-0.45, 1.82, 0.28])                                                 # velocità vettore[km/s]

print(f"rho_sez = {np.linalg.norm(rho_sez)}")

print(f"dot_rho_sez = {np.linalg.norm(dot_rho_sez)}")


#Calcoliamo i JD passati dal 1/01/2000 alle 12 UTC
J2000 = JulianDay.JD(1,1,2000)
JD_0 = JulianDay.JD(D,M,Y)

d = JD_0 - J2000
T = d/36525

#print(d)
#Definizione di alphaG = alphaG_0 + 360.98564736629*d + 0.0003875*T^2
alpha_G_tot = ALPHA_G0 + 360.98564736629*d + 0.0003875*T**2                                      #[deg]
alpha_G = (alpha_G_tot/360-np.floor(alpha_G_tot/360))*360

print(f"alpha_G = {alpha_G}")
#Calcolo di ascensione retta (alpha) e declinazione (delta) di GS

alpha_GS = (alpha_G + LONG_GS)*np.pi/180
delta_GS = (LAT_GS)*np.pi/180  

print(f"alpha_GS = {alpha_GS*180/np.pi} \n delta_GS = {delta_GS*180/np.pi}")

#Definizione di raggio vettore e velocità, rispetto al sistema topico, nel sistema di riferimento inerziale

T_matrix_SEZ_In = np.array( [[np.sin(delta_GS)*np.cos(alpha_GS), -np.sin(alpha_GS), np.cos(delta_GS)*np.cos(alpha_GS)], 
                     [np.sin(delta_GS)*np.sin(alpha_GS), np.cos(alpha_GS), np.cos(delta_GS)*np.sin(alpha_GS)],
                     [- np.cos(delta_GS), 0, np.sin(delta_GS)]])

#print("\n\n\n")

#print(f"Matrice di rotazione = {T_matrix_SEZ_In}")

#print("\n\n\n")
rho_In = np.dot(T_matrix_SEZ_In,rho_sez)
#print(f"rho_IN = {rho_In}")

#Calcolo raggi vettore e velocità relative al sistema di riferimento inerziale
#Raggio vettore

r_GS_In = R_Earth*(np.transpose(np.array([np.cos(alpha_GS)*np.cos(delta_GS), np.cos(delta_GS)*np.sin(alpha_GS), np.sin(delta_GS)])))

print(f"{np.linalg.norm(r_GS_In)}")
print(f"r_GS_In = {r_GS_In}")

r_vect_In = r_GS_In + rho_In
#print("\n")
#print(f"r_vect_In = {r_vect_In}")
#print(f"r_vect_In = {np.linalg.norm(r_vect_In)}")
#print("\n")

# Vettore velocità
r_sez = np.dot(np.transpose(T_matrix_SEZ_In),r_vect_In)
omega_Earth_sez = np.dot(np.linalg.matrix_transpose(T_matrix_SEZ_In),OMEGA_Earth)

v_vect_sez = np.transpose(dot_rho_sez) + np.cross(omega_Earth_sez,r_sez)
v_vect_In = np.dot(T_matrix_SEZ_In,np.transpose(v_vect_sez))

#print("\n")
#print(f"v_vect_In = {v_vect_In}")
#print(f"v_vect_In = {np.linalg.norm(v_vect_In)}")
#print("\n")

#Determinazione dei parametri orbitali classici
I_VERS = np.array([1,0,0])
J_VERS = np.array([0,1,0])
K_VERS = np.array([0,0,1])


h_vect = np.cross(r_vect_In, v_vect_In)
h_mod = np.linalg.norm(h_vect)
h_vers = h_vect/h_mod

#print(f"modulo_h_versore = {np.linalg.norm(h_vers)}")
n_vect = np.cross(np.transpose(K_VERS), h_vect)
n_vers = n_vect/np.linalg.norm(n_vect)
#print(f"modulo_n_versore = {np.linalg.norm(n_vers)}")

#print(f"n_vers = {n_vers}, h_vers = {h_vers}")
      
raan = np.acos(np.dot(I_VERS, n_vers))*180/np.pi

ecc_vect = np.cross(v_vect_In, h_vect)/MU_EARTH - r_vect_In/np.linalg.norm(r_vect_In)
ecc_mod = np.linalg.norm(ecc_vect)
ecc_vers = ecc_vect/ecc_mod
#print(f"mod_ecc_vers = {np.linalg.norm(ecc_vers)}")

omega = np.acos(np.dot(n_vers, ecc_vers))*180/np.pi
 
inc = np.acos(np.dot(K_VERS, h_vers))*180/np.pi

a = h_mod**2/MU_EARTH*(1/(1 - ecc_mod**2))

true_anomaly = np.acos(np.dot(r_vect_In,ecc_vect)/(np.linalg.norm(r_vect_In)*(np.linalg.norm(ecc_vect))))*180/np.pi

print(f"semiasse maggiore = {a} \n eccentricità = {ecc_mod} \n RAAN = {raan} \n omega = {omega} \n inclinazione = {inc} \n anomalia vera = {true_anomaly}")

Energy = (np.linalg.norm(v_vect_In)**2)/2 - MU_EARTH/(np.linalg.norm(r_vect_In))
print(f"Energy = {Energy}")

#print(f"a = {- MU_EARTH/(2*Energy)}")





