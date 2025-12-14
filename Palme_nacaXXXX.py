import numpy as np
from scipy.optimize import fsolve
import os

# ==========================================
# PARAMÈTRES DU PROJET
# ==========================================
R_rotor   = 1.7                   # Rayon (3.4m / 2)
R_hub     = 0.30                  # Moyeu (30 cm)
N_blades  = 3                     # 3 Pales
w         =10.5                   # Vitesse de rotation (rad/s)
v_inf     =3                      # Vitesse du vent (m/s)
TSR       =  w*R_rotor/v_inf      # TSR Global (Lambda)
Alpha_opt = 5.75                  # Angle d'attaque optimal (Pic Finesse)
CL_opt    = 1.071                 # Cl correspondant (Pic Finesse)
CD_opt    = 0.0997                # Cd correspondant (Optionnel, pour être précis dans Eq 8.4)
NACA_Code = "4418"

# Sécurités fabrication
Min_Chord = 0.10
Max_Chord = 0.50

output_dir = "Points_Turbine_Palme_NACA"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# ==========================================
# 1. FONCTIONS MATHÉMATIQUES
# ==========================================
def generate_naca_4digit(code, num_points=100):
    # Génération géométrie NACA (Même fonction que précédemment)
    m = int(code[0])/100.0
    p = int(code[1])/10.0
    t = int(code[2:])/100.0
    x = np.linspace(0, 1, num_points)
    yt = 5 * t * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)
    if m == 0:
        xu, yu, xl, yl = x, yt, x, -yt
    else:
        idx_fwd = x <= p
        idx_aft = x > p
        yc[idx_fwd] = (m/p**2)*(2*p*x[idx_fwd]-x[idx_fwd]**2)
        dyc_dx[idx_fwd] = (2*m/p**2)*(p-x[idx_fwd])
        yc[idx_aft] = (m/(1-p)**2)*((1-2*p)+2*p*x[idx_aft]-x[idx_aft]**2)
        dyc_dx[idx_aft] = (2*m/(1-p)**2)*(p-x[idx_aft])
        theta = np.arctan(dyc_dx)
        xu = x - yt*np.sin(theta)
        yu = yc + yt*np.cos(theta)
        xl = x + yt*np.sin(theta)
        yl = yc - yt*np.cos(theta)
    return np.concatenate((xu[::-1], xl[1:])), np.concatenate((yu[::-1], yl[1:]))

def solve_hansen_polynomial(x_local):
    """
    Résout l'équation 8.1 de Hansen pour trouver 'a' optimal
    16a^3 - 24a^2 + a(9 - 3x^2) - 1 + x^2 = 0
    """
    def equation(a):
        return 16*a**3 - 24*a**2 + a*(9 - 3*x_local**2) - 1 + x_local**2
    
    # On part de a=0.33 comme estimation initiale
    a_solution = fsolve(equation, 0.33)[0]
    return a_solution

x_naca, y_naca = generate_naca_4digit(NACA_Code)

# ==========================================
# 2. BOUCLE DE CALCUL
# ==========================================
r_values = np.linspace(R_hub, R_rotor, 15)

print(f"--- GÉNÉRATION SELON HANSEN EQ(8.4) ---")
print(f"{'Rayon(m)':<10} | {'TSR Loc':<8} | {'a Opt':<8} | {'Twist(°)':<10} | {'Corde(mm)':<10}")

for i, r in enumerate(r_values):
    # 1. TSR Local (x dans le livre)
    x_local = TSR * (r / R_rotor)
    
    # 2. Trouver 'a' optimal (Résolution Eq 8.1)
    a = solve_hansen_polynomial(x_local)
    
    # 3. Trouver 'a_prime' (Eq 4.38 Hansen)
    if abs(4*a - 1) < 1e-5: 
        a_prime = 0
    else:
        a_prime = (1 - 3*a) / (4*a - 1)
        
    # 4. Calcul de l'angle de flux Phi (Eq 8.3 Hansen)
    if (1 + a_prime) * x_local == 0:
        phi_rad = np.pi/2
    else:
        tan_phi = (1 - a) / ((1 + a_prime) * x_local)
        phi_rad = np.arctan(tan_phi)
    
    phi_deg = np.degrees(phi_rad)
    
    # 5. Vrillage (Twist)
    beta_deg = phi_deg - Alpha_opt
    beta_rad = np.radians(beta_deg)
    
    # 6. Corde Optimale (Eq 8.4 Hansen)
    # On utilise F=1 pour le dimensionnement initial
    F = 1.0 
    
    # Cn = Cl * cos(phi) + Cd * sin(phi) (Eq 8.5)
    Cn = CL_opt * np.cos(phi_rad) + CD_opt * np.sin(phi_rad)
    
    # Formule 8.4 : c/R = ...
    numerator = 8 * np.pi * F * a * x_local * (np.sin(phi_rad)**2)
    denominator = (1 - a) * N_blades * TSR * Cn
    
    chord_ratio = numerator / denominator # c/R
    chord = 3*chord_ratio * R_rotor         # c en mètres
    
    
 
    # 7. Construction 3D et Export
    scale = chord
    x_s = (x_naca - 0.25) * scale
    y_s = y_naca * scale
    
    c_rot, s_rot = np.cos(beta_rad), np.sin(beta_rad)
    X_sw = (x_s * c_rot - y_s * s_rot) * 1000 
    Y_sw = (x_s * s_rot + y_s * c_rot) * 1000 
    Z_sw = np.full_like(X_sw, r * 1000) 
    
    filename = f"{output_dir}/Sec_{i:02d}.txt"
    np.savetxt(filename, np.column_stack((X_sw, Y_sw, Z_sw)), fmt='%.4f', delimiter='\t')
    
    print(f"{r:<10.2f} | {x_local:<8.2f} | {a:<8.4f} | {beta_deg:<10.2f} | {chord*1000:<10.0f}")

print(f"\nTerminé ! Fichiers dans '{output_dir}'.")
