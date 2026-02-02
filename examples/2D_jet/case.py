#!/usr/bin/env python3
import math
import json
import cantera as ct
ctfile = "h2o2.yaml"
sol_L = ct.Solution(ctfile)
sol_L.TPX = 400, 101325, "O2:1,N2:1"

pA = 101325
rhoA = sol_L.density
gam = 1.4
c_l = math.sqrt(1.4 * pA / rhoA)

pS = pA
velS = 0.02 
rhoS = rhoA

leng = 1e-2
Ny = 250
Nx = Ny
dx = leng / Nx
L = 10/1000 # 10 mm
time_end = 2 * leng / velS
cfl = 0.8

dt = cfl * dx / c_l
Nt = int(time_end / dt)

eps = 1e-5

# Configuring case dictionary
case = {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -L,
            "x_domain%end": L,
            "y_domain%beg": -L/2,
            "y_domain%end": L/2,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": dt, 
            "t_step_start": 0,
            "t_step_stop": 4000,
            "t_step_save": 400,
            "t_step_print": 1,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
             "chemistry": "T" ,
             "chem_params%diffusion": "F",
             "chem_params%reactions": "F",
  #           "elliptic_smoothing" : "T",
 #            "elliptic_smoothing_iters": 200,
            "bc_x%beg": -8,
            "bc_x%end": -8,
            "bc_y%beg": -15,
            "bc_y%end": -15,
            "num_bc_patches": 1,
            "patch_bc(1)%dir": 2,
            "patch_bc(1)%loc": 1,
            "patch_bc(1)%geometry": 1,
            "patch_bc(1)%type": -3,
            "patch_bc(1)%centroid(1)": 0.0,
            "patch_bc(1)%length(1)": 0.15*L,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 10 * L,
            "patch_icpp(1)%length_y": 10 * L,
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": 101325,
            "patch_icpp(1)%alpha_rho(1)": rhoA,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%hcid": 283,
            "patch_icpp(2)%x_centroid": 0,
            "patch_icpp(2)%y_centroid": L/2,
            "patch_icpp(2)%length_x":  0.15*L,
            "patch_icpp(2)%length_y": 0.6*L,
            "patch_icpp(2)%vel(2)": -velS,
            "patch_icpp(2)%vel(1)": 0,
            "patch_icpp(2)%pres": pS,
            "patch_icpp(2)%alpha_rho(1)": rhoA,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
            "cantera_file": ctfile,
        }
for i in range(len(sol_L.Y)):
      case[f"patch_icpp(1)%Y({i+1})"] = sol_L.Y[i]
      case[f"patch_icpp(2)%Y({i+1})"] = sol_L.Y[i]
    
if __name__ == "__main__":
    print(json.dumps(case))
