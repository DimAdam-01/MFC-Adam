#!/usr/bin/env python3
import json
import argparse
import math

import cantera as ct

parser = argparse.ArgumentParser(
    prog="nD_inert_shocktube",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default='{}', metavar="DICT",
                    help="MFC's toolchain's internal state.")
parser.add_argument("--no-chem", dest='chemistry', default=True, action="store_false",
                    help="Disable chemstry.")

args = parser.parse_args()

ctfile    = 'grigri.yaml'
sol_L     = ct.Solution(ctfile)
sol_L.TPX =  300,  101325, 'H:1'

L    = 0.03
Nx   = 300
Ny   = 300
dx   = L / Nx
dt   = 3e-9
Tend = 0.9e-4

NT         = int(Tend / dt)
SAVE_COUNT = 100
NS         = 4000

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": 0,
    "x_domain%end": +L,
    "y_domain%beg": 0,
    "y_domain%end": +L,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "cyl_coord": "F",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": NT,
    "t_step_save": NS,
    "t_step_print": SAVE_COUNT,
    # Simulation Algorithm
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mixture_err": "F",
    "mpp_lim": "F",
    "time_stepper": 3,
    "avg_state": 1,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_Re_flux": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "bc_x%beg": -8,
    "bc_x%end": -8,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "num_patches": 1,
    "num_fluids": 1,
    "viscous": "F",
    'chemistry'                    : 'T' if not args.chemistry else 'T',
    'chem_params%diffusion'        : 'F',
    'chem_params%reactions'        : 'F',
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # Fluid Parameters (Heavy Gas)
    "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0.0e00,
    #"fluid_pp(1)%Re(1)": 1 / 0.0219,
    # Fluid Parameters (Light Gas)

    # Body Forces
 
    # Water Patch
    "patch_icpp(1)%geometry": 7,
    "patch_icpp(1)%hcid": 207,
    "patch_icpp(1)%x_centroid": L/2,
    "patch_icpp(1)%y_centroid": L / 2,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1e5,
    "patch_icpp(1)%alpha_rho(1)": 1,
 #   "patch_icpp(1)%alpha_rho(2)": eps * 1,
"patch_icpp(1)%alpha(1)": 1 ,
  #  "patch_icpp(1)%alpha(2)": eps,
     'cantera_file'                 : ctfile,
}

print(json.dumps(data))
