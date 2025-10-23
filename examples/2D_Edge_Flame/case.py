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

ctfile    = 'h2o2.yaml'
sol_L     = ct.Solution(ctfile)
sol_L.TPX =  300,  101325, 'H:1'

Lx    = 3.6/1000
Ly    = 7.2/1000
Nx   = 500
Ny   = 1000
dx   = Lx / Nx
dy   = Ly / Ny
dt   = 3e-8
Tend = 0.9e-4

NT         = int(Tend / dt)
SAVE_COUNT = 1
NS         = 4


# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": 0,
    "x_domain%end":Lx,
    "y_domain%beg": 0,
    "y_domain%end": Ly,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "cyl_coord": "F",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": 10000,
    "t_step_save": 300,
    "t_step_print": 1,
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
    "bc_x%beg": -11,
    "bc_x%end": -12,
    "bc_y%beg": -11,
    "bc_y%end": -2,
    "num_patches": 1,
    "num_fluids": 1,
    "viscous": "F",
    'chemistry'                    : 'T' if not args.chemistry else 'T',
    'chem_params%diffusion'        : 'F',
    'chem_params%reactions'        : 'T',
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
 "chem_wrt_T"                  : "T",
    # Fluid Parameters (Heavy Gas)
    "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0.0e00,
    #"fluid_pp(1)%Re(1)": 1 / 0.0219,
    # Fluid Parameters (Light Gas)

    # Body Forces
 
    # Water Patch
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%hcid": 270,
    "patch_icpp(1)%x_centroid": Lx/2,
    "patch_icpp(1)%y_centroid": Ly/2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
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
