#!/usr/bin/env python3 import json import argparse import math
import cantera as ct

#!/usr/bin/env python3
import json
import argparse
import math

parser = argparse.ArgumentParser(
    prog="nD_inert_shocktube",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--mfc", type=json.loads, default='{}', metavar="DICT",
                    help="MFC's toolchain's internal state.")
parser.add_argument("--no-chem", dest='chemistry', default=True, action="store_false",
                    help="Disable chemstry.")

args = parser.parse_args()

ctfile    = 'syngas.yaml'
sol_L     = ct.Solution(ctfile)
sol_L.TPX =  300,  101325, 'H:1'
#!/usr/bin/env python3

# SURROUNDING FLOW
# Nondimensional parameters
Re0 = 320.0  # Reynolds number
M0 = 0.2  # Mach number

# Fluid properties
gamma = 1.4

# Free stream velocity & pressure
u0 = 100
pres0 = 1.0 / (gamma * M0**2)

# Domain size
Lx = 0.01152
Ly = 0.01344
Lz = 0.00768

# Number of grid cells
Nx = 127
Ny = 255
Nz = 63

# Grid spacing
dx = Lx / float(Nx)
dy = Ly / float(Ny)
dz = Lz / float(Nz)

# Time advancement
cfl = 0.4
T = 384.0*10**(-6)
dt = dx/cfl/1400.0
Ntfinal = int(T / dt)
Ntstart = 0
Nfiles = 40
t_save = int(math.ceil((Ntfinal - 0) / float(Nfiles)))
Nt = t_save * Nfiles
t_step_start = Ntstart
t_step_stop = int(Nt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": Lx,
            "y_domain%beg": -Ly / 2.0,
            "y_domain%end": Ly / 2.0,
            "z_domain%beg": 0.0,
            "z_domain%end": Lz,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "dt": dt,
            "t_step_start": t_step_start,
            "t_step_stop": t_step_stop,
            "t_step_save": t_save,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-40,
            "weno_Re_flux": "F",
            "wenoz": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "viscous": "T",
    'chemistry'                    : 'T' if not args.chemistry else 'T',
    'chem_params%diffusion'        : 'T',
   'chem_params%reactions'        : 'T',

            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "cons_vars_wrt": "T",
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "fd_order": 1,
 "chem_wrt_T"                  : "T",          
  # Patch 1
            "patch_icpp(1)%geometry": 9,
	    "patch_icpp(1)%hcid": 370,
            "patch_icpp(1)%x_centroid": Lx / 2.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": Lz / 2.0,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%length_z": Lz,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": 100,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": pres0,
            # Mixing layer
            # Fluids Physical Parameters
            # Surrounding liquid
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": Re0,
                "cantera_file": ctfile,
        }
    )
)
