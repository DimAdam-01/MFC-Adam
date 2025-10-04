#!/usr/bin/env python3
import math
import json

# Numerical setup
Nx = 4095
dx = 2.0 / (1.0 * (Nx + 1))
dt = dx/0.5/700
Nt = float(2.0/dt)
# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0,
            "x_domain%end": 2.0,
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": int(Nt),
            "t_step_save": int(math.ceil(Nt / 10.0)),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1 L (-5 < x < -4.5)
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 1.0,
            "patch_icpp(1)%length_x": 2.0,
            "patch_icpp(1)%vel(1)": 0.5,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)":"0.09253158446106331*(1 + 0.1*sin(x*pi))",
            "patch_icpp(1)%alpha(1)": 1.0,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
        }
    )
)
