#!/usr/bin/env python3
import json
import math
Lx = 0.04347
Ly = 0.01466+3.048/1000
Cav_Lx = 43.47/1000-21.99/1000-14.19/1000
Cav_Ly = 3.048/1000
Cav_CenterX = (Lx - 21.99/1000-14.15/1000)/2
Cav_CenterY = 3.048/2000 
Cav_Center2X = 2*Cav_CenterX+14.15/1000+21.99/2000

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": 0.04347,
            "y_domain%beg": 0.0,
            "y_domain%end": Ly,
            "m": 499,
            "n": 499,
            "p": 0,
            "dt": 5e-09,
            "t_step_start": 0,
            "t_step_stop":60,
            "t_step_save": 20,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "mp_weno": "F",
            # 'recon_type'                   : 1,
            "weno_order": 5,
            "weno_eps": 1e-16,
            #'muscl_order'                  : 2,
            #'muscl_lim'                    : 1,
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -7,
            "bc_x%end": -8,
            "bc_y%beg": -16,
            "bc_y%end": -16,
              "bc_x%grcbc_in": "T",
              "bc_x%vel_in(1)": 50.0,
                 "bc_x%vel_in(2)": 0,
             "bc_x%pres_in": 101325.0,
             "bc_x%alpha_rho_in(1)": 1.00,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
                 "ib": "T",
            "num_ibs": 2,
            # Patch 1: Base
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%hcid": 291,
            "patch_icpp(1)%x_centroid": Lx/2,
            "patch_icpp(1)%y_centroid": Ly/2,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%vel(1)": 0,
            "patch_icpp(1)%vel(2)": 0,
            "patch_icpp(1)%pres": 101325,
            "patch_icpp(1)%alpha_rho(1)": 1.00,
            "patch_icpp(1)%alpha(1)": 1,
            # Patch 1: IBM
             "patch_ib(1)%geometry": 3,
            "patch_ib(1)%x_centroid": 0,
            "patch_ib(1)%y_centroid": 0,
            "patch_ib(1)%length_x": 2*Cav_Lx,
            "patch_ib(1)%length_y": 2*Cav_Ly,
            "patch_ib(1)%slip": "F",
            # Patch 2: IBM
            "patch_ib(2)%geometry": 3,
            "patch_ib(2)%x_centroid": Lx,
            "patch_ib(2)%y_centroid": 0,
            "patch_ib(2)%length_x": 2*21.99/1000,
            "patch_ib(2)%length_y": 2*Cav_Ly,
            "patch_ib(2)%slip": "F",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0e00,
             "viscous": "T",
            "fluid_pp(1)%Re(1)": 100000,
        }
    )
)
