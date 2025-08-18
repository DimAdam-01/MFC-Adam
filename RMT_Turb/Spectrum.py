# Clean von Kármán spectrum (NO exponential cutoff), calibrated numerically to match BOTH:
#   ∫ E(k) dk = 3/2 u'^2   and   L = (∫ E/k dk)/(∫ E dk).
#
# We use the shape F(κ) = κ^4 / (1+κ^2)^(17/6), set
#   E(k) = (a * u'^2)/k0 * F(k/k0)
# and solve constants from closed-form κ-integrals:
#   I0 = ∫ F(κ) dκ,   I1 = ∫ F(κ)/κ dκ
#   a  = (3/2) / I0,
#   k0 = (I1/I0) / L_target.
# This guarantees (in the continuum) the two constraints exactly.
#
# Then we plot and verify on a finite log grid.
import numpy as np
import matplotlib.pyplot as plt

# ---------------------- User parameters ----------------------
U_jet = 100.0            # mean/jet speed [m/s]
intensity = 0.05         # u'/U
u_rms = intensity * U_jet
H = 0.98/1000.0          # characteristic height [m]
L_target = H/3.0         # target integral length [m]

npts_kappa = 400000      # accuracy for κ integrals (shape-only)
npts_k = 200000          # points for plotting/integration in k
kspan = (1e-5, 1e4)      # k-range relative to k0 for verification/plot
# -------------------------------------------------------------

def F_kappa(kappa):
    return kappa**4.0 / (1.0 + kappa**2.0)**(17.0/6.0)

# Compute I0 = ∫ F(κ) dκ and I1 = ∫ F(κ)/κ dκ over a wide log range
tmin, tmax = -20.0, 20.0                      # κ from e^-20 to e^20
t = np.linspace(tmin, tmax, npts_kappa)
kappa = np.exp(t)
Fk = F_kappa(kappa)
I0 = np.trapz(Fk * kappa, t)                  # ∫ F(κ) dκ = ∫ F(e^t) e^t dt
I1 = np.trapz(Fk, t)                          # ∫ F(κ)/κ dκ = ∫ F(e^t) dt

a = (1.5) / I0
k0 = (I1 / I0) / L_target

# Build k-grid and spectrum
kmin = kspan[0] * k0
kmax = kspan[1] * k0
k = np.logspace(np.log10(kmin), np.log10(kmax), npts_k)
E = (a * u_rms**2 / k0) * F_kappa(k / k0)

# Checks on discrete grid
E_int = np.trapz(E, k)
L_spec = np.trapz(E / k, k) / E_int
k_int_target = 2.0*np.pi / L_target
k_int_spec = 2.0*np.pi / L_spec

print("=== Clean von Kármán spectrum (no exponential) ===")
print(f"u_rms                    = {u_rms:.6e} [m/s]  (intensity {intensity:.2%})")
print(f"L_target                 = {L_target:.6e} [m]")
print(f"Solved constants: a      = {a:.10f},  k0 = {k0:.10e} [1/m]")
print("----- Discrete checks (over the plotted k-range) -----")
print(f"∫ E dk                   = {E_int:.8e}  vs 3/2 u'^2 = {(1.5*u_rms**2):.8e}   "
      f"(rel err {(E_int-1.5*u_rms**2)/(1.5*u_rms**2):.3e})")
print(f"L from spectrum          = {L_spec:.8e} [m] vs target {L_target:.8e} [m]   "
      f"(rel err {(L_spec-L_target)/L_target:.3e})")
print(f"k_int target             = {k_int_target:.6e} [1/m];  from spectrum {k_int_spec:.6e} [1/m]")

# Plot
plt.figure(figsize=(10,6))
plt.loglog(k, E, label="E(k) (von Kármán)")
plt.axvline(k0, linestyle=":", label="k0 (solved)")
plt.axvline(k_int_target, linestyle="--", label="k_int = 2π/L_target")
plt.xlabel("k [1/m]"); plt.ylabel("E(k) [m^3/s^2]")
plt.title("Von Kármán spectrum (normalized to u' and L)")
plt.legend(); plt.tight_layout(); plt.show()
