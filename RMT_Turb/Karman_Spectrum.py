#!/usr/bin/env python3
"""
Calibrate an isotropic von Kármán 1D energy spectrum E(k) so that:
  ∫ E(k) dk = (3/2) u_rms^2
  L = (π / (2 u_rms^2)) ∫ E(k)/k dk

Then plot E(k) and mark k = 1/L (vertical) and E(k=1/L) (horizontal).
"""

import math
import argparse
import numpy as np
import matplotlib.pyplot as plt


def E_vk(k, L0, C):
    """Von Kármán 1D spectrum: E(k) = C (k L0)^4 / (1+(k L0)^2)^(17/6)."""
    kl = k * L0
    return C * (kl**4) / ((1.0 + kl**2) ** (17.0 / 6.0))


def calibrate_von_karman(u_rms, L_target, kmin, kmax, nk, iters=40, rtol=1e-8):
    """
    Solve for (L0, C) such that:
      ∫ E dk = 1.5 u_rms^2
      L = (π / (2 u_rms^2)) ∫ E/k dk
    using bisection on L0 and analytic scaling for C.

    Returns: k (array), Ek (array), dk (float), L0 (float), C (float)
    """
    if not (kmax > kmin > 0.0):
        raise ValueError("Require 0 < kmin < kmax.")

    # k-grid (avoid k=0)
    k = np.linspace(kmin, kmax, nk, dtype=float)
    dk = (kmax - kmin) / (nk - 1)

    # Helper: implied integral length for a given L0, assuming C normalized by energy
    # With C=1, define:
    #   I1(L0) = ∫ E(k;L0,1) dk
    #   I2(L0) = ∫ E(k;L0,1)/k dk
    # Energy constraint => C = (1.5 u'^2)/I1
    # Implied L = (π/(2 u'^2)) * C * I2 = (3π/4) * I2/I1  (independent of u_rms magnitude)
    def implied_L(L0):
        E0 = E_vk(k, L0, 1.0)
        I1 = np.trapz(E0, k)
        # Guard near k=0 (we already avoid exact 0, but be safe)
        I2 = np.trapz(E0 / np.maximum(k, 1e-30), k)
        return (3.0 * math.pi / 4.0) * (I2 / max(I1, 1e-30))

    # Bisection bracket for L0 around the desired L (works robustly for VK)
    L0_lo = max(1e-12, 0.3 * L_target)
    L0_hi = 3.0 * L_target

    L_lo = implied_L(L0_lo)
    L_hi = implied_L(L0_hi)

    # If bracket is inverted (rare), expand it
    expand_count = 0
    while (L_lo - L_target) * (L_hi - L_target) > 0 and expand_count < 10:
        # Move ends outward geometrically
        L0_lo *= 0.5
        L0_hi *= 2.0
        L_lo = implied_L(L0_lo)
        L_hi = implied_L(L0_hi)
        expand_count += 1

    # Bisection
    for _ in range(iters):
        L0_mid = 0.5 * (L0_lo + L0_hi)
        L_mid = implied_L(L0_mid)
        if abs(L_mid - L_target) <= rtol * max(L_target, 1.0):
            L0 = L0_mid
            break
        if (L_lo - L_target) * (L_mid - L_target) <= 0:
            L0_hi, L_hi = L0_mid, L_mid
        else:
            L0_lo, L_lo = L0_mid, L_mid
    else:
        L0 = 0.5 * (L0_lo + L0_hi)

    # Final normalization C from energy constraint
    E0 = E_vk(k, L0, 1.0)
    I1 = np.trapz(E0, k)                     # ∫ E0 dk
    C = (1.5 * u_rms * u_rms) / max(I1, 1e-30)
    Ek = C * E0

    return k, Ek, dk, L0, C


def main():
    ap = argparse.ArgumentParser(description="Calibrate and plot von Kármán E(k).")
    ap.add_argument("--u_rms", type=float, default=0.05,
                    help="Per-component RMS fluctuation u' (default: 0.05).")
    ap.add_argument("--L", type=float, default=0.10,
                    help="Integral length scale L (default: 0.10).")
    ap.add_argument("--kmin", type=float, default=None,
                    help="Minimum k (default: 2π/Lbox).")
    ap.add_argument("--kmax", type=float, default=None,
                    help="Maximum k (default: 0.6π/dx).")
    ap.add_argument("--Lbox", type=float, default=1.0,
                    help="Box length for default kmin (default: 1.0).")
    ap.add_argument("--dx", type=float, default=0.01,
                    help="Grid spacing for default kmax (default: 0.01).")
    ap.add_argument("--nk", type=int, default=256,
                    help="Number of k-samples (default: 256).")
    ap.add_argument("--out", type=str, default="vk_spectrum.dat",
                    help="Output file for k and E(k) (default: vk_spectrum.dat).")
    args = ap.parse_args()

    # Defaults for spectral band if not provided
    kmin = args.kmin if args.kmin is not None else 2.0 * math.pi / args.Lbox
    kmax = args.kmax if args.kmax is not None else 0.6 * math.pi / args.dx

    k, Ek, dk, L0, C = calibrate_von_karman(args.u_rms, args.L, kmin, kmax, args.nk)

    print("# Calibrated von Kármán spectrum")
    print(f"u_rms    = {args.u_rms:.6e}")
    print(f"L        = {args.L:.6e}")
    print(f"L0       = {L0:.6e}")
    print(f"C        = {C:.6e}")
    print(f"kmin     = {k[0]:.6e}")
    print(f"kmax     = {k[-1]:.6e}")
    print(f"dk       = {dk:.6e}")
    print(f"∫E dk    = {np.trapz(Ek, k):.6e}  (target {(1.5*args.u_rms*args.u_rms):.6e})")

    # Save spectrum
    np.savetxt(args.out, np.column_stack([k, Ek]), header="k   E(k)  (von Kármán calibrated)")

    # Plot E(k) and mark k = 1/L with vertical and horizontal lines
    k_ref = 1.0 / args.L
    # Find E at k_ref by interpolation (within bounds)
    if k[0] <= k_ref <= k[-1]:
        Ek_ref = np.interp(k_ref, k, Ek)
    else:
        Ek_ref = None

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(k, Ek, lw=2, label="von Kármán E(k) (calibrated)")
    if k_ref is not None and k[0] < k_ref < k[-1]:
        ax.axvline(k_ref, ls="--", lw=1.2, label="k = 1/L")
        ax.axhline(Ek_ref, ls="--", lw=1.0)
        ax.plot([k_ref], [Ek_ref], marker="o", ms=5)
        ax.annotate("kL = 1", xy=(k_ref, Ek_ref), xytext=(1.2*k_ref, 1.5*Ek_ref),
                    arrowprops=dict(arrowstyle="->", lw=1.0))

    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$E(k)$")
    ax.set_title("Calibrated von Kármán Spectrum")
    ax.grid(True, which="both", ls=":")
    ax.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

