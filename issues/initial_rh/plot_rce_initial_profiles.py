"""
plot_rce_initial_profiles.py
============================
Compute and plot the initial atmospheric profiles for the LFRic RCE idealized
case using the relative humidity initialization introduced in lfric_apps
(initial_profiles branch).

The script reproduces the analytic initialization logic from:
  - init_thermo_profile_alg_mod  (theta interpolation)
  - init_rh_profile_alg_mod      (RH → q_v conversion)
  - physics_common_mod::qsaturation (Tetens saturation formula)
  - gungho_init_prognostics_driver_mod  (hydrostatic balance)

All values come directly from:
  rose-stem/app/lfric_atm/opt/rose-app-rce.conf
  rose-stem/app/lfric_atm/rose-app.conf (extrusion, planet)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')   # non-interactive — works in WSL without a display
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d

# ---------------------------------------------------------------------------
# Physical constants  (from namelist:planet)
# ---------------------------------------------------------------------------
Rd     = 287.05          # J kg-1 K-1   dry air gas constant
Rv     = 461.50          # J kg-1 K-1   water vapour gas constant
cp     = 1005.0          # J kg-1 K-1   specific heat of dry air at const. p
g      = 9.80665         # m s-2        gravitational acceleration
p_zero = 100000.0        # Pa           reference pressure
kappa  = Rd / cp         # Rd/cp
epsilon = Rd / Rv        # ~0.622
recip_epsilon = Rv / Rd  # ~1.608  (= 1/epsilon in model notation)

# ---------------------------------------------------------------------------
# Vertical grid — UM L70_50t_20s_80km (from gungho_extrusion_mod.f90)
# ---------------------------------------------------------------------------
domain_height = 80000.0   # m
eta = np.array([
    0.0000000, 0.0002500, 0.0006667, 0.0012500, 0.0020000, 0.0029167,
    0.0040000, 0.0052500, 0.0066667, 0.0082500, 0.0100000, 0.0119167,
    0.0140000, 0.0162500, 0.0186667, 0.0212500, 0.0240000, 0.0269167,
    0.0300000, 0.0332500, 0.0366667, 0.0402500, 0.0440000, 0.0479167,
    0.0520000, 0.0562500, 0.0606667, 0.0652500, 0.0700000, 0.0749167,
    0.0800000, 0.0852500, 0.0906668, 0.0962505, 0.1020017, 0.1079213,
    0.1140113, 0.1202745, 0.1267154, 0.1333406, 0.1401592, 0.1471838,
    0.1544313, 0.1619238, 0.1696895, 0.1777643, 0.1861929, 0.1950307,
    0.2043451, 0.2142178, 0.2247466, 0.2360480, 0.2482597, 0.2615432,
    0.2760868, 0.2921094, 0.3098631, 0.3296378, 0.3517651, 0.3766222,
    0.4046373, 0.4362943, 0.4721379, 0.5127798, 0.5589045, 0.6112759,
    0.6707432, 0.7382500, 0.8148403, 0.9016668, 1.0000000
])
# Wtheta points: interfaces (71 levels including surface and top)
z_wth = eta * domain_height                 # height at Wtheta (theta/q) levels
# W3 (rho) points: midpoints between Wtheta levels
z_w3  = 0.5 * (z_wth[:-1] + z_wth[1:])    # 70 rho levels

# ---------------------------------------------------------------------------
# Theta profile  (from namelist:initial_temperature in rose-app-rce.conf)
# ---------------------------------------------------------------------------
theta_heights = np.array([
    0.0, 800.0, 1200.0, 3500.0, 4100.0, 8200.0, 12500.0, 13500.0,
    14200.0, 16.0e3, 20.0e3, 24.0e3, 28.0e3, 32.0e3, 36.0e3, 40.0e3
])
theta_data = np.array([
    297.0, 297.0, 300.0, 306.5, 311.0, 318.0, 328.5, 333.0,
    340.0, 371.0, 483.0, 610.0, 738.0, 928.0, 1227.0, 1447.0
])

# ---------------------------------------------------------------------------
# Relative humidity profile  (from namelist:initial_vapour in rose-app-rce.conf)
# RH is a fraction: 0 = completely dry, 1 = saturated
# ---------------------------------------------------------------------------
rh_heights = np.array([
    0.0, 680.0, 1300.0, 3500.0, 4150.0, 4850.0, 5200.0, 6100.0, 7.0e3,
    8150.0, 9.5e3, 10.5e3, 11.5e3, 12.25e3, 13.0e3, 14.0e3, 18.0e3, 40.0e3
])
rh_data = np.array([
    0.70, 0.70, 0.65, 0.55, 0.52, 0.48, 0.46, 0.40, 0.34, 0.26,
    0.17, 0.12, 0.08, 0.06, 0.05, 0.05, 0.05, 0.05
])

surface_pressure = 101325.0   # Pa (base conf default)

# ---------------------------------------------------------------------------
# Tetens saturation mixing ratio  (physics_common_mod::qsaturation)
# Arguments: T in K, p in mbar.  Returns q_sat in kg/kg.
# ---------------------------------------------------------------------------
QSA1 =  3.8           # kg/kg * mbar
QSA2 = -17.2693882    # K-1
QSA3 =  35.86         # K
QSA4 =  6.109         # mbar
TK0C = 273.15         # K

def qsaturation(T, p_mbar):
    """Tetens saturation mixing ratio (kg/kg). T in K, p in mbar."""
    denom = p_mbar * np.exp(QSA2 * (T - TK0C) / (T - QSA3)) - QSA4
    # Guard against unphysical values (T near or below QSA3)
    with np.errstate(invalid='ignore', divide='ignore'):
        qs = np.where((T > QSA3) & (denom > 0), QSA1 / denom, 999.0)
    return qs

# ---------------------------------------------------------------------------
# Interpolate profiles onto Wtheta model levels (linear, constant extrapolation)
# Matches profile_interp_kernel_mod behaviour: extrapolate at boundaries.
# ---------------------------------------------------------------------------
def interp_profile(z_out, z_in, data_in):
    """Linear interpolation with constant extrapolation at boundaries."""
    f = interp1d(z_in, data_in, kind='linear', bounds_error=False,
                 fill_value=(data_in[0], data_in[-1]))
    return f(z_out)

theta_wth = interp_profile(z_wth, theta_heights, theta_data)  # K (pot. temp.)
rh_wth    = interp_profile(z_wth, rh_heights,    rh_data)     # fraction

# ---------------------------------------------------------------------------
# Hydrostatic pressure profile
# Integrate dp/dz = -rho * g = -p/(Rd*T) * g using the theta profile.
# This is a dry first-guess (moisture contribution to density is small and
# would require iteration).
# ---------------------------------------------------------------------------
# Start with T profile: T = theta * Exner, Exner = (p/p_zero)^kappa
# Integrate hydrostatically from surface.

nz = len(z_wth)
p_wth    = np.zeros(nz)
exner_wth = np.zeros(nz)

p_wth[0] = surface_pressure
exner_wth[0] = (p_wth[0] / p_zero) ** kappa
T_wth    = np.zeros(nz)
T_wth[0] = theta_wth[0] * exner_wth[0]

T_MIN = 1.0   # K  — floor to avoid divide-by-zero in deep stratosphere
for k in range(1, nz):
    dz   = z_wth[k] - z_wth[k-1]
    T_km = max(T_wth[k-1], T_MIN)
    # Isothermal layer first guess, then midpoint correction
    p_mid  = p_wth[k-1] * np.exp(-g * dz / (Rd * T_km))
    ex_mid = (p_mid / p_zero) ** kappa
    T_mid  = max(0.5 * (theta_wth[k-1] + theta_wth[k]) * ex_mid, T_MIN)
    p_wth[k] = p_wth[k-1] * np.exp(-g * dz / (Rd * T_mid))
    exner_wth[k] = (p_wth[k] / p_zero) ** kappa
    T_wth[k] = theta_wth[k] * exner_wth[k]

# ---------------------------------------------------------------------------
# Saturation and actual mixing ratio
# ---------------------------------------------------------------------------
p_mbar   = p_wth / 100.0                         # Pa -> mbar
qsat_wth = qsaturation(T_wth, p_mbar)            # kg/kg

# Convert RH fraction -> mixing ratio (rh_to_mr_kernel_mod formula)
# mr_v = RH * q_sat / (1 + (1 - RH) * q_sat * recip_epsilon)
mr_v_wth = rh_wth * qsat_wth / (1.0 + (1.0 - rh_wth) * qsat_wth * recip_epsilon)

# Re-derive RH from mr_v (verify round-trip, and for plotting)
# rel_hum = (mr_v / q_sat) * (1 + q_sat * recip_eps) / (1 + mr_v * recip_eps)
rh_check = (mr_v_wth / qsat_wth) * (1.0 + qsat_wth * recip_epsilon) / \
           (1.0 + mr_v_wth * recip_epsilon)

# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------
# Static stability: dtheta/dz  (K/km for plotting)
dtheta_dz = np.gradient(theta_wth, z_wth) * 1000.0   # K/km

# Brunt-Väisälä frequency squared (dry): N² = (g/theta) * d(theta)/dz
N2_dry = (g / theta_wth) * np.gradient(theta_wth, z_wth)   # s-2

# Temperature lapse rate: -dT/dz  (K/km), positive = unstable
lapse = -np.gradient(T_wth, z_wth) * 1000.0   # K/km

# ---------------------------------------------------------------------------
# Restrict plots to troposphere + lower stratosphere for readability
# ---------------------------------------------------------------------------
z_top_km = 25.0   # km
km = 1e3
mask = z_wth <= z_top_km * km   # mask for Wtheta levels

z_km = z_wth[mask] / km

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
plt.rcParams.update({
    'font.family':    'DejaVu Sans',
    'font.size':      11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.minor.width': 0.8,
    'ytick.minor.width': 0.8,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'figure.dpi': 150,
})

PANEL_COLOR = '#f7f7f7'

fig = plt.figure(figsize=(16, 10))
fig.patch.set_facecolor('white')

gs = gridspec.GridSpec(2, 4, figure=fig,
                       hspace=0.38, wspace=0.40,
                       left=0.07, right=0.97,
                       top=0.91, bottom=0.08)

ax1 = fig.add_subplot(gs[:, 0])   # tall — temperature / theta
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 1])
ax4 = fig.add_subplot(gs[0, 2])
ax5 = fig.add_subplot(gs[1, 2])
ax6 = fig.add_subplot(gs[0, 3])
ax7 = fig.add_subplot(gs[1, 3])

axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]
for ax in axes:
    ax.set_facecolor(PANEL_COLOR)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.6, color='grey')
    ax.set_ylabel('Height (km)')
    ax.set_ylim(0, z_top_km)
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))

# ------------------------------------------------------------------
# Panel 1: Theta and Temperature (tall, full height for context)
# ------------------------------------------------------------------
l1, = ax1.plot(theta_wth[mask], z_km, color='#c0392b', lw=2.0, label='θ (K)')
ax1b = ax1.twiny()
ax1b.set_facecolor(PANEL_COLOR)
ax1b.spines[['top']].set_visible(True)
l2, = ax1b.plot(T_wth[mask], z_km, color='#2980b9', lw=2.0,
                linestyle='--', label='T (K)')
ax1.set_xlabel('Potential Temperature θ (K)', color='#c0392b')
ax1b.set_xlabel('Temperature T (K)', color='#2980b9')
ax1.tick_params(axis='x', colors='#c0392b')
ax1b.tick_params(axis='x', colors='#2980b9')
ax1.set_title('Temperature')
# MALR reference
z_ref = np.linspace(0, z_top_km * km, 300)
MALR = 9.8   # K/km (dry adiabatic lapse rate)
T_DALR = T_wth[0] - MALR * z_ref / km
ax1b.plot(T_DALR, z_ref / km, color='#7f8c8d', lw=1.0,
          linestyle=':', label='DALR')
lines = [l1, l2, plt.Line2D([0], [0], color='#7f8c8d', lw=1.0,
                              linestyle=':', label='DALR')]
ax1.legend(lines, [l.get_label() for l in lines],
           loc='upper left', framealpha=0.8)

# ------------------------------------------------------------------
# Panel 2: Input RH profile  (namelist values — what user specifies)
# ------------------------------------------------------------------
ax2.step(rh_data * 100, rh_heights / km, where='post',
         color='#8e44ad', lw=1.5, linestyle=':', label='namelist nodes')
ax2.plot(rh_wth[mask] * 100, z_km, color='#8e44ad', lw=2.0,
         label='interpolated')
ax2.set_xlabel('Relative Humidity (%)')
ax2.set_title('Specified RH Profile')
ax2.xaxis.set_minor_locator(plt.MultipleLocator(5))
ax2.set_xlim(0, 100)
ax2.legend(loc='upper right', framealpha=0.8)

# ------------------------------------------------------------------
# Panel 3: Water vapour mixing ratio
# ------------------------------------------------------------------
ax3.plot(mr_v_wth[mask] * 1e3, z_km, color='#16a085', lw=2.0)
ax3.set_xlabel(r'Water vapour mixing ratio $q_v$ (g kg$^{-1}$)')
ax3.set_title('Water Vapour Mixing Ratio')
ax3.xaxis.set_minor_locator(plt.MultipleLocator(0.5))

# ------------------------------------------------------------------
# Panel 4: Saturation mixing ratio
# ------------------------------------------------------------------
ax4.plot(qsat_wth[mask] * 1e3, z_km, color='#2980b9', lw=2.0,
         label=r'$q_{sat}$')
ax4.plot(mr_v_wth[mask] * 1e3, z_km, color='#16a085', lw=2.0,
         linestyle='--', label=r'$q_v$')
ax4.set_xlabel(r'Mixing ratio (g kg$^{-1}$)')
ax4.set_title(r'$q_{sat}$ and $q_v$')
ax4.legend(loc='upper right', framealpha=0.8)

# ------------------------------------------------------------------
# Panel 5: Static stability  dθ/dz
# ------------------------------------------------------------------
ax5.axvline(0, color='#7f8c8d', lw=0.8, linestyle=':')
ax5.plot(dtheta_dz[mask], z_km, color='#e67e22', lw=2.0)
ax5.set_xlabel(r'$\partial\theta/\partial z$ (K km$^{-1}$)')
ax5.set_title('Static Stability')

# ------------------------------------------------------------------
# Panel 6: N² (dry Brunt-Väisälä)
# ------------------------------------------------------------------
ax6.axvline(0, color='#7f8c8d', lw=0.8, linestyle=':')
ax6.plot(N2_dry[mask] * 1e4, z_km, color='#e74c3c', lw=2.0)
ax6.set_xlabel(r'$N^2_\mathrm{dry}$ ($\times 10^{-4}$ s$^{-2}$)')
ax6.set_title(r'Dry Brunt–Väisälä Freq. Squared')
ax6.xaxis.set_minor_locator(plt.MultipleLocator(0.2))

# ------------------------------------------------------------------
# Panel 7: Temperature lapse rate with reference lines
# ------------------------------------------------------------------
ax7.axvline(9.8, color='#7f8c8d', lw=1.0, linestyle=':', label='DALR (9.8 K/km)')
ax7.axvline(6.5, color='#bdc3c7', lw=1.0, linestyle='--', label='ISA (6.5 K/km)')
ax7.plot(lapse[mask], z_km, color='#2c3e50', lw=2.0, label='RCE profile')
ax7.set_xlabel(r'Lapse rate $-\partial T/\partial z$ (K km$^{-1}$)')
ax7.set_title('Temperature Lapse Rate')
ax7.legend(loc='upper right', framealpha=0.8, fontsize=9)

# ------------------------------------------------------------------
# Global title
# ------------------------------------------------------------------
fig.suptitle(
    'LFRic RCE Idealized Case — Initial Atmospheric Profiles\n'
    r'(UM L70 vertical grid, $p_s = 1013.25$ hPa, '
    r'relative humidity initialization)',
    fontsize=13, fontweight='bold', y=0.98
)

# ------------------------------------------------------------------
# Verify round-trip: max |RH_input - RH_back-computed|
# ------------------------------------------------------------------
max_rh_err = np.max(np.abs(rh_wth[mask] - rh_check[mask]))
print(f"Max round-trip RH error (input vs back-computed): {max_rh_err:.2e}")
print(f"Surface q_sat:   {qsat_wth[0]*1e3:.2f} g/kg")
print(f"Surface q_v:     {mr_v_wth[0]*1e3:.2f} g/kg  (RH={rh_wth[0]*100:.0f}%)")
print(f"Surface T:       {T_wth[0]:.1f} K")
print(f"Surface theta:   {theta_wth[0]:.1f} K")
print(f"Surface pressure:{p_wth[0]/100:.1f} hPa")

outfile = '/mnt/c/Users/bship/code_development/lfric_idealized/rce_initial_profiles.png'
fig.savefig(outfile, dpi=200, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {outfile}")
plt.close(fig)
