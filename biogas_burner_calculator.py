"""
Biogas Burner Design Optimization Calculator
=============================================
Calculates optimal burner dimensions when converting an LPG burner to biogas.
Based on standard burner design equations for combustion engineering.
"""

import math

# ─── Constants ──────────────────────────────────────────────────────────────
SIGMA = 5.67e-8          # Stefan-Boltzmann constant [W/m²·K⁴]
PI = math.pi

# ─── Biogas Properties ──────────────────────────────────────────────────────
BIOGAS_CALORIFIC_VALUE   = 21_000_000   # J/m³  (21 MJ/m³, ~60% CH4)
BIOGAS_DENSITY           = 1.15         # kg/m³ at standard conditions
BIOGAS_SPECIFIC_GRAVITY  = 0.90         # relative to air (SG, dimensionless)
BIOGAS_WOBBE_INDEX       = 22.17        # MJ/m³  (indicative)

# ─── LPG Reference Properties (for comparison) ──────────────────────────────
LPG_CALORIFIC_VALUE      = 93_100_000   # J/m³
LPG_DENSITY              = 2.01         # kg/m³

# ─── Burner Design Defaults ─────────────────────────────────────────────────
DEFAULT_EFFICIENCY       = 0.60         # 60 % thermal efficiency (biogas)
DEFAULT_VELOCITY         = 0.40         # m/s  (adequate port velocity)
DEFAULT_DELTA_P          = 200          # Pa   (pressure drop across orifice)
CD                       = 0.82         # discharge coefficient
EPSILON                  = 0.85         # emissivity of burner surface
T_HOT                    = 900          # K    (burner surface temperature)
T_AMB                    = 300          # K    (ambient temperature)
H_CONV                   = 15           # W/m²·K (convective coefficient)
VC_VO_RATIO              = 3.45         # chamber vol / orifice vol (design ratio)


def heat_transfer(area: float) -> dict:
    """Total heat loss from burner surface (convective + radiative)."""
    q_conv = area * H_CONV * (T_HOT - T_AMB)
    q_rad  = SIGMA * EPSILON * area * (T_HOT**4 - T_AMB**4)
    q_total = q_conv + q_rad
    return {"q_conv_W": q_conv, "q_rad_W": q_rad, "q_total_W": q_total}


def volumetric_flow_rate(heat_input_W: float,
                         calorific_value: float = BIOGAS_CALORIFIC_VALUE,
                         density: float = BIOGAS_DENSITY,
                         efficiency: float = DEFAULT_EFFICIENCY) -> float:
    """V̇ = q / (C · ρ · η)  [m³/s]"""
    return heat_input_W / (calorific_value * density * efficiency)


def orifice_area_from_cd(v_dot: float,
                         delta_p: float = DEFAULT_DELTA_P,
                         rho: float = BIOGAS_DENSITY) -> float:
    """A = V̇ / [CD · sqrt(2·ΔP/ρ)]  [m²]"""
    return v_dot / (CD * math.sqrt(2 * delta_p / rho))


def nozzle_exit_velocity(v_dot: float, d_o: float) -> float:
    """U_o = 4·V̇ / (π·d_o²)  [m/s]"""
    return 4 * v_dot / (PI * d_o**2)


def initial_jet_momentum(u_o: float, rho_o: float, d_o: float) -> float:
    """I = π·U_o²·ρ_o·d_o² / 4  [N]"""
    return PI * u_o**2 * rho_o * d_o**2 / 4


def entrained_flow_ratio(x: float, d_o: float,
                         sg: float = BIOGAS_SPECIFIC_GRAVITY) -> float:
    """V_e / V_o(SG) = 0.32·x / [d_o · sqrt(S)]"""
    return (0.32 * x) / (d_o * math.sqrt(sg))


def mass_entrainment_ratio(x: float, d_o: float,
                           sg: float = BIOGAS_SPECIFIC_GRAVITY) -> float:
    """m_e / m_o = 1 + 0.32·x / [d_o · sqrt(S)]"""
    return 1 + (0.32 * x) / (d_o * math.sqrt(sg))


def centerline_concentration(x: float, d_o: float) -> float:
    """C_o / C_m = 0.22·(x/d_o) − 1.5"""
    return 0.22 * (x / d_o) - 1.5


def velocity_decay(x: float, d_o: float) -> float:
    """V_o / V_m = 0.16·(x/d_o) − 1.5"""
    return 0.16 * (x / d_o) - 1.5


def design_biogas_burner(heat_input_kW: float     = 2.0,
                         num_ports_hint: int       = None,
                         port_velocity: float      = DEFAULT_VELOCITY,
                         delta_p_pa: float         = DEFAULT_DELTA_P,
                         efficiency: float         = DEFAULT_EFFICIENCY,
                         mixing_tube_x_over_do: float = 30) -> dict:
    """
    Master design function.  Returns a full specification dict.

    Parameters
    ----------
    heat_input_kW       : desired thermal output [kW]
    num_ports_hint      : override port count (None = auto-select)
    port_velocity       : port exit velocity [m/s]
    delta_p_pa          : pressure drop across orifice [Pa]
    efficiency          : combustion efficiency (0–1)
    mixing_tube_x_over_do : mixing tube length as multiple of orifice diameter
    """
    q = heat_input_kW * 1000  # W

    # 1. Volumetric flow rate of biogas
    v_dot = volumetric_flow_rate(q, BIOGAS_CALORIFIC_VALUE,
                                 BIOGAS_DENSITY, efficiency)

    # 2. Total port flow area  (A = V̇ / velocity)
    total_area = v_dot / port_velocity

    # 3. Orifice area from discharge-coefficient equation
    a_orifice = orifice_area_from_cd(v_dot, delta_p_pa, BIOGAS_DENSITY)

    # 4. Injector (orifice) diameter
    d_orifice = math.sqrt(4 * a_orifice / PI) * 1000  # mm

    # 5. Port diameter  (from total area)
    #    Choose number of ports: auto rule → keep port dia 3–6 mm
    if num_ports_hint:
        n_ports = num_ports_hint
    else:
        # target 4 mm port dia
        target_d = 0.004  # m
        a_per_port = PI * target_d**2 / 4
        n_ports = max(4, round(total_area / a_per_port))
        # adjust to nearest even number
        n_ports = max(4, n_ports + (n_ports % 2))

    a_per_port = total_area / n_ports
    d_port_m   = math.sqrt(4 * a_per_port / PI)
    d_port_mm  = d_port_m * 1000

    # 6. Mixing tube length  L = (x/do) · do
    d_o_m = d_orifice / 1000
    L_mix  = mixing_tube_x_over_do * d_o_m * 1000  # mm

    # 7. Mixing tube diameter  (Vc/Vo = 3.45)
    #    Treat orifice volume as the throat cross-section × some unit length
    #    Practical: D_mix ≈ sqrt(3.45) * D_orifice  (area ratio interpretation)
    d_mix_mm = math.sqrt(VC_VO_RATIO) * d_orifice

    # 8. Nozzle exit velocity
    u_o = nozzle_exit_velocity(v_dot, d_o_m)

    # 9. Jet momentum
    I = initial_jet_momentum(u_o, BIOGAS_DENSITY, d_o_m)

    # 10. Entrainment at mixing tube exit  (x = L_mix in metres)
    x_m = L_mix / 1000
    Ve_Vo = entrained_flow_ratio(x_m, d_o_m)
    me_mo = mass_entrainment_ratio(x_m, d_o_m)

    # 11. Burner surface area for heat-transfer check (cylinder approx)
    burner_r  = (n_ports * d_port_mm / (2 * PI)) / 1000   # m
    burner_h  = L_mix / 1000                               # m
    A_surf    = 2 * PI * burner_r * burner_h + PI * burner_r**2
    ht        = heat_transfer(A_surf)

    # 12. Primary aeration ratio  (stoichiometric air ~ 5.7 for biogas)
    stoich_air = 5.7   # m³ air / m³ biogas
    primary_aeration = Ve_Vo / stoich_air * 100  # %

    return {
        # ── Inputs ──
        "heat_input_kW"         : heat_input_kW,
        "efficiency"            : efficiency,
        "delta_p_Pa"            : delta_p_pa,
        "port_velocity_ms"      : port_velocity,

        # ── Flow rates ──
        "v_dot_m3_per_s"        : v_dot,
        "v_dot_m3_per_hr"       : v_dot * 3600,

        # ── Orifice / Injector ──
        "orifice_area_m2"       : a_orifice,
        "orifice_diameter_mm"   : d_orifice,

        # ── Ports ──
        "num_ports"             : n_ports,
        "port_diameter_mm"      : d_port_mm,
        "total_port_area_mm2"   : total_area * 1e6,

        # ── Mixing Tube ──
        "mixing_tube_length_mm" : L_mix,
        "mixing_tube_dia_mm"    : d_mix_mm,
        "mixing_tube_L_D_ratio" : L_mix / d_mix_mm,

        # ── Jet / Entrainment ──
        "nozzle_exit_velocity"  : u_o,
        "jet_momentum_N"        : I,
        "Ve_Vo_ratio"           : Ve_Vo,
        "mass_entrainment_ratio": me_mo,
        "primary_aeration_pct"  : primary_aeration,

        # ── Heat Transfer ──
        "burner_surface_area_m2": A_surf,
        "q_conv_W"              : ht["q_conv_W"],
        "q_rad_W"               : ht["q_rad_W"],
        "q_total_loss_W"        : ht["q_total_W"],
    }


def print_report(spec: dict):
    w = 60
    print("=" * w)
    print(" BIOGAS BURNER DESIGN SPECIFICATION ".center(w))
    print("=" * w)

    sections = [
        ("INPUT CONDITIONS", [
            ("Heat Input",              f"{spec['heat_input_kW']:.2f} kW"),
            ("Thermal Efficiency",      f"{spec['efficiency']*100:.1f} %"),
            ("Pressure Drop (orifice)", f"{spec['delta_p_Pa']:.0f} Pa"),
            ("Port Exit Velocity",      f"{spec['port_velocity_ms']:.2f} m/s"),
        ]),
        ("FUEL FLOW", [
            ("Volumetric Flow Rate",    f"{spec['v_dot_m3_per_s']*1e6:.3f} mL/s  "
                                        f"({spec['v_dot_m3_per_hr']:.4f} m³/hr)"),
        ]),
        ("ORIFICE / INJECTOR", [
            ("Orifice Area",            f"{spec['orifice_area_m2']*1e6:.4f} mm²"),
            ("Orifice Diameter",        f"{spec['orifice_diameter_mm']:.3f} mm"),
        ]),
        ("BURNER PORTS", [
            ("Number of Ports",         f"{spec['num_ports']}"),
            ("Port Diameter",           f"{spec['port_diameter_mm']:.2f} mm"),
            ("Total Port Area",         f"{spec['total_port_area_mm2']:.4f} mm²"),
        ]),
        ("MIXING TUBE", [
            ("Length",                  f"{spec['mixing_tube_length_mm']:.2f} mm"),
            ("Diameter",                f"{spec['mixing_tube_dia_mm']:.2f} mm"),
            ("L/D Ratio",               f"{spec['mixing_tube_L_D_ratio']:.1f}"),
        ]),
        ("JET & ENTRAINMENT", [
            ("Nozzle Exit Velocity",    f"{spec['nozzle_exit_velocity']:.3f} m/s"),
            ("Jet Momentum",            f"{spec['jet_momentum_N']:.6f} N"),
            ("Ve/Vo Ratio",             f"{spec['Ve_Vo_ratio']:.3f}"),
            ("Mass Entrainment Ratio",  f"{spec['mass_entrainment_ratio']:.3f}"),
            ("Primary Aeration",        f"{spec['primary_aeration_pct']:.1f} %"),
        ]),
        ("HEAT TRANSFER (surface losses)", [
            ("Burner Surface Area",     f"{spec['burner_surface_area_m2']*1e4:.2f} cm²"),
            ("Convective Loss",         f"{spec['q_conv_W']:.2f} W"),
            ("Radiative Loss",          f"{spec['q_rad_W']:.2f} W"),
            ("Total Surface Loss",      f"{spec['q_total_loss_W']:.2f} W"),
        ]),
    ]

    for title, rows in sections:
        print(f"\n  {title}")
        print("  " + "-" * (w - 4))
        for label, value in rows:
            print(f"  {label:<30} {value}")

    print("\n" + "=" * w)
    print(" DESIGN GUIDELINES ".center(w))
    print("=" * w)
    print("""
  • Port velocity 0.3–0.5 m/s is optimal for biogas flames.
  • Primary aeration 40–60 % of stoichiometric gives stable flame.
  • L/D ratio of mixing tube should be 15–25 for good mixing.
  • Increase number of ports or reduce port dia if velocity is low.
  • Biogas Wobbe Index ≈ 22 MJ/m³ (vs LPG ≈ 72 MJ/m³); larger
    orifice is required (approx 3× area compared to LPG).
""")
    print("=" * w)


# ─── CLI Entry Point ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("\nBIOGAS BURNER DESIGN CALCULATOR")
    print("-" * 40)
    print("Press Enter to accept default values.\n")

    try:
        q  = float(input("  Heat input [kW]         (default 2.0) : ") or 2.0)
        np = input(      "  Number of ports         (default auto): ").strip()
        pv = float(input("  Port velocity [m/s]     (default 0.4) : ") or 0.4)
        dp = float(input("  Orifice pressure drop [Pa](default 200): ") or 200)
        ef = float(input("  Combustion efficiency   (default 0.60): ") or 0.60)

        num_ports = int(np) if np else None
        spec = design_biogas_burner(
            heat_input_kW        = q,
            num_ports_hint       = num_ports,
            port_velocity        = pv,
            delta_p_pa           = dp,
            efficiency           = ef,
        )
        print_report(spec)

    except ValueError as e:
        print(f"\nInput error: {e}")
    except Exception as e:
        print(f"\nCalculation error: {e}")
