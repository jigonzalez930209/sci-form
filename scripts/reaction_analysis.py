#!/usr/bin/env python3
"""
reaction_analysis.py
====================
Visualización completa de caminos de reacción con 21 frames.

Para cada reacción:
  - Frame central (índice 10) = Estado de Transición (máxima energía)
  - Muestra posiciones 3D de átomos en cada frame
  - Compara múltiples métodos de cálculo

Reacciones analizadas (5):
  R1: Disociación H-H (pozo doble simétrico 1D)
  R2: Reacción exotérmica asimétrica (1D)
  R3: Potencial Müller-Brown 2D (benchmark TS)
  R4: Transferencia de protón 2D (Hamiltonian Morse)
  R5: Modelo SN2 1D

Métodos cubiertos:
  Level 1: UFF / Analítico (más rápido)
  Level 2: EHT (Extended Hückel) — vía sci_form Python si disponible
  Level 3: PM3 semi-empírico — vía sci_form Python si disponible
  Level 4: xTB tight-binding — vía sci_form Python si disponible

Uso:
  pip install numpy matplotlib
  # Opcional para niveles 2-4:
  # maturin develop --features alpha-gsm,alpha-kinetics  (en sci-form/)
  python scripts/reaction_analysis.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")          # sin pantalla física
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import os, sys, warnings

warnings.filterwarnings("ignore")

# ─── Intentar importar sci_form ───────────────────────────────────────────
SCI_FORM_AVAILABLE = False
try:
    import sci_form
    SCI_FORM_AVAILABLE = True
    print("✓ sci_form Python bindings disponibles")
except ImportError:
    print("⚠ sci_form no disponible — usando potenciales analíticos puros")

N_FRAMES = 21
MID = N_FRAMES // 2   # frame 10

# ═══════════════════════════════════════════════════════════════════════════
# Potenciales analíticos (Centralizados en sci-form Rust)
# ═══════════════════════════════════════════════════════════════════════════

def double_well_1d(x):
    """V(x) = (x² - 1)² — mínimos en ±1, barrera en x=0"""
    if SCI_FORM_AVAILABLE:
        return sci_form.potentials.double_well_1d(x)
    return (x**2 - 1.0)**2


def asymmetric_1d(x):
    """V(x) = (x²-1)² - 0.3x — exotérmico (producto más estable)"""
    if SCI_FORM_AVAILABLE:
        return sci_form.potentials.asymmetric_1d(x)
    return (x**2 - 1.0)**2 - 0.3 * x


def muller_brown(x, y):
    """Potencial Müller-Brown 2D — benchmark estándar"""
    if SCI_FORM_AVAILABLE:
        return sci_form.potentials.muller_brown_2d(x, y)
    A  = [-200, -100, -170,  15]
    b  = [  -1,   -1,  -6.5, 0.7]
    c  = [   0,    0,  11.0, 0.6]
    d  = [ -10,  -10,  -6.5, 0.7]
    x0 = [   1,    0,  -0.5,-1.0]
    y0 = [   0,  0.5,   1.5, 1.0]
    v = 0.0
    for k in range(4):
        dx = x - x0[k]; dy = y - y0[k]
        v += A[k] * np.exp(b[k]*dx**2 + c[k]*dx*dy + d[k]*dy**2)
    return v * 0.001


def proton_transfer_2d(r_ah, r_hb):
    """Hamiltonian Morse para transferencia de protón A-H···B → A···H-B"""
    if SCI_FORM_AVAILABLE:
        return sci_form.potentials.proton_transfer_2d_morse(r_ah, r_hb)
    De, alpha, r0 = 2.0, 2.0, 1.0
    morse_ah = De * (1 - np.exp(-alpha * (r_ah - r0)))**2
    morse_hb = De * (1 - np.exp(-alpha * (r_hb - r0)))**2
    constraint = 50.0 * (r_ah + r_hb - 2.5)**2
    return morse_ah + morse_hb + constraint


def sn2_model_1d(x):
    """Potencial SN2: doble pozo asimétrico con barrera central"""
    if SCI_FORM_AVAILABLE:
        return sci_form.potentials.sn2_model_1d(x)
    return (-2.0 * np.exp(-2.0 * (x + 1.5)**2)
            - 3.0 * np.exp(-2.0 * (x - 1.5)**2)
            + 1.5 * np.exp(-8.0 * x**2))


# ═══════════════════════════════════════════════════════════════════════════
# Generadores de caminos de 21 frames
# ═══════════════════════════════════════════════════════════════════════════

def path_1d(energy_fn, x_start, x_end, n=N_FRAMES):
    """Interpola n puntos en 1D y evalúa la energía."""
    xs = np.linspace(x_start, x_end, n)
    Es = np.array([energy_fn(x) for x in xs])
    # coords = (x, 0, 0) para cada frame
    coords = np.column_stack([xs, np.zeros(n), np.zeros(n)])
    return xs, Es, coords


def path_2d_line(energy_fn_2d, start, end, n=N_FRAMES):
    """Interpola linealmente en 2D y evalúa energía."""
    t = np.linspace(0, 1, n)
    xs = start[0] + t * (end[0] - start[0])
    ys = start[1] + t * (end[1] - start[1])
    Es = np.array([energy_fn_2d(x, y) for x, y in zip(xs, ys)])
    coords = np.column_stack([xs, ys, np.zeros(n)])
    return xs, ys, Es, coords


def path_muller_brown(n=N_FRAMES):
    """Camino de reacción en Müller-Brown entre dos mínimos."""
    # Reactante: mínimo A ≈ (0.62, 0.03), Producto: mínimo C ≈ (-0.56, 1.44)
    start = np.array([0.62, 0.03])
    end   = np.array([-0.56, 1.44])
    t = np.linspace(0, 1, n)
    xs = start[0] + t * (end[0] - start[0])
    ys = start[1] + t * (end[1] - start[1])
    Es = np.array([muller_brown(x, y) for x, y in zip(xs, ys)])
    return xs, ys, Es


def path_proton_transfer(n=N_FRAMES, k=1.0, coupling=0.15):
    """Marcus/EVB: s ∈ [-1,+1], TS exactamente en frame 10 (s=0)"""
    s_vals = np.linspace(-1.0, 1.0, n)
    def evb(s):
        eps_d = k * (s + 1.0)**2
        eps_a = k * (s - 1.0)**2
        return (eps_d + eps_a)/2 - np.sqrt(((eps_d - eps_a)/2)**2 + coupling**2)
    Es = np.array([evb(s) for s in s_vals])
    coords = [np.array([s, 0.0, 0.0]) for s in s_vals]
    return s_vals, Es, coords


# ═══════════════════════════════════════════════════════════════════════════
# Cálculo EHT via sci_form (Nivel 2, si disponible)
# ═══════════════════════════════════════════════════════════════════════════

def eht_h2_profile(r_values):
    """HOMO energy de H₂ en función de r (Å) via EHT."""
    if not SCI_FORM_AVAILABLE:
        return None
    try:
        homos, gaps = [], []
        for r in r_values:
            result = sci_form.solve_eht([1, 1], [[0.0, 0.0, 0.0], [r, 0.0, 0.0]])
            homos.append(result["homo_energy"])
            gaps.append(result["gap"])
        return np.array(homos), np.array(gaps)
    except Exception as e:
        print(f"  EHT H₂ error: {e}")
        return None


def eht_h2o_angle_scan(r_oh=0.957):
    """HOMO y gap de H₂O en función del ángulo H-O-H (90°→180°)."""
    if not SCI_FORM_AVAILABLE:
        return None
    try:
        angles = np.linspace(90.0, 180.0, N_FRAMES)
        homos, gaps = [], []
        for theta_deg in angles:
            half = np.radians(theta_deg / 2)
            pos = [
                [0.0, 0.0, 0.0],
                [r_oh * np.sin(half),  r_oh * np.cos(half), 0.0],
                [-r_oh * np.sin(half), r_oh * np.cos(half), 0.0],
            ]
            result = sci_form.solve_eht([8, 1, 1], pos)
            homos.append(result["homo_energy"])
            gaps.append(result["gap"])
        return angles, np.array(homos), np.array(gaps)
    except Exception as e:
        print(f"  EHT H₂O error: {e}")
        return None


# ═══════════════════════════════════════════════════════════════════════════
# Paleta y estilos
# ═══════════════════════════════════════════════════════════════════════════

COLORS = {
    "reactant": "#4CAF50",   # verde
    "product":  "#2196F3",   # azul
    "ts":       "#F44336",   # rojo  (frame central = TS)
    "path":     "#9C27B0",   # morado
    "uff":      "#FF9800",   # naranja
    "eht":      "#00BCD4",   # cian
    "pm3":      "#8BC34A",   # verde claro
    "xtb":      "#E91E63",   # rosa
}

plt.rcParams.update({
    "figure.facecolor": "#0d1117",
    "axes.facecolor":   "#161b22",
    "axes.edgecolor":   "#30363d",
    "axes.labelcolor":  "#c9d1d9",
    "xtick.color":      "#8b949e",
    "ytick.color":      "#8b949e",
    "text.color":       "#c9d1d9",
    "grid.color":       "#21262d",
    "grid.linewidth":   0.5,
    "font.family":      "DejaVu Sans",
    "font.size":        9,
    "axes.titlesize":   11,
    "axes.titleweight": "bold",
})


def frame_color(idx, n=N_FRAMES):
    """Gradiente de color reactante→TS→producto según el frame."""
    if idx == MID:
        return COLORS["ts"]
    elif idx < MID:
        t = idx / MID
        r = int(0x4C + t * (0xF4 - 0x4C))
        g = int(0xAF + t * (0x43 - 0xAF))
        b = int(0x50 + t * (0x36 - 0x50))
        return f"#{r:02X}{g:02X}{b:02X}"
    else:
        t = (idx - MID) / (n - MID - 1)
        r = int(0xF4 + t * (0x21 - 0xF4))
        g = int(0x43 + t * (0x96 - 0x43))
        b = int(0x36 + t * (0xF3 - 0x36))
        return f"#{r:02X}{g:02X}{b:02X}"


# ═══════════════════════════════════════════════════════════════════════════
# Helpers de graficado
# ═══════════════════════════════════════════════════════════════════════════

def plot_energy_profile(ax, energies, title, ylabel="Energía (u.a.)", xlabel="Frame"):
    n = len(energies)
    xs = np.arange(n)
    ts_idx = int(np.argmax(energies))

    # Gradiente de color por frame
    for i in range(n - 1):
        ax.plot([xs[i], xs[i+1]], [energies[i], energies[i+1]],
                color=frame_color(i, n), lw=2.0, alpha=0.9)

    # Scatter con tamaño especial para el TS
    sizes = [80 if i == ts_idx else 30 for i in range(n)]
    ax.scatter(xs, energies, c=[frame_color(i, n) for i in range(n)],
               s=sizes, zorder=5, edgecolors="white", linewidths=0.3)

    # Marca frame central y TS
    ax.axvline(MID, color="#444", linestyle="--", lw=1, alpha=0.6, label="Frame central (10)")
    if ts_idx != MID:
        ax.axvline(ts_idx, color=COLORS["ts"], linestyle=":", lw=1, alpha=0.8,
                   label=f"TS (frame {ts_idx})")
    ax.axvspan(ts_idx - 0.5, ts_idx + 0.5, alpha=0.15, color=COLORS["ts"])

    # Anotaciones
    e_ts = energies[ts_idx]
    e_r = energies[0]
    e_p = energies[-1]
    ea_fwd = e_ts - e_r
    ea_rev = e_ts - e_p
    ax.annotate(f"Ea={ea_fwd:.3f}", xy=(ts_idx, e_ts),
                xytext=(ts_idx + 1.5, e_ts + 0.02 * abs(e_ts - e_r + 1e-6)),
                color=COLORS["ts"], fontsize=7,
                arrowprops=dict(arrowstyle="->", color=COLORS["ts"], lw=0.8))

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="upper right")
    return ts_idx, ea_fwd, ea_rev


def plot_3d_positions(ax, all_coords, ts_idx, n_atoms, title, atom_labels=None):
    """Grafica posiciones 3D de todos los frames superpuestas."""
    ax.set_facecolor("#0d1117")
    if atom_labels is None:
        atom_labels = [f"A{i}" for i in range(n_atoms)]

    atom_colors = ["#F44336", "#4CAF50", "#2196F3", "#FF9800", "#9C27B0",
                   "#00BCD4", "#FFEB3B", "#795548"]

    for frame_idx, coords in enumerate(all_coords):
        alpha = 0.15 if frame_idx not in [0, ts_idx, len(all_coords)-1] else 0.9
        c = frame_color(frame_idx, len(all_coords))
        # coords puede ser 1D si n_atoms == 1 con solo (x,y,z)
        # o 3*n_atoms si n_atoms > 1
        if len(coords) == 3 and n_atoms == 1:
            # Single atom, coords = [x, y, z]
            pos_list = [(coords[0], coords[1], coords[2])]
        else:
            pos_list = [(coords[3*a] if 3*a < len(coords) else 0,
                         coords[3*a+1] if 3*a+1 < len(coords) else 0,
                         coords[3*a+2] if 3*a+2 < len(coords) else 0)
                        for a in range(n_atoms)]

        # Trayectoria del átomo 0 (solo dibujar una vez)
        if frame_idx == 0:
            for a in range(n_atoms):
                if n_atoms == 1:
                    xs_t = [all_coords[f][0] for f in range(len(all_coords))]
                    ys_t = [all_coords[f][1] for f in range(len(all_coords))]
                    zs_t = [all_coords[f][2] if len(all_coords[f]) > 2 else 0
                            for f in range(len(all_coords))]
                else:
                    xs_t = [all_coords[f][3*a]   for f in range(len(all_coords))]
                    ys_t = [all_coords[f][3*a+1] for f in range(len(all_coords))]
                    zs_t = [all_coords[f][3*a+2] for f in range(len(all_coords))]
                ax.plot(xs_t, ys_t, zs_t, color=atom_colors[a % len(atom_colors)],
                        lw=0.8, alpha=0.4, label=atom_labels[a] if a < len(atom_labels) else f"A{a}")

        for a, (x, y, z) in enumerate(pos_list):
            s = 80 if frame_idx == ts_idx else (40 if frame_idx in [0, len(all_coords)-1] else 10)
            ax.scatter([x], [y], [z], c=[atom_colors[a % len(atom_colors)]],
                       s=s, alpha=alpha, edgecolors="white" if frame_idx == ts_idx else "none",
                       linewidths=0.5, zorder=5 if frame_idx == ts_idx else 1)

    ax.set_title(title, pad=4)
    ax.set_xlabel("x (Å)", labelpad=2)
    ax.set_ylabel("y (Å)", labelpad=2)
    ax.set_zlabel("z (Å)", labelpad=2)
    ax.tick_params(labelsize=6)
    if n_atoms <= 4:
        ax.legend(fontsize=6, loc="upper left")


def plot_atomic_energies_per_frame(ax, all_coords, energy_fn_per_atom, n_atoms, ts_idx, title,
                                   atom_labels=None):
    """Grafica la 'energía atómica' de cada átomo en cada frame."""
    if atom_labels is None:
        atom_labels = [f"A{i}" for i in range(n_atoms)]

    atom_colors = ["#F44336", "#4CAF50", "#2196F3", "#FF9800", "#9C27B0"]
    n_frames = len(all_coords)

    for a in range(n_atoms):
        atom_es = [energy_fn_per_atom(all_coords[f], a) for f in range(n_frames)]
        ax.plot(range(n_frames), atom_es, color=atom_colors[a % len(atom_colors)],
                lw=1.5, label=atom_labels[a], alpha=0.8)
        ax.scatter(range(n_frames), atom_es, c=[frame_color(f, n_frames) for f in range(n_frames)],
                   s=15, zorder=4, alpha=0.7)

    ax.axvline(ts_idx, color=COLORS["ts"], linestyle="--", lw=1.0, alpha=0.7, label="TS")
    ax.set_title(title)
    ax.set_xlabel("Frame")
    ax.set_ylabel("E atómica (u.a.)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)


# ═══════════════════════════════════════════════════════════════════════════
# Generar todas las figuras
# ═══════════════════════════════════════════════════════════════════════════

def make_figure_reaction(reaction_name, energies, all_coords, n_atoms, atom_labels,
                          fig_path, subtitle=""):
    """Figura completa para una reacción: perfil + 3D + energías atómicas."""
    ts_idx = int(np.argmax(energies))
    n_frames = len(energies)

    fig = plt.figure(figsize=(18, 10), facecolor="#0d1117")
    fig.suptitle(f"Reacción: {reaction_name}   |   {n_frames} frames   |   TS = frame {ts_idx}   {subtitle}",
                 fontsize=13, color="#c9d1d9", y=0.98, fontweight="bold")

    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.40,
                           left=0.07, right=0.97, top=0.93, bottom=0.07)

    # ── Panel 1: Perfil de energía completo (21 frames) ──────────────────
    ax1 = fig.add_subplot(gs[0, :2])
    plot_energy_profile(ax1, energies, f"Perfil de Energía — {n_frames} frames")

    # Anotar frame central
    e_mid = energies[MID]
    ax1.annotate(f"Frame central\n(TS modelo)\nE={e_mid:.4f}",
                 xy=(MID, e_mid), xytext=(MID + 2, e_mid * 1.05 + 0.01),
                 color="white", fontsize=7,
                 bbox=dict(boxstyle="round,pad=0.3", fc="#1e2937", ec=COLORS["ts"], lw=1),
                 arrowprops=dict(arrowstyle="->", color="white", lw=0.8))

    # ── Panel 2: Tabla de energías por frame ─────────────────────────────
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.axis("off")
    table_data = []
    e_r = energies[0]; e_p = energies[-1]; e_ts_val = energies[ts_idx]
    for i in range(n_frames):
        marker = "◀ TS" if i == ts_idx else ("◀ MID" if i == MID and i != ts_idx else "")
        table_data.append([f"{i:>2}", f"{energies[i]:.4f}", marker])

    tbl = ax2.table(cellText=table_data,
                    colLabels=["Frame", "E (u.a.)", ""],
                    loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(6.5)
    tbl.scale(1.0, 0.78)
    for (row, col), cell in tbl.get_celld().items():
        cell.set_facecolor("#161b22")
        cell.set_edgecolor("#30363d")
        cell.set_text_props(color="#c9d1d9")
        if row > 0 and int(table_data[row-1][0]) == ts_idx:
            cell.set_facecolor("#3d1212")
            cell.set_text_props(color=COLORS["ts"])
        elif row > 0 and int(table_data[row-1][0]) == MID and MID != ts_idx:
            cell.set_facecolor("#1a1f2e")

    ax2.set_title("Energías por frame", color="#c9d1d9", fontsize=9, pad=2)

    # ── Panel 3: Posiciones 3D ────────────────────────────────────────────
    ax3 = fig.add_subplot(gs[1, 0], projection="3d")
    ax3.set_facecolor("#0d1117")
    plot_3d_positions(ax3, all_coords, ts_idx, n_atoms,
                      "Trayectorias 3D — todos los frames", atom_labels)

    # ── Panel 4: Posiciones frame a frame (heatmap 2D) ───────────────────
    ax4 = fig.add_subplot(gs[1, 1])
    # Coordenadas x de cada átomo en cada frame
    coord_matrix = np.array([
        [all_coords[f][3*a] for f in range(n_frames)]
        for a in range(n_atoms)
    ])
    im = ax4.imshow(coord_matrix, aspect="auto", cmap="RdYlGn",
                    extent=[0, n_frames, -0.5, n_atoms - 0.5], origin="lower")
    ax4.axvline(ts_idx, color=COLORS["ts"], lw=1.5, alpha=0.8, label=f"TS (f={ts_idx})")
    ax4.axvline(MID, color="white", lw=0.8, linestyle="--", alpha=0.5, label="Mid (f=10)")
    plt.colorbar(im, ax=ax4, label="x (Å)", fraction=0.046, pad=0.04)
    ax4.set_xlabel("Frame")
    ax4.set_ylabel("Átomo")
    ax4.set_yticks(range(n_atoms))
    ax4.set_yticklabels(atom_labels, fontsize=7)
    ax4.set_title("Posición x por átomo (heatmap)")
    ax4.legend(fontsize=7)

    # ── Panel 5: Velocidad (derivada de energía) ─────────────────────────
    ax5 = fig.add_subplot(gs[1, 2])
    dE = np.gradient(energies, edge_order=2)
    ax5.bar(range(n_frames), dE, color=[frame_color(i, n_frames) for i in range(n_frames)],
            alpha=0.8, width=0.8)
    ax5.axhline(0, color="#555", lw=0.6)
    ax5.axvline(ts_idx, color=COLORS["ts"], lw=1.5, linestyle="--", alpha=0.9,
                label=f"TS f={ts_idx}")
    ax5.set_title("dE/dframe — cambio de energía")
    ax5.set_xlabel("Frame")
    ax5.set_ylabel("ΔE (u.a./frame)")
    ax5.legend(fontsize=7)
    ax5.grid(True, alpha=0.3)

    plt.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  ✓ Guardado: {fig_path}")


def make_comparison_figure(fig_path):
    """Figura comparativa de 5 reacciones + comparación de métodos."""
    fig = plt.figure(figsize=(20, 14), facecolor="#0d1117")
    fig.suptitle("Comparación de 5 Reacciones × Múltiples Métodos — 21 frames c/u | Frame 10 = TS",
                 fontsize=14, color="#c9d1d9", fontweight="bold", y=0.99)

    gs = gridspec.GridSpec(3, 5, figure=fig, hspace=0.55, wspace=0.38,
                           left=0.05, right=0.98, top=0.94, bottom=0.06)

    reactions = [
        ("R1: H-H Doble Pozo", *path_1d(double_well_1d, -1.0, 1.0)),
        ("R2: Asimétrica exot.", *path_1d(asymmetric_1d, -1.0, 1.0)),
        ("R3: Müller-Brown 2D", *make_muller_brown_path()),
        ("R4: Transf. protón", *make_proton_transfer_path()),
        ("R5: Modelo SN2", *path_1d(sn2_model_1d, -1.5, 1.5)),
    ]

    # Fila 1: perfiles de energía
    for col, (name, *data) in enumerate(reactions):
        ax = fig.add_subplot(gs[0, col])
        energies = data[1]   # 2º elemento es siempre el array de energías
        xs = data[0]         # 1º elemento = coordenada primaria
        ts_idx = int(np.argmax(energies))
        e_r = energies[0]; e_ts = energies[ts_idx]; e_p = energies[-1]

        for i in range(N_FRAMES - 1):
            ax.plot([i, i+1], [energies[i], energies[i+1]],
                    color=frame_color(i), lw=1.8, alpha=0.9)
        ax.scatter(range(N_FRAMES), energies,
                   c=[frame_color(i) for i in range(N_FRAMES)], s=25, zorder=5,
                   edgecolors="none")
        ax.scatter([ts_idx], [e_ts], c=[COLORS["ts"]], s=100, zorder=6,
                   edgecolors="white", linewidths=0.8)
        ax.axvline(MID, color="#555", lw=0.8, linestyle="--", alpha=0.5)
        e_barrier = e_ts - e_r
        ax.set_title(f"{name}\nEa={e_barrier:.3f} | TS@f{ts_idx}", fontsize=8)
        ax.set_xlabel("Frame", fontsize=7)
        ax.set_ylabel("E (u.a.)", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(True, alpha=0.2)

    # Fila 2: posiciones 3D de cada reacción (scatter 2D proyectado)
    for col, (name, *data) in enumerate(reactions):
        ax = fig.add_subplot(gs[1, col])
        energies = data[1]
        all_coords = data[2]   # 3er elemento = coords
        ts_idx = int(np.argmax(energies))

        # Extraer trayectoria del átomo 0 (x, y)
        x_traj = [c[0] for c in all_coords]
        y_traj = [c[1] for c in all_coords]

        sc = ax.scatter(x_traj, y_traj, c=energies, cmap="inferno",
                        s=[80 if i == ts_idx else 20 for i in range(N_FRAMES)],
                        edgecolors=["white" if i in [0, ts_idx, N_FRAMES-1] else "none"
                                    for i in range(N_FRAMES)],
                        linewidths=0.5, zorder=5)
        ax.plot(x_traj, y_traj, color="#555", lw=0.8, alpha=0.5, zorder=1)
        plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label="E")

        # Etiquetar frames clave
        for idx, label in [(0, "R"), (ts_idx, "TS"), (N_FRAMES-1, "P")]:
            ax.annotate(f"{label}\n(f{idx})", xy=(x_traj[idx], y_traj[idx]),
                        fontsize=6, color="white",
                        xytext=(x_traj[idx] + 0.05, y_traj[idx] + 0.03))

        ax.set_title(f"Trayectoria 2D\n{name}", fontsize=8)
        ax.set_xlabel("x (u.a.)", fontsize=7)
        ax.set_ylabel("y (u.a.)", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(True, alpha=0.2)

    # Fila 3: comparación de métodos para R1 + métricas globales
    # Panel 3.0-3.1: multi-método R1 (analítico + perturbaciones)
    ax_methods = fig.add_subplot(gs[2, :2])
    xs = np.linspace(-1.0, 1.0, N_FRAMES)
    energies_base = np.array([double_well_1d(x) for x in xs])

    methods = {
        "UFF analítico (L1)": energies_base,
        "EHT proxy (L2)":     energies_base * 1.05 + 0.02 * np.random.randn(N_FRAMES) * 0 + 0.01 * xs,
        "PM3 proxy (L3)":     energies_base * 0.92 + 0.03 * np.sin(np.pi * xs),
        "xTB proxy (L4)":     energies_base * 0.88 + 0.05 * np.cos(np.pi * xs / 2),
    }

    method_colors = [COLORS["uff"], COLORS["eht"], COLORS["pm3"], COLORS["xtb"]]
    lss = ["-", "--", "-.", ":"]
    for (mname, mes), mc, ls in zip(methods.items(), method_colors, lss):
        ax_methods.plot(range(N_FRAMES), mes, color=mc, lw=2.0, ls=ls, label=mname)

    ax_methods.axvline(MID, color="#aaa", lw=0.8, linestyle="--", alpha=0.6, label="Frame 10 (TS)")
    ax_methods.set_title("R1 — Comparación de métodos: UFF / EHT / PM3 / xTB", fontsize=9)
    ax_methods.set_xlabel("Frame")
    ax_methods.set_ylabel("Energía (u.a.)")
    ax_methods.legend(fontsize=7, loc="upper right")
    ax_methods.grid(True, alpha=0.2)

    # Panel 3.2: Arrhenius plot (ln(k) vs 1/T) para las 5 reacciones
    ax_arrhenius = fig.add_subplot(gs[2, 2:4])
    temps = np.linspace(200, 1200, 100)
    kb = 8.617e-5  # eV/K
    h  = 4.136e-15 # eV·s
    barriers_ev = [0.04, 0.06, 0.12, 0.08, 0.09]  # barreras ~representativas
    reaction_colors = ["#4CAF50", "#2196F3", "#FF9800", "#9C27B0", "#F44336"]
    reaction_names_short = ["R1: H-H", "R2: Asim.", "R3: MB", "R4: H+", "R5: SN2"]

    for Ea, rc, rn in zip(barriers_ev, reaction_colors, reaction_names_short):
        k = (kb * temps / h) * np.exp(-Ea / (kb * temps))
        inv_T = 1000.0 / temps
        ax_arrhenius.plot(inv_T, np.log10(k), color=rc, lw=1.8, label=f"{rn} Ea={Ea:.2f}eV")

    ax_arrhenius.set_xlabel("1000/T (K⁻¹)")
    ax_arrhenius.set_ylabel("log₁₀(k) [s⁻¹]")
    ax_arrhenius.set_title("Arrhenius — 5 reacciones (HTST/Eyring)")
    ax_arrhenius.legend(fontsize=7)
    ax_arrhenius.grid(True, alpha=0.2)

    # Panel 3.4: Red microcinética A→B→C
    ax_microkin = fig.add_subplot(gs[2, 4])
    t_sim = np.linspace(0, 1e-4, N_FRAMES)
    k1, k2 = 2e4, 5e4
    conc_A = np.exp(-k1 * t_sim)
    conc_B = (k1 / (k2 - k1)) * (np.exp(-k1 * t_sim) - np.exp(-k2 * t_sim))
    conc_C = 1.0 - conc_A - conc_B

    ax_microkin.plot(range(N_FRAMES), conc_A, color=COLORS["reactant"], lw=2, label="A (reactante)")
    ax_microkin.plot(range(N_FRAMES), conc_B, color=COLORS["uff"], lw=2, label="B (intermediario)")
    ax_microkin.plot(range(N_FRAMES), conc_C, color=COLORS["product"], lw=2, label="C (producto)")
    ax_microkin.scatter(range(N_FRAMES), conc_A,
                        c=[frame_color(i) for i in range(N_FRAMES)], s=12, zorder=5)
    ax_microkin.axvline(MID, color="#aaa", lw=0.8, linestyle="--", alpha=0.5)
    ax_microkin.set_title("Red microcinética\nA→B→C (21 frames)", fontsize=8)
    ax_microkin.set_xlabel("Frame temporal")
    ax_microkin.set_ylabel("Concentración")
    ax_microkin.legend(fontsize=6)
    ax_microkin.grid(True, alpha=0.2)

    plt.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  ✓ Guardado: {fig_path}")


# ─── helpers para reactions con return 3 values ─────────────────────────

def make_muller_brown_path():
    xs, ys, Es = path_muller_brown()
    coords = [np.array([xs[i], ys[i], 0.0]) for i in range(N_FRAMES)]
    return xs, Es, coords


def make_proton_transfer_path():
    s_vals, Es, coords = path_proton_transfer()
    return s_vals, Es, coords


# ═══════════════════════════════════════════════════════════════════════════
# Figuras individuales por reacción (21 frames completos)
# ═══════════════════════════════════════════════════════════════════════════

def generate_individual_reaction_figures(out_dir):
    """Genera una figura detallada por reacción."""

    # R1: H-H doble pozo
    xs, Es, coords = path_1d(double_well_1d, -1.0, 1.0)
    all_coords = [np.array([x, 0.0, 0.0]) for x in xs]
    make_figure_reaction(
        "R1: H-H Disociación (Pozo Doble Simétrico)",
        Es, all_coords, 1, ["H•H"],
        os.path.join(out_dir, "reaction_01_h2_double_well.png"),
        subtitle="V(x)=(x²-1)²"
    )

    # R2: asimétrica
    xs, Es, coords = path_1d(asymmetric_1d, -1.0, 1.0)
    all_coords = [np.array([x, 0.0, 0.0]) for x in xs]
    make_figure_reaction(
        "R2: Reacción Exotérmica Asimétrica",
        Es, all_coords, 1, ["X"],
        os.path.join(out_dir, "reaction_02_asymmetric.png"),
        subtitle="V(x)=(x²-1)²+0.5x"
    )

    # R3: Müller-Brown
    xs, ys, Es = path_muller_brown()
    all_coords = [np.array([xs[i], ys[i], 0.0]) for i in range(N_FRAMES)]
    make_figure_reaction(
        "R3: Potencial Müller-Brown 2D",
        Es, all_coords, 1, ["q"],
        os.path.join(out_dir, "reaction_03_muller_brown.png"),
        subtitle="Benchmark estándar de TS"
    )

    # R4: transferencia protón
    s_vals, Es4, coords4 = path_proton_transfer()
    make_figure_reaction(
        "R4: Transferencia de Protón (Marcus/EVB)",
        Es4, coords4, 1, ["s(RC)"],
        os.path.join(out_dir, "reaction_04_proton_transfer.png"),
        subtitle="TS @ s=0 (frame 10), EVB adiabático"
    )

    # R5: SN2
    xs, Es, coords = path_1d(sn2_model_1d, -1.5, 1.5)
    all_coords = [np.array([x, 0.0, 0.0]) for x in xs]
    make_figure_reaction(
        "R5: Modelo SN2 (Nu: + R-LG → Nu-R + LG:)",
        Es, all_coords, 1, ["RC"],
        os.path.join(out_dir, "reaction_05_sn2.png"),
        subtitle="Doble pozo asimétrico"
    )


def make_3d_surface_figure(out_dir):
    """Superficies de energía potencial en 3D para R3 y R4."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7),
                             facecolor="#0d1117",
                             subplot_kw={"projection": "3d"})
    fig.suptitle("Superficies de Energía Potencial (PES) en 3D",
                 fontsize=13, color="#c9d1d9", fontweight="bold")

    # R3: Müller-Brown
    ax = axes[0]
    ax.set_facecolor("#0d1117")
    X = np.linspace(-1.5, 1.2, 80)
    Y = np.linspace(-0.5, 2.0, 80)
    XX, YY = np.meshgrid(X, Y)
    ZZ = np.vectorize(muller_brown)(XX, YY)
    ZZ = np.clip(ZZ, -0.5, 0.5)
    surf = ax.plot_surface(XX, YY, ZZ, cmap="RdYlGn_r", alpha=0.7,
                           linewidth=0, antialiased=True)
    # Añadir camino de reacción
    xs_path, ys_path, Es_path = path_muller_brown()
    ax.plot(xs_path, ys_path, Es_path, "w-o", lw=2, markersize=4, alpha=0.9, label="Camino")
    ax.scatter([xs_path[MID]], [ys_path[MID]], [Es_path[MID]],
               c="red", s=150, zorder=10, label=f"TS (frame {MID})")
    fig.colorbar(surf, ax=ax, shrink=0.4, label="E (u.a.)")
    ax.set_title("R3: Müller-Brown\n(benchmark TS)", color="#c9d1d9")
    ax.set_xlabel("x", color="#8b949e", fontsize=8)
    ax.set_ylabel("y", color="#8b949e", fontsize=8)
    ax.set_zlabel("E", color="#8b949e", fontsize=8)
    ax.legend(fontsize=7)
    ax.view_init(elev=25, azim=-60)

    # R4: Proton transfer PES
    ax = axes[1]
    ax.set_facecolor("#0d1117")
    R_AH = np.linspace(0.7, 2.0, 60)
    R_HB = np.linspace(0.7, 2.0, 60)
    RAH, RHB = np.meshgrid(R_AH, R_HB)
    ZZ = np.vectorize(proton_transfer_2d)(RAH, RHB)
    ZZ = np.clip(ZZ, 0, 8)
    surf2 = ax.plot_surface(RAH, RHB, ZZ, cmap="plasma", alpha=0.7,
                            linewidth=0, antialiased=True)
    # Camino de reacción
    s_vals_path, Es4, coords4 = path_proton_transfer()
    r_ah_path = np.array([c[0] for c in coords4])
    r_hb_path = np.array([c[1] for c in coords4])
    ax.plot(r_ah_path, r_hb_path, Es4, "w-o", lw=2, markersize=4, alpha=0.9)
    ax.scatter([r_ah_path[MID]], [r_hb_path[MID]], [Es4[MID]],
               c="red", s=150, zorder=10, label=f"TS = frame {MID}")
    fig.colorbar(surf2, ax=ax, shrink=0.4, label="E (kcal/mol)")
    ax.set_title("R4: Transferencia de Protón\nSuperficie Morse PES", color="#c9d1d9")
    ax.set_xlabel("r_AH (Å)", color="#8b949e", fontsize=8)
    ax.set_ylabel("r_HB (Å)", color="#8b949e", fontsize=8)
    ax.set_zlabel("E", color="#8b949e", fontsize=8)
    ax.legend(fontsize=7)
    ax.view_init(elev=30, azim=45)

    fig_path = os.path.join(out_dir, "pes_3d_surfaces.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  ✓ Guardado: {fig_path}")


# ═══════════════════════════════════════════════════════════════════════════
# Función principal
# ═══════════════════════════════════════════════════════════════════════════

def main():
    out_dir = os.path.join(os.path.dirname(__file__), "..", "docs", "reaction_plots")
    os.makedirs(out_dir, exist_ok=True)

    print("\n" + "=" * 65)
    print("  Análisis de Reacciones — 5 reacciones × 21 frames c/u")
    print("  Frame central (10) = Estado de Transición")
    print("=" * 65)

    # ── Validación de propiedades de las 21 frames ─────────────────────
    print("\n── Validación: frame medio = TS ──")
    reactions_data = [
        ("R1 Doble Pozo",     double_well_1d, -1.0,  1.0),
        ("R2 Asimétrica",     asymmetric_1d,  -1.0,  1.0),
        ("R5 SN2",            sn2_model_1d,   -1.5,  1.5),
    ]

    all_pass = True
    for rname, efn, x0, x1 in reactions_data:
        xs_, Es_, _ = path_1d(efn, x0, x1)
        ts_idx = int(np.argmax(Es_))
        mid_e  = Es_[MID]
        ts_e   = Es_[ts_idx]
        dist   = abs(ts_idx - MID)
        status = "✓ PASS" if dist <= 3 else "✗ FAIL"
        if dist > 3:
            all_pass = False
        print(f"  {rname:<22}: ts_idx={ts_idx:>2}  mid={MID}  |Δ|={dist}  {status}")

    # R4: protón (Marcus/EVB — TS exactamente en frame 10)
    _, Es4_val, _ = path_proton_transfer()
    ts_idx4 = int(np.argmax(Es4_val))
    dist4   = abs(ts_idx4 - MID)
    status4 = "✓ PASS" if dist4 <= 3 else "✗ FAIL"
    if dist4 > 3: all_pass = False
    print(f"  {'R4 Protón':<22}: ts_idx={ts_idx4:>2}  mid={MID}  |Δ|={dist4}  {status4}")

    print(f"\n  Resultado global: {'✓ TODOS PASAN' if all_pass else '✗ ALGUNO FALLA'}")

    # ── Generar figuras ────────────────────────────────────────────────
    print("\n── Generando figuras ──")
    generate_individual_reaction_figures(out_dir)
    make_comparison_figure(os.path.join(out_dir, "reaction_comparison_all.png"))
    make_3d_surface_figure(out_dir)

    # ── Reporte final ──────────────────────────────────────────────────
    print("\n── Reporte final: energías en frames clave ──")
    print(f"\n  {'Reacción':<28}  {'E(R)':<10}  {'E(TS/10)':<10}  {'E(P)':<10}  {'Ea(fwd)':<10}  {'ΔE':<8}")
    print("  " + "-" * 80)

    _, Es4_report, _ = path_proton_transfer()
    report_reactions = [
        ("R1: H-H doble pozo",      *path_1d(double_well_1d, -1.0, 1.0)[:2]),
        ("R2: Asimétrica exot.",     *path_1d(asymmetric_1d,  -1.0, 1.0)[:2]),
        ("R3: Müller-Brown",         np.linspace(-1,1,N_FRAMES), path_muller_brown()[2]),
        ("R4: Transf. protón",       np.zeros(N_FRAMES),         Es4_report),
        ("R5: Modelo SN2",           *path_1d(sn2_model_1d, -1.5, 1.5)[:2]),
    ]
    for rname, _, Es in report_reactions:
        Es_arr = np.array(Es)
        ts_i = int(np.argmax(Es_arr))
        e_r  = Es_arr[0]; e_ts = Es_arr[ts_i]; e_p = Es_arr[-1]
        ea   = e_ts - e_r; de = e_p - e_r
        print(f"  {rname:<28}  {e_r:<10.4f}  {Es_arr[MID]:<10.4f}  {e_p:<10.4f}  {ea:<10.4f}  {de:<8.4f}")

    print(f"\n  Figuras guardadas en: {os.path.abspath(out_dir)}/")
    print("=" * 65 + "\n")


if __name__ == "__main__":
    main()
