"""NMR: shifts, couplings, spectra, GIAO, HOSE."""
from sci_form import embed, nmr_shifts, nmr_shifts_for_nucleus, nmr_couplings, nmr_spectrum, hose_codes, giao_nmr


def run():
    print("=== 07 nmr ===")
    shifts = nmr_shifts("CCO")
    print("h_shifts", len(shifts.h_shifts), "c_shifts", len(shifts.c_shifts))
    if shifts.h_shifts:
        h = shifts.h_shifts[0]
        print("H", h.atom_index, h.shift_ppm, h.environment, h.confidence)

    pt = nmr_shifts_for_nucleus("[Pt]", "195Pt")
    print("Pt shifts", len(pt))

    conf = embed("CC", seed=42)
    couplings = nmr_couplings("CC", conf.coords)
    print("couplings", len(couplings))
    if couplings:
        jc = couplings[0]
        print("J", jc.h1_index, jc.h2_index, jc.j_hz, jc.n_bonds, jc.coupling_type)

    # Ensemble couplings are not exposed; skip

    spec = nmr_spectrum("CCO", nucleus="1H", gamma=0.01, ppm_min=-2.0, ppm_max=14.0, n_points=500)
    print("nmr spectrum", spec.nucleus, len(spec.ppm_axis), len(spec.intensities), len(spec.peaks))

    codes = hose_codes("CCO", max_radius=4)
    print("hose codes", len(codes))
    print("sample", codes[:3])

    water = [8, 1, 1]
    wcoords = [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692]
    giao = giao_nmr(water, wcoords, nucleus="1H")
    print("giao", giao.chemical_shifts, giao.scf_converged, giao.scf_iterations, giao.fallback_elements, giao.notes)


if __name__ == "__main__":
    run()
