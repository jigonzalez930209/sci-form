"""Spectroscopy: UV-vis sTDA, vibrational/IR, DOS already elsewhere."""
from sci_form import embed, stda_uvvis, vibrational_analysis, ir_spectrum

from common import require_ok


def run():
    print("=== 06 spectroscopy ===")
    conf = require_ok(embed("c1ccccc1", seed=42), "embed benzene")

    spec = stda_uvvis(conf.elements, conf.coords, sigma=0.3, e_min=1.0, e_max=8.0, n_points=200)
    print("wavelengths", len(spec.wavelengths_nm), "absorptivity", len(spec.absorptivity), "excitations", len(spec.excitations))
    if spec.excitations:
        ex = spec.excitations[0]
        print("exc", ex.wavelength_nm, ex.oscillator_strength, ex.from_mo, ex.to_mo)

    vib = vibrational_analysis(conf.elements, conf.coords, method="eht")
    print("vib atoms", vib.n_atoms, "real modes", vib.n_real_modes, "zpve", vib.zpve_ev)
    if vib.modes:
        m = vib.modes[0]
        print("mode", m.frequency_cm1, m.ir_intensity, m.is_real, len(m.displacement))

    ir = ir_spectrum(conf.elements, conf.coords, method="eht")
    print("IR wavenumbers", len(ir.wavenumbers), "peaks", len(ir.peaks), "gamma", ir.gamma)
    if ir.peaks:
        pk = ir.peaks[0]
        print("peak", pk)

    # Broadened IR helper is not exposed in this wheel; reuse ir spectrum instead
    print("IR peaks", len(ir.peaks), "wavenumbers", len(ir.wavenumbers))


if __name__ == "__main__":
    run()
