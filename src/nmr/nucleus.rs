//! NMR-active nuclides and the reference data attached to them.
//!
//! This module defines [`NmrNucleus`], the crate-wide handle on "which atom
//! is being observed", together with the tabulated reference data that the
//! rest of the `nmr` module needs:
//!
//! * the Larmor frequency of each nuclide ([`NmrNucleus::default_frequency_mhz`],
//!   in MHz),
//! * the natural (unbroadened) resonance line width
//!   ([`NmrNucleus::default_fwhm_hz`], in Hz),
//! * the empirical chemical-shift window and centre of each nuclide
//!   ([`NmrNucleus::empirical_range`], [`NmrNucleus::empirical_center`], in ppm),
//! * a reference shielding constant ([`NmrNucleus::reference_shielding`], in ppm)
//!   and the free-atom diamagnetic contribution it is built from
//!   ([`NmrNucleus::free_atom_diamagnetic`], in ppm).
//!
//! The frequencies and line widths are the conventional values quoted for a
//! standard high-field spectrometer; they are used to convert between Hz and
//! ppm when a spectrum is simulated (see `nmr::spectrum`), not to predict the
//! field-dependent shift of a resonance.
//!
//! The nuclide set covers the isotopes that are actually observed in routine
//! structure determination, plus a catalogue of NMR-active metal isotopes
//! that are not worth a dedicated enum variant.

use serde::{Deserialize, Serialize};

/// An NMR-active nuclide.
///
/// Each variant identifies one isotope with a non-zero nuclear spin. The
/// variant itself carries no payload: the physical reference data (Larmor
/// frequency in MHz, natural line width in Hz, chemical-shift window and
/// centre in ppm, reference shielding in ppm) is exposed through the
/// associated methods rather than stored per instance.
///
/// The ordering of the variants is the declaration order, so the derived
/// `Ord` implementation sorts by that order. It is *not* the periodic-table
/// order: ⁵⁵Mn is declared last, whereas [`NmrNucleus::ALL`] lists the same
/// nuclides in element order.
///
/// # Example
///
/// ```
/// use sci_form::nmr::nucleus::NmrNucleus;
///
/// let nucleus = NmrNucleus::parse("13C").unwrap();
/// assert_eq!(nucleus.atomic_number(), 6);
/// assert_eq!(nucleus.mass_number(), 13);
/// assert_eq!(nucleus.canonical(), "13C");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
// The variant set may grow in a future release; downstream crates must
// match with a wildcard arm.
#[non_exhaustive]
pub enum NmrNucleus {
    /// Protium (¹H): spin ½, the dominant target of ¹H NMR.
    H1,
    /// Deuterium (²H): spin 1 and therefore quadrupolar; the standard ¹H lock solvent.
    H2,
    /// Tritium (³H): spin ½, a radioactive hydrogen label.
    H3,
    /// Helium-3 (³He): spin ½, used for low-temperature polarised targets.
    He3,
    /// Lithium-6 (⁶Li): spin 1 and quadrupolar.
    Li6,
    /// Lithium-7 (⁷Li): spin 3/2 and quadrupolar.
    Li7,
    /// Beryllium-9 (⁹Be): spin 3/2 and quadrupolar.
    Be9,
    /// Boron-10 (¹⁰B): spin 3 and strongly quadrupolar.
    B10,
    /// Boron-11 (¹¹B): spin 3/2 and quadrupolar; the usual boron nucleus.
    B11,
    /// Carbon-13 (¹³C): spin ½, the standard carbon nucleus.
    C13,
    /// Nitrogen-14 (¹⁴N): spin 1 and quadrupolar; the abundant nitrogen isotope.
    N14,
    /// Nitrogen-15 (¹⁵N): spin ½, the enrichment-friendly nitrogen nucleus.
    N15,
    /// Oxygen-17 (¹⁷O): spin 5/2 and quadrupolar.
    O17,
    /// Fluorine-19 (¹⁹F): spin ½ and 100 % naturally abundant.
    F19,
    /// Sodium-23 (²³Na): spin 3/2 and quadrupolar.
    Na23,
    /// Magnesium-25 (²⁵Mg): spin 5/2 and quadrupolar.
    Mg25,
    /// Aluminium-27 (²⁷Al): spin 5/2 and quadrupolar.
    Al27,
    /// Silicon-29 (²⁹Si): spin ½.
    Si29,
    /// Phosphorus-31 (³¹P): spin ½ and 100 % naturally abundant.
    P31,
    /// Sulfur-33 (³³S): spin 3/2 and quadrupolar.
    S33,
    /// Chlorine-35 (³⁵Cl): spin 3/2 and quadrupolar.
    Cl35,
    /// Chlorine-37 (³⁷Cl): spin 3/2 and quadrupolar.
    Cl37,
    /// Potassium-39 (³⁹K): spin 3/2 and quadrupolar.
    K39,
    /// Potassium-40 (⁴⁰K): spin 4 and strongly quadrupolar; the long-lived radioactive isotope.
    K40,
    /// Potassium-41 (⁴¹K): spin 3/2 and quadrupolar.
    K41,
    /// Calcium-43 (⁴³Ca): spin 7/2 and quadrupolar.
    Ca43,
    /// Scandium-45 (⁴⁵Sc): spin 7/2 and quadrupolar.
    Sc45,
    /// Titanium-47 (⁴⁷Ti): spin 5/2 and quadrupolar.
    Ti47,
    /// Titanium-49 (⁴⁹Ti): spin 7/2 and quadrupolar.
    Ti49,
    /// Vanadium-50 (⁵⁰V): spin 6 and strongly quadrupolar.
    V50,
    /// Vanadium-51 (⁵¹V): spin 7/2 and quadrupolar.
    V51,
    /// Chromium-53 (⁵³Cr): spin 3/2 and quadrupolar.
    Cr53,
    /// Iron-57 (⁵⁷Fe): spin ½.
    Fe57,
    /// Cobalt-59 (⁵⁹Co): spin 7/2 and quadrupolar.
    Co59,
    /// Nickel-61 (⁶¹Ni): spin 3/2 and quadrupolar.
    Ni61,
    /// Copper-63 (⁶³Cu): spin 3/2 and quadrupolar.
    Cu63,
    /// Copper-65 (⁶⁵Cu): spin 3/2 and quadrupolar.
    Cu65,
    /// Zinc-67 (⁶⁷Zn): spin 5/2 and quadrupolar.
    Zn67,
    /// Gallium-69 (⁶⁹Ga): spin 3/2 and quadrupolar.
    Ga69,
    /// Gallium-71 (⁷¹Ga): spin 3/2 and quadrupolar.
    Ga71,
    /// Germanium-73 (⁷³Ge): spin 9/2 and quadrupolar.
    Ge73,
    /// Arsenic-75 (⁷⁵As): spin 3/2 and quadrupolar.
    As75,
    /// Selenium-77 (⁷⁷Se): spin ½.
    Se77,
    /// Bromine-79 (⁷⁹Br): spin 3/2 and quadrupolar.
    Br79,
    /// Bromine-81 (⁸¹Br): spin 3/2 and quadrupolar.
    Br81,
    /// Rubidium-85 (⁸⁵Rb): spin 5/2 and quadrupolar.
    Rb85,
    /// Rubidium-87 (⁸⁷Rb): spin 3/2 and quadrupolar.
    Rb87,
    /// Strontium-87 (⁸⁷Sr): spin 9/2 and quadrupolar.
    Sr87,
    /// Zirconium-91 (⁹¹Zr): spin 5/2 and quadrupolar.
    Zr91,
    /// Niobium-93 (⁹³Nb): spin 9/2 and quadrupolar.
    Nb93,
    /// Molybdenum-95 (⁹⁵Mo): spin 5/2 and quadrupolar.
    Mo95,
    /// Molybdenum-97 (⁹⁷Mo): spin 5/2 and quadrupolar.
    Mo97,
    /// Ruthenium-99 (⁹⁹Ru): spin 5/2 and quadrupolar.
    Ru99,
    /// Ruthenium-101 (¹⁰¹Ru): spin 5/2 and quadrupolar.
    Ru101,
    /// Rhodium-103 (¹⁰³Rh): spin ½.
    Rh103,
    /// Palladium-105 (¹⁰⁵Pd): spin 5/2 and quadrupolar.
    Pd105,
    /// Silver-107 (¹⁰⁷Ag): spin ½.
    Ag107,
    /// Silver-109 (¹⁰⁹Ag): spin ½.
    Ag109,
    /// Cadmium-111 (¹¹¹Cd): spin ½.
    Cd111,
    /// Cadmium-113 (¹¹³Cd): spin ½.
    Cd113,
    /// Indium-113 (¹¹³In): spin 9/2 and quadrupolar.
    In113,
    /// Indium-115 (¹¹⁵In): spin 9/2 and quadrupolar.
    In115,
    /// Tin-115 (¹¹⁵Sn): spin ½.
    Sn115,
    /// Tin-117 (¹¹⁷Sn): spin ½.
    Sn117,
    /// Tin-119 (¹¹⁹Sn): spin ½.
    Sn119,
    /// Antimony-121 (¹²¹Sb): spin 5/2 and quadrupolar.
    Sb121,
    /// Antimony-123 (¹²³Sb): spin 7/2 and quadrupolar.
    Sb123,
    /// Tellurium-123 (¹²³Te): spin ½.
    Te123,
    /// Tellurium-125 (¹²⁵Te): spin ½.
    Te125,
    /// Iodine-127 (¹²⁷I): spin 5/2 and quadrupolar.
    I127,
    /// Xenon-129 (¹²⁹Xe): spin ½.
    Xe129,
    /// Xenon-131 (¹³¹Xe): spin 3/2 and quadrupolar.
    Xe131,
    /// Cesium-133 (¹³³Cs): spin 7/2 and quadrupolar.
    Cs133,
    /// Barium-135 (¹³⁵Ba): spin 3/2 and quadrupolar.
    Ba135,
    /// Barium-137 (¹³⁷Ba): spin 3/2 and quadrupolar.
    Ba137,
    /// Tungsten-183 (¹⁸³W): spin ½.
    W183,
    /// Platinum-195 (¹⁹⁵Pt): spin ½.
    Pt195,
    /// Gold-197 (¹⁹⁷Au): spin 3/2 and quadrupolar.
    Au197,
    /// Mercury-199 (¹⁹⁹Hg): spin ½.
    Hg199,
    /// Mercury-201 (²⁰¹Hg): spin 3/2 and quadrupolar.
    Hg201,
    /// Thallium-203 (²⁰³Tl): spin ½.
    Tl203,
    /// Thallium-205 (²⁰⁵Tl): spin ½.
    Tl205,
    /// Lead-207 (²⁰⁷Pb): spin ½.
    Pb207,
    /// Bismuth-209 (²⁰⁹Bi): spin 9/2 and quadrupolar.
    Bi209,
    /// Manganese-55 (⁵⁵Mn): spin 5/2 and quadrupolar.
    Mn55,
    /// NMR-active metal isotope stored in `METAL_CATALOG`.
    ///
    /// The payload is the zero-based index into the built-in metal catalogue
    /// (the `METAL_CATALOG` table), not an atomic number. Use
    /// [`NmrNucleus::atomic_number`] and [`NmrNucleus::mass_number`] to
    /// recover the element and isotope.
    Catalog(u8),
}

struct MetalCatalogEntry {
    z: u8,
    mass: u16,
    symbol: &'static str,
    label: &'static str,
    pretty: &'static str,
    freq_mhz: f64,
    /// Nuclear spin I. Values above 1/2 are quadrupolar.
    spin_i: f64,
    range: (f64, f64),
    center: f64,
    fwhm_hz: f64,
}

/// Default NMR-active isotope for metals that are not already a unit variant.
/// Cerium is omitted: every natural Ce isotope has I = 0.
const METAL_CATALOG: &[MetalCatalogEntry] = &[
    MetalCatalogEntry {
        z: 39,
        mass: 89,
        symbol: "Y",
        label: "89Y",
        pretty: "⁸⁹Y",
        freq_mhz: 19.6,
        spin_i: 0.5,
        range: (-1500.0, 1500.0),
        center: 0.0,
        fwhm_hz: 25.0,
    },
    MetalCatalogEntry {
        z: 43,
        mass: 99,
        symbol: "Tc",
        label: "99Tc",
        pretty: "⁹⁹Tc",
        freq_mhz: 90.1,
        spin_i: 4.5,
        range: (-5000.0, 8000.0),
        center: 0.0,
        fwhm_hz: 400.0,
    },
    MetalCatalogEntry {
        z: 57,
        mass: 139,
        symbol: "La",
        label: "139La",
        pretty: "¹³⁹La",
        freq_mhz: 56.5,
        spin_i: 3.5,
        range: (-1500.0, 1500.0),
        center: 0.0,
        fwhm_hz: 400.0,
    },
    MetalCatalogEntry {
        z: 59,
        mass: 141,
        symbol: "Pr",
        label: "141Pr",
        pretty: "¹⁴¹Pr",
        freq_mhz: 11.7,
        spin_i: 2.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 60,
        mass: 143,
        symbol: "Nd",
        label: "143Nd",
        pretty: "¹⁴³Nd",
        freq_mhz: 22.0,
        spin_i: 3.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 61,
        mass: 147,
        symbol: "Pm",
        label: "147Pm",
        pretty: "¹⁴⁷Pm",
        freq_mhz: 22.0,
        spin_i: 3.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 62,
        mass: 149,
        symbol: "Sm",
        label: "149Sm",
        pretty: "¹⁴⁹Sm",
        freq_mhz: 14.0,
        spin_i: 3.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 63,
        mass: 151,
        symbol: "Eu",
        label: "151Eu",
        pretty: "¹⁵¹Eu",
        freq_mhz: 10.0,
        spin_i: 2.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 64,
        mass: 155,
        symbol: "Gd",
        label: "155Gd",
        pretty: "¹⁵⁵Gd",
        freq_mhz: 16.0,
        spin_i: 1.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 800.0,
    },
    MetalCatalogEntry {
        z: 65,
        mass: 159,
        symbol: "Tb",
        label: "159Tb",
        pretty: "¹⁵⁹Tb",
        freq_mhz: 24.0,
        spin_i: 1.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 66,
        mass: 163,
        symbol: "Dy",
        label: "163Dy",
        pretty: "¹⁶³Dy",
        freq_mhz: 18.0,
        spin_i: 2.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 67,
        mass: 165,
        symbol: "Ho",
        label: "165Ho",
        pretty: "¹⁶⁵Ho",
        freq_mhz: 8.7,
        spin_i: 3.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 68,
        mass: 167,
        symbol: "Er",
        label: "167Er",
        pretty: "¹⁶⁷Er",
        freq_mhz: 11.5,
        spin_i: 3.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 69,
        mass: 169,
        symbol: "Tm",
        label: "169Tm",
        pretty: "¹⁶⁹Tm",
        freq_mhz: 33.0,
        spin_i: 0.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 40.0,
    },
    MetalCatalogEntry {
        z: 70,
        mass: 171,
        symbol: "Yb",
        label: "171Yb",
        pretty: "¹⁷¹Yb",
        freq_mhz: 70.0,
        spin_i: 0.5,
        range: (-20000.0, 20000.0),
        center: 0.0,
        fwhm_hz: 40.0,
    },
    MetalCatalogEntry {
        z: 71,
        mass: 175,
        symbol: "Lu",
        label: "175Lu",
        pretty: "¹⁷⁵Lu",
        freq_mhz: 48.0,
        spin_i: 3.5,
        range: (-2000.0, 2000.0),
        center: 0.0,
        fwhm_hz: 400.0,
    },
    MetalCatalogEntry {
        z: 72,
        mass: 177,
        symbol: "Hf",
        label: "177Hf",
        pretty: "¹⁷⁷Hf",
        freq_mhz: 12.8,
        spin_i: 3.5,
        range: (-3000.0, 3000.0),
        center: 0.0,
        fwhm_hz: 400.0,
    },
    MetalCatalogEntry {
        z: 73,
        mass: 181,
        symbol: "Ta",
        label: "181Ta",
        pretty: "¹⁸¹Ta",
        freq_mhz: 48.0,
        spin_i: 3.5,
        range: (-4000.0, 4000.0),
        center: 0.0,
        fwhm_hz: 600.0,
    },
    MetalCatalogEntry {
        z: 75,
        mass: 187,
        symbol: "Re",
        label: "187Re",
        pretty: "¹⁸⁷Re",
        freq_mhz: 91.0,
        spin_i: 2.5,
        range: (-5000.0, 5000.0),
        center: 0.0,
        fwhm_hz: 400.0,
    },
    MetalCatalogEntry {
        z: 76,
        mass: 187,
        symbol: "Os",
        label: "187Os",
        pretty: "¹⁸⁷Os",
        freq_mhz: 9.2,
        spin_i: 0.5,
        range: (-6000.0, 4000.0),
        center: 0.0,
        fwhm_hz: 30.0,
    },
    MetalCatalogEntry {
        z: 77,
        mass: 193,
        symbol: "Ir",
        label: "193Ir",
        pretty: "¹⁹³Ir",
        freq_mhz: 8.0,
        spin_i: 1.5,
        range: (-6000.0, 6000.0),
        center: 0.0,
        fwhm_hz: 400.0,
    },
    MetalCatalogEntry {
        z: 87,
        mass: 223,
        symbol: "Fr",
        label: "223Fr",
        pretty: "²²³Fr",
        freq_mhz: 88.0,
        spin_i: 1.5,
        range: (-500.0, 500.0),
        center: 0.0,
        fwhm_hz: 200.0,
    },
    MetalCatalogEntry {
        z: 88,
        mass: 223,
        symbol: "Ra",
        label: "223Ra",
        pretty: "²²³Ra",
        freq_mhz: 28.0,
        spin_i: 1.5,
        range: (-500.0, 500.0),
        center: 0.0,
        fwhm_hz: 200.0,
    },
    MetalCatalogEntry {
        z: 89,
        mass: 227,
        symbol: "Ac",
        label: "227Ac",
        pretty: "²²⁷Ac",
        freq_mhz: 56.0,
        spin_i: 1.5,
        range: (-2000.0, 2000.0),
        center: 0.0,
        fwhm_hz: 400.0,
    },
    MetalCatalogEntry {
        z: 90,
        mass: 229,
        symbol: "Th",
        label: "229Th",
        pretty: "²²⁹Th",
        freq_mhz: 14.0,
        spin_i: 2.5,
        range: (-2000.0, 2000.0),
        center: 0.0,
        fwhm_hz: 400.0,
    },
    MetalCatalogEntry {
        z: 91,
        mass: 231,
        symbol: "Pa",
        label: "231Pa",
        pretty: "²³¹Pa",
        freq_mhz: 50.0,
        spin_i: 1.5,
        range: (-4000.0, 4000.0),
        center: 0.0,
        fwhm_hz: 400.0,
    },
    MetalCatalogEntry {
        z: 92,
        mass: 235,
        symbol: "U",
        label: "235U",
        pretty: "²³⁵U",
        freq_mhz: 7.6,
        spin_i: 3.5,
        range: (-4000.0, 4000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 93,
        mass: 237,
        symbol: "Np",
        label: "237Np",
        pretty: "²³⁷Np",
        freq_mhz: 45.0,
        spin_i: 2.5,
        range: (-6000.0, 6000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
    MetalCatalogEntry {
        z: 94,
        mass: 239,
        symbol: "Pu",
        label: "239Pu",
        pretty: "²³⁹Pu",
        freq_mhz: 15.0,
        spin_i: 0.5,
        range: (-4000.0, 4000.0),
        center: 0.0,
        fwhm_hz: 40.0,
    },
    MetalCatalogEntry {
        z: 95,
        mass: 243,
        symbol: "Am",
        label: "243Am",
        pretty: "²⁴³Am",
        freq_mhz: 24.0,
        spin_i: 2.5,
        range: (-8000.0, 8000.0),
        center: 0.0,
        fwhm_hz: 500.0,
    },
];

fn catalog_entry(index: u8) -> &'static MetalCatalogEntry {
    &METAL_CATALOG[index as usize]
}

impl NmrNucleus {
    /// Every nuclide that has a dedicated enum variant.
    ///
    /// The array holds the 85 non-catalogue nuclides in periodic-table order
    /// (element by element, ascending atomic number); the catalogue metals
    /// are not included because they are addressed through
    /// [`NmrNucleus::Catalog`]. Use [`NmrNucleus::supported_labels`] for the
    /// full list of human-readable labels including the catalogue entries.
    pub const ALL: [NmrNucleus; 85] = [
        NmrNucleus::H1,
        NmrNucleus::H2,
        NmrNucleus::H3,
        NmrNucleus::He3,
        NmrNucleus::Li6,
        NmrNucleus::Li7,
        NmrNucleus::Be9,
        NmrNucleus::B10,
        NmrNucleus::B11,
        NmrNucleus::C13,
        NmrNucleus::N14,
        NmrNucleus::N15,
        NmrNucleus::O17,
        NmrNucleus::F19,
        NmrNucleus::Na23,
        NmrNucleus::Mg25,
        NmrNucleus::Al27,
        NmrNucleus::Si29,
        NmrNucleus::P31,
        NmrNucleus::S33,
        NmrNucleus::Cl35,
        NmrNucleus::Cl37,
        NmrNucleus::K39,
        NmrNucleus::K40,
        NmrNucleus::K41,
        NmrNucleus::Ca43,
        NmrNucleus::Sc45,
        NmrNucleus::Ti47,
        NmrNucleus::Ti49,
        NmrNucleus::V50,
        NmrNucleus::V51,
        NmrNucleus::Cr53,
        NmrNucleus::Mn55,
        NmrNucleus::Fe57,
        NmrNucleus::Co59,
        NmrNucleus::Ni61,
        NmrNucleus::Cu63,
        NmrNucleus::Cu65,
        NmrNucleus::Zn67,
        NmrNucleus::Ga69,
        NmrNucleus::Ga71,
        NmrNucleus::Ge73,
        NmrNucleus::As75,
        NmrNucleus::Se77,
        NmrNucleus::Br79,
        NmrNucleus::Br81,
        NmrNucleus::Rb85,
        NmrNucleus::Rb87,
        NmrNucleus::Sr87,
        NmrNucleus::Zr91,
        NmrNucleus::Nb93,
        NmrNucleus::Mo95,
        NmrNucleus::Mo97,
        NmrNucleus::Ru99,
        NmrNucleus::Ru101,
        NmrNucleus::Rh103,
        NmrNucleus::Pd105,
        NmrNucleus::Ag107,
        NmrNucleus::Ag109,
        NmrNucleus::Cd111,
        NmrNucleus::Cd113,
        NmrNucleus::In113,
        NmrNucleus::In115,
        NmrNucleus::Sn115,
        NmrNucleus::Sn117,
        NmrNucleus::Sn119,
        NmrNucleus::Sb121,
        NmrNucleus::Sb123,
        NmrNucleus::Te123,
        NmrNucleus::Te125,
        NmrNucleus::I127,
        NmrNucleus::Xe129,
        NmrNucleus::Xe131,
        NmrNucleus::Cs133,
        NmrNucleus::Ba135,
        NmrNucleus::Ba137,
        NmrNucleus::W183,
        NmrNucleus::Pt195,
        NmrNucleus::Au197,
        NmrNucleus::Hg199,
        NmrNucleus::Hg201,
        NmrNucleus::Tl203,
        NmrNucleus::Tl205,
        NmrNucleus::Pb207,
        NmrNucleus::Bi209,
    ];

    /// Parse a nuclide from a textual alias.
    ///
    /// The alias is normalised by dropping every non-ASCII-alphanumeric
    /// character and lower-casing what remains, so `"13C"`, `"C13"` and
    /// `"carbon"` all resolve. Accepted spellings are the mass number in
    /// either order, the element name for the common isotopes, and the
    /// catalogue metal labels (both `"89Y"` and `"Y89"`).
    ///
    /// Note that the normalisation is ASCII-only: a Unicode label such as
    /// `"¹³C"` has its superscript mass number stripped and is *not*
    /// accepted. Use [`NmrNucleus::unicode_label`] to *print* a Unicode
    /// label, not to parse one.
    ///
    /// # Errors
    ///
    /// Returns `Err` with a human-readable message listing the supported
    /// labels when the alias does not resolve to a known nuclide.
    ///
    /// # Example
    ///
    /// ```
    /// use sci_form::nmr::nucleus::NmrNucleus;
    ///
    /// assert_eq!(NmrNucleus::parse("proton").unwrap(), NmrNucleus::H1);
    /// assert!(NmrNucleus::parse("not-an-isotope").is_err());
    /// ```
    pub fn parse(alias: &str) -> Result<Self, String> {
        let normalized: String = alias
            .chars()
            .filter(|character| character.is_ascii_alphanumeric())
            .flat_map(|character| character.to_lowercase())
            .collect();

        let nucleus = match normalized.as_str() {
            "1h" | "h1" | "proton" | "hydrogen" => Self::H1,
            "2h" | "h2" | "deuterium" => Self::H2,
            "3h" | "h3" | "tritium" => Self::H3,
            "3he" | "he3" | "helium" => Self::He3,
            "6li" | "li6" => Self::Li6,
            "7li" | "li7" | "lithium" => Self::Li7,
            "9be" | "be9" | "beryllium" => Self::Be9,
            "10b" | "b10" => Self::B10,
            "11b" | "b11" | "boron" => Self::B11,
            "13c" | "c13" | "carbon" => Self::C13,
            "14n" | "n14" => Self::N14,
            "15n" | "n15" | "nitrogen" => Self::N15,
            "17o" | "o17" | "oxygen" => Self::O17,
            "19f" | "f19" | "fluorine" => Self::F19,
            "23na" | "na23" | "sodium" => Self::Na23,
            "25mg" | "mg25" | "magnesium" => Self::Mg25,
            "27al" | "al27" | "aluminium" | "aluminum" => Self::Al27,
            "29si" | "si29" | "silicon" => Self::Si29,
            "31p" | "p31" | "phosphorus" => Self::P31,
            "33s" | "s33" | "sulfur" | "sulphur" => Self::S33,
            "35cl" | "cl35" | "chlorine" => Self::Cl35,
            "37cl" | "cl37" => Self::Cl37,
            "39k" | "k39" | "potassium" => Self::K39,
            "40k" | "k40" => Self::K40,
            "41k" | "k41" => Self::K41,
            "43ca" | "ca43" | "calcium" => Self::Ca43,
            "45sc" | "sc45" | "scandium" => Self::Sc45,
            "47ti" | "ti47" => Self::Ti47,
            "49ti" | "ti49" | "titanium" => Self::Ti49,
            "50v" | "v50" => Self::V50,
            "51v" | "v51" | "vanadium" => Self::V51,
            "53cr" | "cr53" | "chromium" => Self::Cr53,
            "55mn" | "mn55" | "manganese" => Self::Mn55,
            "57fe" | "fe57" | "iron" => Self::Fe57,
            "59co" | "co59" | "cobalt" => Self::Co59,
            "61ni" | "ni61" | "nickel" => Self::Ni61,
            "63cu" | "cu63" | "copper" => Self::Cu63,
            "65cu" | "cu65" => Self::Cu65,
            "67zn" | "zn67" | "zinc" => Self::Zn67,
            "69ga" | "ga69" | "gallium" => Self::Ga69,
            "71ga" | "ga71" => Self::Ga71,
            "73ge" | "ge73" | "germanium" => Self::Ge73,
            "75as" | "as75" | "arsenic" | "arsenico" => Self::As75,
            "77se" | "se77" | "selenium" => Self::Se77,
            "79br" | "br79" | "bromine" => Self::Br79,
            "81br" | "br81" => Self::Br81,
            "85rb" | "rb85" => Self::Rb85,
            "87rb" | "rb87" | "rubidium" => Self::Rb87,
            "87sr" | "sr87" | "strontium" => Self::Sr87,
            "91zr" | "zr91" | "zirconium" => Self::Zr91,
            "93nb" | "nb93" | "niobium" => Self::Nb93,
            "95mo" | "mo95" | "molybdenum" => Self::Mo95,
            "97mo" | "mo97" => Self::Mo97,
            "99ru" | "ru99" => Self::Ru99,
            "101ru" | "ru101" | "ruthenium" => Self::Ru101,
            "103rh" | "rh103" | "rhodium" => Self::Rh103,
            "105pd" | "pd105" | "palladium" => Self::Pd105,
            "107ag" | "ag107" | "silver" => Self::Ag107,
            "109ag" | "ag109" => Self::Ag109,
            "111cd" | "cd111" | "cadmium" => Self::Cd111,
            "113cd" | "cd113" => Self::Cd113,
            "113in" | "in113" => Self::In113,
            "115in" | "in115" | "indium" => Self::In115,
            "115sn" | "sn115" => Self::Sn115,
            "117sn" | "sn117" => Self::Sn117,
            "119sn" | "sn119" | "tin" => Self::Sn119,
            "121sb" | "sb121" | "antimony" => Self::Sb121,
            "123sb" | "sb123" => Self::Sb123,
            "123te" | "te123" => Self::Te123,
            "125te" | "te125" | "tellurium" => Self::Te125,
            "127i" | "i127" | "iodine" => Self::I127,
            "129xe" | "xe129" | "xenon" => Self::Xe129,
            "131xe" | "xe131" => Self::Xe131,
            "133cs" | "cs133" | "cesium" | "caesium" => Self::Cs133,
            "135ba" | "ba135" => Self::Ba135,
            "137ba" | "ba137" | "barium" => Self::Ba137,
            "183w" | "w183" | "tungsten" | "wolfram" => Self::W183,
            "195pt" | "pt195" | "platinum" => Self::Pt195,
            "197au" | "au197" | "gold" => Self::Au197,
            "199hg" | "hg199" | "mercury" => Self::Hg199,
            "201hg" | "hg201" => Self::Hg201,
            "203tl" | "tl203" => Self::Tl203,
            "205tl" | "tl205" | "thallium" => Self::Tl205,
            "207pb" | "pb207" | "lead" => Self::Pb207,
            "209bi" | "bi209" | "bismuth" => Self::Bi209,
            _ => {
                if let Some(index) = METAL_CATALOG.iter().position(|entry| {
                    let label = entry.label.to_ascii_lowercase();
                    let swapped = format!(
                        "{}{}",
                        &label[label.len() - entry.symbol.len()..],
                        &label[..label.len() - entry.symbol.len()]
                    );
                    normalized == label || normalized == swapped
                }) {
                    return Ok(Self::Catalog(index as u8));
                }
                return Err(format!(
                    "Unknown nucleus '{}'. Supported: {}",
                    alias,
                    Self::supported_labels().join(", ")
                ));
            }
        };

        Ok(nucleus)
    }

    /// All accepted label spellings, one per supported nuclide.
    ///
    /// The list contains the canonical labels of the 85 dedicated variants
    /// followed by the catalogue metal labels. It is the same list that the
    /// error message of [`NmrNucleus::parse`] quotes.
    pub fn supported_labels() -> Vec<&'static str> {
        let mut labels: Vec<&'static str> = Self::ALL
            .iter()
            .map(|nucleus| nucleus.canonical())
            .collect();
        labels.extend(METAL_CATALOG.iter().map(|entry| entry.label));
        labels
    }

    /// Proton number (Z) of the element carrying this nuclide.
    ///
    /// Dimensionless. For catalogue metals the value is read from the
    /// catalogue entry, so it is always the atomic number and not the mass
    /// number.
    pub fn atomic_number(self) -> u8 {
        match self {
            Self::H1 | Self::H2 | Self::H3 => 1,
            Self::He3 => 2,
            Self::Li6 | Self::Li7 => 3,
            Self::Be9 => 4,
            Self::B10 | Self::B11 => 5,
            Self::C13 => 6,
            Self::N14 | Self::N15 => 7,
            Self::O17 => 8,
            Self::F19 => 9,
            Self::Na23 => 11,
            Self::Mg25 => 12,
            Self::Al27 => 13,
            Self::Si29 => 14,
            Self::P31 => 15,
            Self::S33 => 16,
            Self::Cl35 | Self::Cl37 => 17,
            Self::K39 | Self::K40 | Self::K41 => 19,
            Self::Ca43 => 20,
            Self::Sc45 => 21,
            Self::Ti47 | Self::Ti49 => 22,
            Self::V50 | Self::V51 => 23,
            Self::Cr53 => 24,
            Self::Mn55 => 25,
            Self::Fe57 => 26,
            Self::Co59 => 27,
            Self::Ni61 => 28,
            Self::Cu63 | Self::Cu65 => 29,
            Self::Zn67 => 30,
            Self::Ga69 | Self::Ga71 => 31,
            Self::Ge73 => 32,
            Self::As75 => 33,
            Self::Se77 => 34,
            Self::Br79 | Self::Br81 => 35,
            Self::Rb85 | Self::Rb87 => 37,
            Self::Sr87 => 38,
            Self::Zr91 => 40,
            Self::Nb93 => 41,
            Self::Mo95 | Self::Mo97 => 42,
            Self::Ru99 | Self::Ru101 => 44,
            Self::Rh103 => 45,
            Self::Pd105 => 46,
            Self::Ag107 | Self::Ag109 => 47,
            Self::Cd111 | Self::Cd113 => 48,
            Self::In113 | Self::In115 => 49,
            Self::Sn115 | Self::Sn117 | Self::Sn119 => 50,
            Self::Sb121 | Self::Sb123 => 51,
            Self::Te123 | Self::Te125 => 52,
            Self::I127 => 53,
            Self::Xe129 | Self::Xe131 => 54,
            Self::Cs133 => 55,
            Self::Ba135 | Self::Ba137 => 56,
            Self::W183 => 74,
            Self::Pt195 => 78,
            Self::Au197 => 79,
            Self::Hg199 | Self::Hg201 => 80,
            Self::Tl203 | Self::Tl205 => 81,
            Self::Pb207 => 82,
            Self::Bi209 => 83,
            Self::Catalog(index) => catalog_entry(index).z,
        }
    }

    /// Mass number (A) of this nuclide.
    ///
    /// Dimensionless: the integer nucleon count, e.g. 13 for ¹³C. It is the
    /// nominal mass, not the exact isotopic mass in Da.
    pub fn mass_number(self) -> u16 {
        match self {
            Self::H1 => 1,
            Self::H2 => 2,
            Self::H3 => 3,
            Self::He3 => 3,
            Self::Li6 => 6,
            Self::Li7 => 7,
            Self::Be9 => 9,
            Self::B10 => 10,
            Self::B11 => 11,
            Self::C13 => 13,
            Self::N14 => 14,
            Self::N15 => 15,
            Self::O17 => 17,
            Self::F19 => 19,
            Self::Na23 => 23,
            Self::Mg25 => 25,
            Self::Al27 => 27,
            Self::Si29 => 29,
            Self::P31 => 31,
            Self::S33 => 33,
            Self::Cl35 => 35,
            Self::Cl37 => 37,
            Self::K39 => 39,
            Self::K40 => 40,
            Self::K41 => 41,
            Self::Ca43 => 43,
            Self::Sc45 => 45,
            Self::Ti47 => 47,
            Self::Ti49 => 49,
            Self::V50 => 50,
            Self::V51 => 51,
            Self::Cr53 => 53,
            Self::Mn55 => 55,
            Self::Fe57 => 57,
            Self::Co59 => 59,
            Self::Ni61 => 61,
            Self::Cu63 => 63,
            Self::Cu65 => 65,
            Self::Zn67 => 67,
            Self::Ga69 => 69,
            Self::Ga71 => 71,
            Self::Ge73 => 73,
            Self::As75 => 75,
            Self::Se77 => 77,
            Self::Br79 => 79,
            Self::Br81 => 81,
            Self::Rb85 => 85,
            Self::Rb87 => 87,
            Self::Sr87 => 87,
            Self::Zr91 => 91,
            Self::Nb93 => 93,
            Self::Mo95 => 95,
            Self::Mo97 => 97,
            Self::Ru99 => 99,
            Self::Ru101 => 101,
            Self::Rh103 => 103,
            Self::Pd105 => 105,
            Self::Ag107 => 107,
            Self::Ag109 => 109,
            Self::Cd111 => 111,
            Self::Cd113 => 113,
            Self::In113 => 113,
            Self::In115 => 115,
            Self::Sn115 => 115,
            Self::Sn117 => 117,
            Self::Sn119 => 119,
            Self::Sb121 => 121,
            Self::Sb123 => 123,
            Self::Te123 => 123,
            Self::Te125 => 125,
            Self::I127 => 127,
            Self::Xe129 => 129,
            Self::Xe131 => 131,
            Self::Cs133 => 133,
            Self::Ba135 => 135,
            Self::Ba137 => 137,
            Self::W183 => 183,
            Self::Pt195 => 195,
            Self::Au197 => 197,
            Self::Hg199 => 199,
            Self::Hg201 => 201,
            Self::Tl203 => 203,
            Self::Tl205 => 205,
            Self::Pb207 => 207,
            Self::Bi209 => 209,
            Self::Catalog(index) => catalog_entry(index).mass,
        }
    }

    /// Canonical ASCII label of the nuclide, e.g. `"13C"` or `"89Y"`.
    ///
    /// The mass number precedes the element symbol. This is the label used
    /// in the SMILES/SMARTS-style atom map notations and in the error message
    /// of [`NmrNucleus::parse`].
    pub fn canonical(self) -> &'static str {
        match self {
            Self::H1 => "1H",
            Self::H2 => "2H",
            Self::H3 => "3H",
            Self::He3 => "3He",
            Self::Li6 => "6Li",
            Self::Li7 => "7Li",
            Self::Be9 => "9Be",
            Self::B10 => "10B",
            Self::B11 => "11B",
            Self::C13 => "13C",
            Self::N14 => "14N",
            Self::N15 => "15N",
            Self::O17 => "17O",
            Self::F19 => "19F",
            Self::Na23 => "23Na",
            Self::Mg25 => "25Mg",
            Self::Al27 => "27Al",
            Self::Si29 => "29Si",
            Self::P31 => "31P",
            Self::S33 => "33S",
            Self::Cl35 => "35Cl",
            Self::Cl37 => "37Cl",
            Self::K39 => "39K",
            Self::K40 => "40K",
            Self::K41 => "41K",
            Self::Ca43 => "43Ca",
            Self::Sc45 => "45Sc",
            Self::Ti47 => "47Ti",
            Self::Ti49 => "49Ti",
            Self::V50 => "50V",
            Self::V51 => "51V",
            Self::Cr53 => "53Cr",
            Self::Mn55 => "55Mn",
            Self::Fe57 => "57Fe",
            Self::Co59 => "59Co",
            Self::Ni61 => "61Ni",
            Self::Cu63 => "63Cu",
            Self::Cu65 => "65Cu",
            Self::Zn67 => "67Zn",
            Self::Ga69 => "69Ga",
            Self::Ga71 => "71Ga",
            Self::Ge73 => "73Ge",
            Self::As75 => "75As",
            Self::Se77 => "77Se",
            Self::Br79 => "79Br",
            Self::Br81 => "81Br",
            Self::Rb85 => "85Rb",
            Self::Rb87 => "87Rb",
            Self::Sr87 => "87Sr",
            Self::Zr91 => "91Zr",
            Self::Nb93 => "93Nb",
            Self::Mo95 => "95Mo",
            Self::Mo97 => "97Mo",
            Self::Ru99 => "99Ru",
            Self::Ru101 => "101Ru",
            Self::Rh103 => "103Rh",
            Self::Pd105 => "105Pd",
            Self::Ag107 => "107Ag",
            Self::Ag109 => "109Ag",
            Self::Cd111 => "111Cd",
            Self::Cd113 => "113Cd",
            Self::In113 => "113In",
            Self::In115 => "115In",
            Self::Sn115 => "115Sn",
            Self::Sn117 => "117Sn",
            Self::Sn119 => "119Sn",
            Self::Sb121 => "121Sb",
            Self::Sb123 => "123Sb",
            Self::Te123 => "123Te",
            Self::Te125 => "125Te",
            Self::I127 => "127I",
            Self::Xe129 => "129Xe",
            Self::Xe131 => "131Xe",
            Self::Cs133 => "133Cs",
            Self::Ba135 => "135Ba",
            Self::Ba137 => "137Ba",
            Self::W183 => "183W",
            Self::Pt195 => "195Pt",
            Self::Au197 => "197Au",
            Self::Hg199 => "199Hg",
            Self::Hg201 => "201Hg",
            Self::Tl203 => "203Tl",
            Self::Tl205 => "205Tl",
            Self::Pb207 => "207Pb",
            Self::Bi209 => "209Bi",
            Self::Catalog(index) => catalog_entry(index).label,
        }
    }

    /// Unicode superscript label of the nuclide, e.g. `"¹³C"`.
    ///
    /// Cosmetically identical to [`NmrNucleus::canonical`] but safe to print
    /// directly into a spectrum axis label or a report. The element symbol is
    /// kept in plain ASCII; only the mass number is superscripted.
    pub fn unicode_label(self) -> &'static str {
        match self {
            Self::H1 => "¹H",
            Self::H2 => "²H",
            Self::H3 => "³H",
            Self::He3 => "³He",
            Self::Li6 => "⁶Li",
            Self::Li7 => "⁷Li",
            Self::Be9 => "⁹Be",
            Self::B10 => "¹⁰B",
            Self::B11 => "¹¹B",
            Self::C13 => "¹³C",
            Self::N14 => "¹⁴N",
            Self::N15 => "¹⁵N",
            Self::O17 => "¹⁷O",
            Self::F19 => "¹⁹F",
            Self::Na23 => "²³Na",
            Self::Mg25 => "²⁵Mg",
            Self::Al27 => "²⁷Al",
            Self::Si29 => "²⁹Si",
            Self::P31 => "³¹P",
            Self::S33 => "³³S",
            Self::Cl35 => "³⁵Cl",
            Self::Cl37 => "³⁷Cl",
            Self::K39 => "³⁹K",
            Self::K40 => "⁴⁰K",
            Self::K41 => "⁴¹K",
            Self::Ca43 => "⁴³Ca",
            Self::Sc45 => "⁴⁵Sc",
            Self::Ti47 => "⁴⁷Ti",
            Self::Ti49 => "⁴⁹Ti",
            Self::V50 => "⁵⁰V",
            Self::V51 => "⁵¹V",
            Self::Cr53 => "⁵³Cr",
            Self::Mn55 => "⁵⁵Mn",
            Self::Fe57 => "⁵⁷Fe",
            Self::Co59 => "⁵⁹Co",
            Self::Ni61 => "⁶¹Ni",
            Self::Cu63 => "⁶³Cu",
            Self::Cu65 => "⁶⁵Cu",
            Self::Zn67 => "⁶⁷Zn",
            Self::Ga69 => "⁶⁹Ga",
            Self::Ga71 => "⁷¹Ga",
            Self::Ge73 => "⁷³Ge",
            Self::As75 => "⁷⁵As",
            Self::Se77 => "⁷⁷Se",
            Self::Br79 => "⁷⁹Br",
            Self::Br81 => "⁸¹Br",
            Self::Rb85 => "⁸⁵Rb",
            Self::Rb87 => "⁸⁷Rb",
            Self::Sr87 => "⁸⁷Sr",
            Self::Zr91 => "⁹¹Zr",
            Self::Nb93 => "⁹³Nb",
            Self::Mo95 => "⁹⁵Mo",
            Self::Mo97 => "⁹⁷Mo",
            Self::Ru99 => "⁹⁹Ru",
            Self::Ru101 => "¹⁰¹Ru",
            Self::Rh103 => "¹⁰³Rh",
            Self::Pd105 => "¹⁰⁵Pd",
            Self::Ag107 => "¹⁰⁷Ag",
            Self::Ag109 => "¹⁰⁹Ag",
            Self::Cd111 => "¹¹¹Cd",
            Self::Cd113 => "¹¹³Cd",
            Self::In113 => "¹¹³In",
            Self::In115 => "¹¹⁵In",
            Self::Sn115 => "¹¹⁵Sn",
            Self::Sn117 => "¹¹⁷Sn",
            Self::Sn119 => "¹¹⁹Sn",
            Self::Sb121 => "¹²¹Sb",
            Self::Sb123 => "¹²³Sb",
            Self::Te123 => "¹²³Te",
            Self::Te125 => "¹²⁵Te",
            Self::I127 => "¹²⁷I",
            Self::Xe129 => "¹²⁹Xe",
            Self::Xe131 => "¹³¹Xe",
            Self::Cs133 => "¹³³Cs",
            Self::Ba135 => "¹³⁵Ba",
            Self::Ba137 => "¹³⁷Ba",
            Self::W183 => "¹⁸³W",
            Self::Pt195 => "¹⁹⁵Pt",
            Self::Au197 => "¹⁹⁷Au",
            Self::Hg199 => "¹⁹⁹Hg",
            Self::Hg201 => "²⁰¹Hg",
            Self::Tl203 => "²⁰³Tl",
            Self::Tl205 => "²⁰⁵Tl",
            Self::Pb207 => "²⁰⁷Pb",
            Self::Bi209 => "²⁰⁹Bi",
            Self::Catalog(index) => catalog_entry(index).pretty,
        }
    }

    /// Element symbol of the nuclide, e.g. `"C"` for ¹³C.
    ///
    /// The symbol carries no mass number. An element outside the supported
    /// range yields the placeholder `"X"`.
    pub fn element_symbol(self) -> &'static str {
        if let Self::Catalog(index) = self {
            return catalog_entry(index).symbol;
        }
        match self.atomic_number() {
            1 => "H",
            2 => "He",
            3 => "Li",
            4 => "Be",
            5 => "B",
            6 => "C",
            7 => "N",
            8 => "O",
            9 => "F",
            11 => "Na",
            12 => "Mg",
            13 => "Al",
            14 => "Si",
            15 => "P",
            16 => "S",
            17 => "Cl",
            19 => "K",
            20 => "Ca",
            21 => "Sc",
            22 => "Ti",
            23 => "V",
            24 => "Cr",
            25 => "Mn",
            26 => "Fe",
            27 => "Co",
            28 => "Ni",
            29 => "Cu",
            30 => "Zn",
            31 => "Ga",
            32 => "Ge",
            33 => "As",
            34 => "Se",
            35 => "Br",
            37 => "Rb",
            38 => "Sr",
            40 => "Zr",
            41 => "Nb",
            42 => "Mo",
            44 => "Ru",
            45 => "Rh",
            46 => "Pd",
            47 => "Ag",
            48 => "Cd",
            49 => "In",
            50 => "Sn",
            51 => "Sb",
            52 => "Te",
            53 => "I",
            54 => "Xe",
            55 => "Cs",
            56 => "Ba",
            74 => "W",
            78 => "Pt",
            79 => "Au",
            80 => "Hg",
            81 => "Tl",
            82 => "Pb",
            83 => "Bi",
            _ => "X",
        }
    }

    /// Whether this is one of the two nuclei the shift predictor is tuned for.
    ///
    /// True for ¹H and ¹³C only. Those two carry dedicated additivity rules
    /// and receive a higher prediction confidence; every other nucleus falls
    /// back to the generic relative model.
    pub fn is_primary_target(self) -> bool {
        matches!(self, Self::H1 | Self::C13)
    }

    /// Whether the nucleus is treated as quadrupolar by the shift predictor.
    ///
    /// Quadrupolar nuclei (nuclear spin `I > 1/2`) have no sharp resonance:
    /// the quadrupolar interaction broadens them, so the predictor only
    /// returns a rough relative estimate for them and flags the lower
    /// confidence accordingly.
    ///
    /// For catalogue metals the answer follows the tabulated nuclear spin
    /// (`I > 1/2`). For the dedicated variants it follows a hard-coded
    /// whitelist of nuclides that are treated as non-quadrupolar; every
    /// other variant is reported as quadrupolar.
    ///
    /// Caveat: ⁷³Ge has `I = 9/2` and is physically quadrupolar, but the
    /// whitelist classifies it as non-quadrupolar. The catalogue metals have
    /// no such exception.
    pub fn is_quadrupolar(self) -> bool {
        if let Self::Catalog(index) = self {
            return catalog_entry(index).spin_i > 0.5;
        }
        !matches!(
            self,
            Self::H1
                | Self::H3
                | Self::He3
                | Self::C13
                | Self::N15
                | Self::F19
                | Self::Si29
                | Self::P31
                | Self::Se77
                | Self::Ge73
                | Self::Sn115
                | Self::Sn117
                | Self::Sn119
                | Self::Te123
                | Self::Te125
                | Self::Xe129
                | Self::Tl203
                | Self::Tl205
                | Self::Pb207
                | Self::Pt195
                | Self::Hg199
                | Self::Rh103
                | Self::Ag107
                | Self::Ag109
                | Self::Cd111
                | Self::Cd113
        )
    }

    /// The nuclide that is normally observed for a given element.
    ///
    /// The mapping picks the most sensitive NMR-active isotope of the
    /// element (e.g. ¹³C for carbon, ¹⁵N for nitrogen, ⁷Li for lithium).
    /// Elements without any NMR-active isotope — and elements outside the
    /// supported range — yield `None`.
    ///
    /// # Example
    ///
    /// ```
    /// use sci_form::nmr::nucleus::NmrNucleus;
    ///
    /// assert_eq!(NmrNucleus::default_for_element(6), Some(NmrNucleus::C13));
    /// assert_eq!(NmrNucleus::default_for_element(2), Some(NmrNucleus::He3));
    /// assert_eq!(NmrNucleus::default_for_element(58), None); // Ce: all isotopes I = 0
    /// ```
    pub fn default_for_element(atomic_number: u8) -> Option<Self> {
        let nucleus = match atomic_number {
            1 => Self::H1,
            2 => Self::He3,
            3 => Self::Li7,
            4 => Self::Be9,
            5 => Self::B11,
            6 => Self::C13,
            7 => Self::N15,
            8 => Self::O17,
            9 => Self::F19,
            11 => Self::Na23,
            12 => Self::Mg25,
            13 => Self::Al27,
            14 => Self::Si29,
            15 => Self::P31,
            16 => Self::S33,
            17 => Self::Cl35,
            19 => Self::K39,
            20 => Self::Ca43,
            21 => Self::Sc45,
            22 => Self::Ti49,
            23 => Self::V51,
            24 => Self::Cr53,
            25 => Self::Mn55,
            26 => Self::Fe57,
            27 => Self::Co59,
            28 => Self::Ni61,
            29 => Self::Cu63,
            30 => Self::Zn67,
            31 => Self::Ga69,
            32 => Self::Ge73,
            33 => Self::As75,
            34 => Self::Se77,
            35 => Self::Br79,
            37 => Self::Rb87,
            38 => Self::Sr87,
            40 => Self::Zr91,
            41 => Self::Nb93,
            42 => Self::Mo95,
            44 => Self::Ru101,
            45 => Self::Rh103,
            46 => Self::Pd105,
            47 => Self::Ag107,
            48 => Self::Cd111,
            49 => Self::In115,
            50 => Self::Sn119,
            51 => Self::Sb121,
            52 => Self::Te125,
            53 => Self::I127,
            54 => Self::Xe129,
            55 => Self::Cs133,
            56 => Self::Ba137,
            74 => Self::W183,
            78 => Self::Pt195,
            79 => Self::Au197,
            80 => Self::Hg199,
            81 => Self::Tl205,
            82 => Self::Pb207,
            83 => Self::Bi209,
            _ => {
                return METAL_CATALOG
                    .iter()
                    .position(|entry| entry.z == atomic_number)
                    .map(|index| Self::Catalog(index as u8));
            }
        };

        Some(nucleus)
    }

    /// Reference Larmor frequency of the nuclide, in MHz.
    ///
    /// The value is the frequency the nuclide resonates at on a standard
    /// high-field spectrometer (400 MHz for ¹H, 100.6 MHz for ¹³C, …). It is
    /// used to convert between Hz and ppm: a line width of `w` Hz corresponds
    /// to `w / frequency_mhz` ppm. For catalogue metals the tabulated
    /// frequency of the entry is returned.
    pub fn default_frequency_mhz(self) -> f64 {
        match self {
            Self::H1 => 400.0,
            Self::H2 => 61.4,
            Self::H3 => 426.0,
            Self::He3 => 304.6,
            Self::Li6 => 58.9,
            Self::Li7 => 155.5,
            Self::Be9 => 56.3,
            Self::B10 => 42.9,
            Self::B11 => 128.4,
            Self::C13 => 100.6,
            Self::N14 => 28.9,
            Self::N15 => 40.6,
            Self::O17 => 54.2,
            Self::F19 => 376.5,
            Self::Na23 => 105.8,
            Self::Mg25 => 24.5,
            Self::Al27 => 104.2,
            Self::Si29 => 79.5,
            Self::P31 => 162.0,
            Self::S33 => 30.7,
            Self::Cl35 => 39.2,
            Self::Cl37 => 32.7,
            Self::K39 => 18.7,
            Self::K40 => 2.2,
            Self::K41 => 10.9,
            Self::Ca43 => 26.9,
            Self::Sc45 => 97.2,
            Self::Ti47 => 22.6,
            Self::Ti49 => 23.2,
            Self::V50 => 42.1,
            Self::V51 => 105.2,
            Self::Cr53 => 22.6,
            Self::Mn55 => 100.0,
            Self::Fe57 => 54.2,
            Self::Co59 => 95.0,
            Self::Ni61 => 35.2,
            Self::Cu63 => 106.3,
            Self::Cu65 => 113.9,
            Self::Zn67 => 25.1,
            Self::Ga69 => 122.1,
            Self::Ga71 => 152.5,
            Self::Ge73 => 14.1,
            Self::As75 => 52.2,
            Self::Se77 => 76.3,
            Self::Br79 => 100.2,
            Self::Br81 => 108.3,
            Self::Rb85 => 130.9,
            Self::Rb87 => 186.9,
            Self::Sr87 => 74.2,
            Self::Zr91 => 24.5,
            Self::Nb93 => 98.2,
            Self::Mo95 => 26.1,
            Self::Mo97 => 26.6,
            Self::Ru99 => 17.3,
            Self::Ru101 => 22.2,
            Self::Rh103 => 15.8,
            Self::Pd105 => 17.2,
            Self::Ag107 => 18.6,
            Self::Ag109 => 21.4,
            Self::Cd111 => 88.8,
            Self::Cd113 => 93.3,
            Self::In113 => 88.6,
            Self::In115 => 92.9,
            Self::Sn115 => 149.0,
            Self::Sn117 => 149.3,
            Self::Sn119 => 149.2,
            Self::Sb121 => 102.7,
            Self::Sb123 => 54.3,
            Self::Te123 => 126.0,
            Self::Te125 => 126.3,
            Self::I127 => 80.0,
            Self::Xe129 => 110.7,
            Self::Xe131 => 32.9,
            Self::Cs133 => 52.4,
            Self::Ba135 => 28.0,
            Self::Ba137 => 42.0,
            Self::W183 => 16.7,
            Self::Pt195 => 85.8,
            Self::Au197 => 8.7,
            Self::Hg199 => 74.7,
            Self::Hg201 => 22.1,
            Self::Catalog(index) => catalog_entry(index).freq_mhz,
            Self::Tl203 => 96.6,
            Self::Tl205 => 100.2,
            Self::Pb207 => 83.7,
            Self::Bi209 => 84.2,
        }
    }

    /// Natural full width at half maximum of the resonance line, in Hz.
    ///
    /// This is the unbroadened line width of the nuclide: it is a tabulated
    /// estimate that is typically one to a few Hz for spin-½ nuclei and one
    /// to two orders of magnitude larger for quadrupolar ones (which are
    /// broadened by the quadrupolar interaction). The value is used as the
    /// default Lorentzian broadening when a spectrum is simulated without an
    /// explicit line width.
    pub fn default_fwhm_hz(self) -> f64 {
        match self {
            Self::H1 => 1.0,
            Self::C13 | Self::F19 | Self::Si29 | Self::Ge73 => 5.0,
            Self::P31 | Self::N15 | Self::Se77 | Self::Te123 | Self::Te125 => 12.0,
            Self::H3 | Self::He3 | Self::Xe129 => 15.0,
            Self::Pt195 | Self::Hg199 | Self::Pb207 | Self::W183 => 25.0,
            Self::Tl203 | Self::Tl205 | Self::Cd111 | Self::Cd113 => 30.0,
            Self::Ag107 | Self::Ag109 | Self::Sn115 | Self::Sn117 | Self::Sn119 => 35.0,
            Self::H2 => 30.0,
            Self::Li6 | Self::Li7 => 15.0,
            Self::B11 => 35.0,
            Self::B10 => 70.0,
            Self::N14 => 90.0,
            Self::Na23 | Self::Mg25 | Self::Ca43 | Self::Sr87 => 80.0,
            Self::Al27 | Self::Sc45 => 120.0,
            Self::K39 | Self::K41 | Self::Rb87 | Self::Cs133 | Self::Ba135 | Self::Ba137 => 140.0,
            Self::K40 | Self::Rb85 => 220.0,
            Self::O17 | Self::S33 => 180.0,
            Self::Cl35 | Self::Cl37 => 260.0,
            Self::Br79 | Self::Br81 => 600.0,
            Self::I127 => 1500.0,
            Self::Xe131 | Self::Hg201 => 120.0,
            Self::V50 | Self::V51 | Self::Cu63 | Self::Cu65 | Self::Ga69 | Self::Ga71 => 180.0,
            Self::Ti47 | Self::Ti49 | Self::Zn67 | Self::As75 | Self::Sb121 | Self::Sb123 => 150.0,
            Self::Cr53 | Self::Mn55 | Self::Co59 | Self::Nb93 => 260.0,
            Self::Mo95 | Self::Mo97 | Self::Ru99 | Self::Ru101 | Self::In113 | Self::In115 => 200.0,
            Self::Bi209 | Self::Au197 => 320.0,
            Self::Catalog(index) => catalog_entry(index).fwhm_hz,
            _ if self.is_quadrupolar() => 120.0,
            _ => 15.0,
        }
    }

    /// Empirical chemical-shift window of the nuclide, in ppm.
    ///
    /// The returned pair is `(low, high)` and brackets the chemical shifts
    /// that are physically plausible for the nuclide, referenced to the usual
    /// reference compound (TMS for ¹H/¹³C, CFCl₃ for ¹⁹F, H₃PO₄ for ³¹P,
    /// liquid NH₃ for ¹⁴N/¹⁵N, …). The generic shift model clamps its
    /// predictions to this window.
    pub fn empirical_range(self) -> (f64, f64) {
        match self {
            Self::H1 | Self::H2 | Self::H3 => (-2.0, 14.0),
            Self::C13 => (-10.0, 230.0),
            Self::F19 => (-260.0, 100.0),
            Self::P31 => (-250.0, 250.0),
            Self::N14 | Self::N15 => (-400.0, 700.0),
            Self::B10 | Self::B11 => (-120.0, 120.0),
            Self::Si29 => (-200.0, 200.0),
            Self::Se77 => (-1000.0, 2000.0),
            Self::O17 => (-100.0, 1200.0),
            Self::S33 => (-500.0, 800.0),
            Self::Cl35 | Self::Cl37 => (-1000.0, 1500.0),
            Self::Br79 | Self::Br81 => (-2500.0, 5000.0),
            Self::I127 => (-4000.0, 9000.0),
            Self::Te123 | Self::Te125 => (-2000.0, 3000.0),
            Self::Xe129 | Self::Xe131 => (-500.0, 500.0),
            Self::Li6 | Self::Li7 | Self::Na23 => (-80.0, 80.0),
            Self::K39 | Self::K40 | Self::K41 => (-200.0, 200.0),
            Self::Be9 | Self::Mg25 | Self::Ca43 | Self::Sr87 => (-100.0, 200.0),
            Self::Rb85 | Self::Rb87 | Self::Cs133 | Self::Ba135 | Self::Ba137 => (-150.0, 250.0),
            Self::He3 => (-50.0, 50.0),
            Self::Al27 => (-50.0, 120.0),
            Self::Sc45 | Self::Ti47 | Self::Ti49 | Self::V50 | Self::V51 => (-2500.0, 2500.0),
            Self::Cr53 | Self::Mn55 | Self::Fe57 | Self::Co59 | Self::Ni61 => (-5000.0, 10000.0),
            Self::Cu63 | Self::Cu65 | Self::Zn67 | Self::Ga69 | Self::Ga71 => (-1500.0, 1500.0),
            Self::Ge73 | Self::As75 => (-1500.0, 1500.0),
            Self::Zr91 | Self::Nb93 | Self::Mo95 | Self::Mo97 => (-4000.0, 4000.0),
            Self::Ru99 | Self::Ru101 | Self::Rh103 | Self::Pd105 => (-6000.0, 6000.0),
            Self::Ag107 | Self::Ag109 => (-2000.0, 2000.0),
            Self::Cd111 | Self::Cd113 => (-1200.0, 1200.0),
            Self::In113 | Self::In115 => (-1000.0, 1500.0),
            Self::Sn115 | Self::Sn117 | Self::Sn119 => (-2500.0, 1500.0),
            Self::Sb121 | Self::Sb123 => (-3000.0, 3000.0),
            Self::W183 => (-4000.0, 2000.0),
            Self::Pt195 => (-7000.0, 2000.0),
            Self::Au197 => (-2000.0, 8000.0),
            Self::Hg199 | Self::Hg201 => (-4000.0, 4000.0),
            Self::Tl203 | Self::Tl205 => (-3000.0, 2000.0),
            Self::Pb207 => (-6000.0, 6000.0),
            Self::Bi209 => (-8000.0, 8000.0),
            Self::Catalog(index) => catalog_entry(index).range,
        }
    }

    /// Empirical centre of the chemical-shift window, in ppm.
    ///
    /// This is the chemical shift that a featureless environment would be
    /// assigned for the nuclide; the generic relative shift model starts
    /// from this value and adds the environment-dependent corrections. It
    /// uses the same reference compounds as [`NmrNucleus::empirical_range`].
    pub fn empirical_center(self) -> f64 {
        match self {
            Self::H1 | Self::H2 | Self::H3 => 2.5,
            Self::C13 => 80.0,
            Self::N14 | Self::N15 => 120.0,
            Self::O17 => 180.0,
            Self::F19 => -120.0,
            Self::P31 => 10.0,
            Self::B10 | Self::B11 => 15.0,
            Self::Si29 => -20.0,
            Self::Se77 => 120.0,
            Self::S33 => 80.0,
            Self::Cl35 | Self::Cl37 | Self::Br79 | Self::Br81 | Self::I127 => 0.0,
            Self::Te123 | Self::Te125 => 100.0,
            Self::Al27 => 35.0,
            Self::Sc45 | Self::Ti47 | Self::Ti49 | Self::V50 | Self::V51 => 0.0,
            Self::Cr53 | Self::Mn55 | Self::Fe57 | Self::Co59 | Self::Ni61 => 1500.0,
            Self::Cu63 | Self::Cu65 | Self::Zn67 | Self::Ga69 | Self::Ga71 => 0.0,
            Self::Ge73 | Self::As75 => 100.0,
            Self::Zr91 | Self::Nb93 | Self::Mo95 | Self::Mo97 => 0.0,
            Self::Ru99 | Self::Ru101 | Self::Rh103 | Self::Pd105 => -250.0,
            Self::Ag107 | Self::Ag109 => 0.0,
            Self::Cd111 | Self::Cd113 => 0.0,
            Self::In113 | Self::In115 => 250.0,
            Self::Sn115 | Self::Sn117 | Self::Sn119 => -300.0,
            Self::Sb121 | Self::Sb123 => -150.0,
            Self::W183 => -800.0,
            Self::Pt195 => -3200.0,
            Self::Au197 => 1200.0,
            Self::Hg199 | Self::Hg201 => -200.0,
            Self::Tl203 | Self::Tl205 => 150.0,
            Self::Pb207 => 500.0,
            Self::Bi209 => 1200.0,
            Self::Catalog(index) => catalog_entry(index).center,
            Self::He3
            | Self::Li6
            | Self::Li7
            | Self::Be9
            | Self::Na23
            | Self::Mg25
            | Self::K39
            | Self::K40
            | Self::K41
            | Self::Ca43
            | Self::Rb85
            | Self::Rb87
            | Self::Sr87
            | Self::Cs133
            | Self::Ba135
            | Self::Ba137
            | Self::Xe129
            | Self::Xe131 => 0.0,
        }
    }

    /// Reference shielding constant of the nuclide, in ppm.
    ///
    /// The shielding constant `σ` is the reference value for the nuclide
    /// (31.7 ppm for ¹H, 188.0 ppm for ¹³C, …). Nuclides without a tabulated
    /// value fall back to 85 % of the free-atom diamagnetic contribution
    /// ([`NmrNucleus::free_atom_diamagnetic`]). Larger values mean stronger
    /// shielding and, for a fixed reference compound, an upfield shift.
    pub fn reference_shielding(self) -> f64 {
        match self {
            Self::H1 | Self::H2 | Self::H3 => 31.7,
            Self::C13 => 188.0,
            Self::N14 | Self::N15 => 244.0,
            Self::O17 => 324.0,
            Self::F19 => 328.0,
            Self::P31 => 328.0,
            Self::B10 | Self::B11 => 105.0,
            Self::Si29 => 420.0,
            Self::Se77 => 860.0,
            Self::S33 => 1000.0,
            _ => self.free_atom_diamagnetic() * 0.85,
        }
    }

    /// Free-atom diamagnetic shielding contribution, in ppm.
    ///
    /// The tabulated diamagnetic shielding of the isolated neutral atom of
    /// the element (17.75 ppm for H, 260.7 ppm for C, 325.5 ppm for N, …).
    /// It is the ingredient from which the parametric reference shielding of
    /// the heavier elements is built. Elements outside the table fall back to
    /// `17.0 × Z` ppm.
    pub fn free_atom_diamagnetic(self) -> f64 {
        match self.atomic_number() {
            1 => 17.75,
            2 => 59.93,
            3 => 97.0,
            4 => 140.0,
            5 => 185.0,
            6 => 260.7,
            7 => 325.5,
            8 => 398.0,
            9 => 479.0,
            11 => 610.0,
            12 => 690.0,
            13 => 780.0,
            14 => 872.0,
            15 => 993.0,
            16 => 1118.0,
            17 => 1253.0,
            atomic_number => atomic_number as f64 * 17.0,
        }
    }
}
