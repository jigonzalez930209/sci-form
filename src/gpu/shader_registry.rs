//! Centralized WGSL shader registry.
//!
//! All GPU compute shaders are catalogued here with metadata for
//! validation, tier classification, and dispatch parameter guidance.

use super::memory_budget::GpuMemoryBudget;

/// GPU compute tier classification (from algorithm analysis).
///
/// Higher tiers benefit more from GPU acceleration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GpuTier {
    /// Massive speedup: O(N⁴) two-electron integrals, O(grid × N_basis) orbital grids.
    Tier1,
    /// Significant speedup: O(N²) matrix builds (overlap, Fock, Coulomb), O(voxels) marching cubes.
    Tier2,
    /// Moderate speedup: O(N²) pairwise (D4 dispersion, EEQ Coulomb, KPM Chebyshev).
    Tier3,
    /// CPU-preferred: SCF loop control, DIIS, eigensolve (latency-bound).
    Tier4,
}

impl GpuTier {
    /// Whether GPU dispatch is recommended for this tier.
    pub fn gpu_recommended(&self) -> bool {
        matches!(self, GpuTier::Tier1 | GpuTier::Tier2)
    }
}

impl std::fmt::Display for GpuTier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GpuTier::Tier1 => write!(f, "Tier 1 (massive speedup)"),
            GpuTier::Tier2 => write!(f, "Tier 2 (significant speedup)"),
            GpuTier::Tier3 => write!(f, "Tier 3 (moderate speedup)"),
            GpuTier::Tier4 => write!(f, "Tier 4 (CPU preferred)"),
        }
    }
}

/// Descriptor for a registered GPU shader.
#[derive(Debug, Clone)]
pub struct ShaderDescriptor {
    /// Human-readable name.
    pub name: &'static str,
    /// GPU tier classification.
    pub tier: GpuTier,
    /// Workgroup size [x, y, z].
    pub workgroup_size: [u32; 3],
    /// Entry point function name.
    pub entry_point: &'static str,
    /// Number of storage bindings required.
    pub storage_bindings: u32,
    /// Number of uniform bindings required.
    pub uniform_bindings: u32,
    /// Brief description of what this shader computes.
    pub description: &'static str,
}

impl ShaderDescriptor {
    /// Total bindings (storage + uniform).
    pub fn total_bindings(&self) -> u32 {
        self.storage_bindings + self.uniform_bindings
    }

    /// Pre-flight check against memory budget limits.
    pub fn check_against_budget(&self, budget: &GpuMemoryBudget) -> Result<(), String> {
        let total = self.storage_bindings + self.uniform_bindings;
        if total > budget.limits.max_storage_buffers_per_stage + 4 {
            return Err(format!(
                "Shader '{}' needs {} bindings, budget allows {}",
                self.name, total, budget.limits.max_storage_buffers_per_stage
            ));
        }
        let invocations = self.workgroup_size[0] as u64
            * self.workgroup_size[1] as u64
            * self.workgroup_size[2] as u64;
        if invocations > budget.limits.max_invocations_per_workgroup as u64 {
            return Err(format!(
                "Shader '{}' workgroup has {} invocations, max {}",
                self.name, invocations, budget.limits.max_invocations_per_workgroup
            ));
        }
        Ok(())
    }
}

// ─── Shader catalogue ──────────────────────────────────────────────

/// Vector addition — smoke-test / validation shader.
pub const SHADER_VECTOR_ADD: ShaderDescriptor = ShaderDescriptor {
    name: "vector_add",
    tier: GpuTier::Tier4,
    workgroup_size: [64, 1, 1],
    entry_point: "main",
    storage_bindings: 3, // lhs, rhs, out
    uniform_bindings: 1, // params
    description: "Element-wise vector addition (GPU smoke test)",
};

/// Orbital grid evaluation: ψ_i(r) = Σ_μ C_{μi} φ_μ(r).
pub const SHADER_ORBITAL_GRID: ShaderDescriptor = ShaderDescriptor {
    name: "orbital_grid",
    tier: GpuTier::Tier1,
    workgroup_size: [8, 8, 4],
    entry_point: "main",
    storage_bindings: 4, // basis, mo_coeffs, primitives, output
    uniform_bindings: 1, // grid_params
    description: "MO wavefunction on 3D grid (GPU Tier 1: O(grid × N_basis))",
};

/// Marching cubes isosurface extraction.
pub const SHADER_MARCHING_CUBES: ShaderDescriptor = ShaderDescriptor {
    name: "marching_cubes",
    tier: GpuTier::Tier2,
    workgroup_size: [4, 4, 4],
    entry_point: "main",
    storage_bindings: 5, // scalar_field, edge_table, tri_table, vertices, tri_count
    uniform_bindings: 1, // mc_params
    description: "Isosurface extraction via marching cubes (GPU Tier 2: O(voxels))",
};

/// ESP grid: V(r) = Σ_A Z_A/|r-R_A| - Σ_μν P_μν ∫ φ_μ(r') φ_ν(r')/|r-r'| dr'.
pub const SHADER_ESP_GRID: ShaderDescriptor = ShaderDescriptor {
    name: "esp_grid",
    tier: GpuTier::Tier1,
    workgroup_size: [8, 8, 4],
    entry_point: "main",
    storage_bindings: 4, // atoms, density, basis_info, output
    uniform_bindings: 1, // grid_params
    description: "Electrostatic potential on 3D grid (GPU Tier 1: O(grid × N²))",
};

/// D4 dispersion pairwise energy.
pub const SHADER_D4_DISPERSION: ShaderDescriptor = ShaderDescriptor {
    name: "d4_dispersion",
    tier: GpuTier::Tier3,
    workgroup_size: [16, 16, 1],
    entry_point: "main",
    storage_bindings: 3, // positions, d4_params, output_energies
    uniform_bindings: 1, // config
    description: "D4 pairwise dispersion (GPU Tier 3: O(N²))",
};

/// EEQ Coulomb matrix.
pub const SHADER_EEQ_COULOMB: ShaderDescriptor = ShaderDescriptor {
    name: "eeq_coulomb",
    tier: GpuTier::Tier3,
    workgroup_size: [16, 16, 1],
    entry_point: "main",
    storage_bindings: 3, // positions, radii, coulomb_matrix
    uniform_bindings: 1, // config
    description: "EEQ damped Coulomb matrix gamma_ij (GPU Tier 3: O(N²))",
};

/// Electron density grid: ρ(r) = Σ_μν P_μν φ_μ(r) φ_ν(r).
pub const SHADER_DENSITY_GRID: ShaderDescriptor = ShaderDescriptor {
    name: "density_grid",
    tier: GpuTier::Tier1,
    workgroup_size: [8, 8, 4],
    entry_point: "main",
    storage_bindings: 4, // basis, density_matrix, primitives, output
    uniform_bindings: 1, // grid_params
    description: "Electron density on 3D grid (GPU Tier 1: O(grid × N²))",
};

/// Two-electron repulsion integrals: (μν|λσ).
pub const SHADER_TWO_ELECTRON: ShaderDescriptor = ShaderDescriptor {
    name: "two_electron_eri",
    tier: GpuTier::Tier1,
    workgroup_size: [64, 1, 1],
    entry_point: "main",
    storage_bindings: 4, // basis, primitives, quartets, output
    uniform_bindings: 1, // params
    description: "Two-electron repulsion integrals (GPU Tier 1: O(N⁴))",
};

/// Fock matrix build: F = H + G(P).
pub const SHADER_FOCK_BUILD: ShaderDescriptor = ShaderDescriptor {
    name: "fock_build",
    tier: GpuTier::Tier1,
    workgroup_size: [16, 16, 1],
    entry_point: "main",
    storage_bindings: 4, // h_core, density, eris, output
    uniform_bindings: 1, // params
    description: "Fock matrix construction G(P) (GPU Tier 1: O(N⁴))",
};

/// One-electron matrices: S, T, V.
pub const SHADER_ONE_ELECTRON: ShaderDescriptor = ShaderDescriptor {
    name: "one_electron",
    tier: GpuTier::Tier2,
    workgroup_size: [16, 16, 1],
    entry_point: "main",
    storage_bindings: 4, // basis, primitives, atoms, output
    uniform_bindings: 1, // params
    description: "One-electron matrices S,T,V (GPU Tier 2: O(N²))",
};

/// SCC-DFTB gamma matrix.
pub const SHADER_GAMMA_MATRIX: ShaderDescriptor = ShaderDescriptor {
    name: "gamma_matrix",
    tier: GpuTier::Tier3,
    workgroup_size: [16, 16, 1],
    entry_point: "main",
    storage_bindings: 3, // eta, positions, output
    uniform_bindings: 1, // params
    description: "SCC-DFTB gamma matrix (GPU Tier 3: O(N²) pairwise Coulomb)",
};

/// ALPB Born radii.
pub const SHADER_ALPB_BORN_RADII: ShaderDescriptor = ShaderDescriptor {
    name: "alpb_born_radii",
    tier: GpuTier::Tier3,
    workgroup_size: [64, 1, 1],
    entry_point: "main",
    storage_bindings: 3, // positions, rho, output
    uniform_bindings: 1, // params
    description: "ALPB Born radii (GPU Tier 3: O(N²) descreening)",
};

/// CPM Coulomb matrix.
pub const SHADER_CPM_COULOMB: ShaderDescriptor = ShaderDescriptor {
    name: "cpm_coulomb",
    tier: GpuTier::Tier3,
    workgroup_size: [16, 16, 1],
    entry_point: "main",
    storage_bindings: 2, // positions, output
    uniform_bindings: 1, // params
    description: "CPM Coulomb matrix J_ij (GPU Tier 3: O(N²) pairwise electrostatics)",
};

/// All registered shaders.
pub const ALL_SHADERS: &[&ShaderDescriptor] = &[
    &SHADER_VECTOR_ADD,
    &SHADER_ORBITAL_GRID,
    &SHADER_MARCHING_CUBES,
    &SHADER_ESP_GRID,
    &SHADER_D4_DISPERSION,
    &SHADER_EEQ_COULOMB,
    &SHADER_DENSITY_GRID,
    &SHADER_TWO_ELECTRON,
    &SHADER_FOCK_BUILD,
    &SHADER_ONE_ELECTRON,
    &SHADER_GAMMA_MATRIX,
    &SHADER_ALPB_BORN_RADII,
    &SHADER_CPM_COULOMB,
];

/// Look up a shader by name.
pub fn find_shader(name: &str) -> Option<&'static ShaderDescriptor> {
    ALL_SHADERS.iter().find(|s| s.name == name).copied()
}

/// List all shaders in a given tier.
pub fn shaders_by_tier(tier: GpuTier) -> Vec<&'static ShaderDescriptor> {
    ALL_SHADERS
        .iter()
        .filter(|s| s.tier == tier)
        .copied()
        .collect()
}

/// Generate a summary report of all registered shaders.
pub fn shader_catalogue_report() -> String {
    let mut report = String::from("GPU Shader Catalogue\n====================\n\n");
    for tier in &[
        GpuTier::Tier1,
        GpuTier::Tier2,
        GpuTier::Tier3,
        GpuTier::Tier4,
    ] {
        let shaders = shaders_by_tier(*tier);
        if shaders.is_empty() {
            continue;
        }
        report.push_str(&format!("{tier}\n"));
        for s in &shaders {
            report.push_str(&format!(
                "  {} — wg[{},{},{}], {} bindings — {}\n",
                s.name,
                s.workgroup_size[0],
                s.workgroup_size[1],
                s.workgroup_size[2],
                s.total_bindings(),
                s.description,
            ));
        }
        report.push('\n');
    }
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tier_display() {
        assert_eq!(GpuTier::Tier1.to_string(), "Tier 1 (massive speedup)");
        assert!(GpuTier::Tier1.gpu_recommended());
        assert!(!GpuTier::Tier4.gpu_recommended());
    }

    #[test]
    fn test_shader_lookup() {
        let s = find_shader("orbital_grid").unwrap();
        assert_eq!(s.tier, GpuTier::Tier1);
        assert_eq!(s.workgroup_size, [8, 8, 4]);
    }

    #[test]
    fn test_shader_lookup_missing() {
        assert!(find_shader("nonexistent").is_none());
    }

    #[test]
    fn test_shaders_by_tier() {
        let t1 = shaders_by_tier(GpuTier::Tier1);
        assert!(t1.len() >= 3); // orbital_grid, esp_grid, density_grid
        assert!(t1.iter().all(|s| s.tier == GpuTier::Tier1));
    }

    #[test]
    fn test_budget_check_passes() {
        let budget = GpuMemoryBudget::webgpu_default();
        assert!(SHADER_ORBITAL_GRID.check_against_budget(&budget).is_ok());
        assert!(SHADER_MARCHING_CUBES.check_against_budget(&budget).is_ok());
    }

    #[test]
    fn test_catalogue_report() {
        let report = shader_catalogue_report();
        assert!(report.contains("orbital_grid"));
        assert!(report.contains("Tier 1"));
        assert!(report.contains("Tier 3"));
    }

    #[test]
    fn test_total_bindings() {
        assert_eq!(SHADER_ORBITAL_GRID.total_bindings(), 5);
        assert_eq!(SHADER_VECTOR_ADD.total_bindings(), 4);
    }

    #[test]
    fn test_all_shaders_registered() {
        assert_eq!(ALL_SHADERS.len(), 13);
    }
}
