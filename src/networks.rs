//! Rule-based reaction-network discovery and kMC integration.
//!
//! The discovery engine applies user-supplied rules to a set of canonical
//! SMILES, deduplicates species and records elementary edges. It is deliberately
//! rule-based: a rule owns the chemistry, while the engine handles graph
//! closure, canonicalisation, depth limits and conversion to the existing
//! Gillespie kMC solver.
//!
//! # Units
//!
//! Barriers are in eV, prefactors in s⁻¹ and counts are integer populations.
//! The kMC time is in seconds.

use std::collections::{HashMap, HashSet};

use serde::{Deserialize, Serialize};

use crate::error::SciFormError;
use crate::graph::Molecule;
use crate::kinetics::{ElementaryStep, ThermodynamicState};
use crate::kmc::{KmcConfig, KmcSpecies, KmcStep};

/// A discovery rule.
pub trait ReactionRule {
    /// Stable rule name.
    fn name(&self) -> &str;
    /// Number of reactant species consumed by the rule.
    fn reactant_count(&self) -> usize;
    /// Activation barrier in eV.
    fn barrier_ev(&self) -> f64;
    /// Optional prefactor in s⁻¹.
    fn prefactor_s_inv(&self) -> Option<f64>;

    /// Apply the rule to canonical reactant SMILES. Returns one or more
    /// canonical product sets.
    fn apply(&self, reactants: &[&str]) -> Vec<Vec<String>>;
}

/// A rule with explicit reactant and product templates.
///
/// This is the simplest useful rule: if the canonical reactant set matches,
/// return the canonical product set. It is intentionally explicit so a caller
/// can encode a known mechanism or a curated SMARTS-to-product mapping.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TemplateRule {
    /// Rule name.
    pub name: String,
    /// Reactant canonical SMILES.
    pub reactants: Vec<String>,
    /// Product canonical SMILES.
    pub products: Vec<String>,
    /// Activation barrier in eV.
    pub barrier_ev: f64,
    /// Optional prefactor in s⁻¹.
    pub prefactor_s_inv: Option<f64>,
}

impl TemplateRule {
    /// Build and canonicalize a template rule.
    pub fn new(
        name: &str,
        reactants: &[&str],
        products: &[&str],
        barrier_ev: f64,
        prefactor_s_inv: Option<f64>,
    ) -> Result<Self, SciFormError> {
        if name.trim().is_empty() {
            return Err(SciFormError::validation("rule name must not be empty"));
        }
        if reactants.is_empty() || products.is_empty() {
            return Err(SciFormError::validation(
                "a rule needs at least one reactant and one product",
            ));
        }
        if !barrier_ev.is_finite() || barrier_ev < 0.0 {
            return Err(SciFormError::validation(
                "barrier_ev must be finite and >= 0",
            ));
        }
        if let Some(prefactor) = prefactor_s_inv {
            if !prefactor.is_finite() || prefactor <= 0.0 {
                return Err(SciFormError::validation(
                    "prefactor_s_inv must be finite and > 0",
                ));
            }
        }
        let reactants = reactants
            .iter()
            .map(|smiles| canonical(smiles))
            .collect::<Result<Vec<_>, _>>()?;
        let products = products
            .iter()
            .map(|smiles| canonical(smiles))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            name: name.to_string(),
            reactants,
            products,
            barrier_ev,
            prefactor_s_inv,
        })
    }
}

impl ReactionRule for TemplateRule {
    fn name(&self) -> &str {
        &self.name
    }

    fn reactant_count(&self) -> usize {
        self.reactants.len()
    }

    fn barrier_ev(&self) -> f64 {
        self.barrier_ev
    }

    fn prefactor_s_inv(&self) -> Option<f64> {
        self.prefactor_s_inv
    }

    fn apply(&self, reactants: &[&str]) -> Vec<Vec<String>> {
        if reactants.len() != self.reactants.len() {
            return Vec::new();
        }
        let mut expected = self.reactants.clone();
        let mut actual: Vec<String> = reactants.iter().map(|value| value.to_string()).collect();
        expected.sort();
        actual.sort();
        if expected == actual {
            vec![self.products.clone()]
        } else {
            Vec::new()
        }
    }
}

/// One discovered species.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkSpecies {
    /// Stable 0-based index.
    pub index: usize,
    /// Canonical SMILES.
    pub smiles: String,
}

/// One directed elementary edge.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkEdge {
    /// Reactant species indices.
    pub reactants: Vec<usize>,
    /// Product species indices.
    pub products: Vec<usize>,
    /// Rule that generated the edge.
    pub rule: String,
    /// Activation barrier in eV.
    pub barrier_ev: f64,
    /// Optional prefactor in s⁻¹.
    pub prefactor_s_inv: Option<f64>,
}

/// Rule-based reaction network.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionNetwork {
    /// Discovered species.
    pub species: Vec<NetworkSpecies>,
    /// Directed elementary edges.
    pub edges: Vec<NetworkEdge>,
    /// Discovery notes and skipped products.
    pub notes: Vec<String>,
}

/// Discovery configuration.
pub struct NetworkDiscoveryConfig {
    /// Initial species SMILES.
    pub initial_species: Vec<String>,
    /// Rules to apply.
    pub rules: Vec<Box<dyn ReactionRule>>,
    /// Maximum BFS depth.
    pub max_depth: usize,
    /// Maximum number of species.
    pub max_species: usize,
}

/// Serialisable discovery configuration using template rules.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TemplateNetworkConfig {
    /// Initial species SMILES.
    pub initial_species: Vec<String>,
    /// Explicit template rules.
    pub rules: Vec<TemplateRule>,
    /// Maximum BFS depth.
    pub max_depth: usize,
    /// Maximum number of species.
    pub max_species: usize,
}

impl TemplateNetworkConfig {
    /// Discover the network from the serialisable template rules.
    pub fn discover(&self) -> Result<ReactionNetwork, SciFormError> {
        let rules: Vec<Box<dyn ReactionRule>> = self
            .rules
            .iter()
            .cloned()
            .map(|rule| Box::new(rule) as Box<dyn ReactionRule>)
            .collect();
        discover_reaction_network(&NetworkDiscoveryConfig {
            initial_species: self.initial_species.clone(),
            rules,
            max_depth: self.max_depth,
            max_species: self.max_species,
        })
    }
}

/// Discover a reaction network by BFS over the supplied rules.
pub fn discover_reaction_network(
    config: &NetworkDiscoveryConfig,
) -> Result<ReactionNetwork, SciFormError> {
    if config.initial_species.is_empty() {
        return Err(SciFormError::validation(
            "at least one initial species is required",
        ));
    }
    if config.rules.is_empty() {
        return Err(SciFormError::validation("at least one rule is required"));
    }
    if config.max_depth == 0 || config.max_species == 0 {
        return Err(SciFormError::validation(
            "max_depth and max_species must be >= 1",
        ));
    }

    let mut species: Vec<NetworkSpecies> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut edges: Vec<NetworkEdge> = Vec::new();
    let mut edge_keys: HashSet<(Vec<usize>, Vec<usize>, String)> = HashSet::new();
    let mut notes = Vec::new();

    for smiles in &config.initial_species {
        let canonical = canonical(smiles)?;
        if !index.contains_key(&canonical) {
            let id = species.len();
            index.insert(canonical.clone(), id);
            species.push(NetworkSpecies {
                index: id,
                smiles: canonical,
            });
        }
    }

    for _depth in 0..config.max_depth {
        let frontier = species.len();
        let mut new_species = false;
        for i in 0..frontier {
            for rule in &config.rules {
                if rule.reactant_count() == 1 {
                    let products = rule.apply(&[&species[i].smiles]);
                    for product_set in products {
                        new_species |= add_edge(
                            &product_set,
                            vec![i],
                            rule.as_ref(),
                            &mut species,
                            &mut index,
                            &mut edges,
                            &mut edge_keys,
                            &mut notes,
                            config.max_species,
                        )?;
                    }
                } else if rule.reactant_count() == 2 {
                    for j in i..frontier {
                        let products = rule.apply(&[&species[i].smiles, &species[j].smiles]);
                        for product_set in products {
                            new_species |= add_edge(
                                &product_set,
                                vec![i, j],
                                rule.as_ref(),
                                &mut species,
                                &mut index,
                                &mut edges,
                                &mut edge_keys,
                                &mut notes,
                                config.max_species,
                            )?;
                        }
                    }
                }
            }
        }
        if !new_species {
            break;
        }
    }

    Ok(ReactionNetwork {
        species,
        edges,
        notes,
    })
}

fn add_edge(
    products: &[String],
    reactants: Vec<usize>,
    rule: &dyn ReactionRule,
    species: &mut Vec<NetworkSpecies>,
    index: &mut HashMap<String, usize>,
    edges: &mut Vec<NetworkEdge>,
    edge_keys: &mut HashSet<(Vec<usize>, Vec<usize>, String)>,
    notes: &mut Vec<String>,
    max_species: usize,
) -> Result<bool, SciFormError> {
    let mut product_indices = Vec::with_capacity(products.len());
    let mut added = false;
    for product in products {
        let canonical = canonical(product)?;
        let id = if let Some(id) = index.get(&canonical) {
            *id
        } else if species.len() >= max_species {
            notes.push(format!(
                "species limit reached; skipped product {canonical}"
            ));
            continue;
        } else {
            let id = species.len();
            index.insert(canonical.clone(), id);
            species.push(NetworkSpecies {
                index: id,
                smiles: canonical,
            });
            added = true;
            id
        };
        product_indices.push(id);
    }
    if product_indices.is_empty() {
        return Ok(added);
    }
    let key = (
        reactants.clone(),
        product_indices.clone(),
        rule.name().to_string(),
    );
    if edge_keys.insert(key) {
        edges.push(NetworkEdge {
            reactants,
            products: product_indices,
            rule: rule.name().to_string(),
            barrier_ev: rule.barrier_ev(),
            prefactor_s_inv: rule.prefactor_s_inv(),
        });
    }
    Ok(added)
}

impl ReactionNetwork {
    /// Convert the network into the existing Gillespie kMC configuration.
    ///
    /// `initial_counts` must have one value per species. Every edge becomes an
    /// irreversible step with the rule's barrier and prefactor. The reaction
    /// free energy is set to zero because the discovery layer has no
    /// thermodynamic data.
    pub fn to_kmc_config(
        &self,
        initial_counts: &[u64],
        state: ThermodynamicState,
        max_events: usize,
        max_time_s: f64,
        seed: u64,
    ) -> Result<KmcConfig, SciFormError> {
        if initial_counts.len() != self.species.len() {
            return Err(SciFormError::validation(
                "initial_counts length must match the species count",
            ));
        }
        let species = self
            .species
            .iter()
            .zip(initial_counts)
            .map(|(species, &count)| KmcSpecies {
                id: species.smiles.clone(),
                initial_count: count,
            })
            .collect();
        let steps = self
            .edges
            .iter()
            .enumerate()
            .map(|(index, edge)| {
                let reactants = edge.reactants.iter().map(|&id| (id, 1)).collect::<Vec<_>>();
                let products = edge.products.iter().map(|&id| (id, 1)).collect::<Vec<_>>();
                KmcStep {
                    step: ElementaryStep {
                        step_id: format!("edge-{index}-{}", edge.rule),
                        activation_free_energy_ev: edge.barrier_ev,
                        reaction_free_energy_ev: 0.0,
                        prefactor_s_inv: edge.prefactor_s_inv,
                    },
                    reactants,
                    products,
                }
            })
            .collect();
        Ok(KmcConfig {
            species,
            steps,
            state,
            max_events,
            max_time_s,
            seed,
            frame_stride: 100,
        })
    }
}

fn canonical(smiles: &str) -> Result<String, SciFormError> {
    let molecule = Molecule::from_smiles(smiles)?;
    Ok(crate::canonical_smiles::to_canonical_smiles(&molecule))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kinetics::ThermodynamicState;
    use crate::kmc::run_kinetic_monte_carlo;

    fn linear_rules() -> Vec<Box<dyn ReactionRule>> {
        vec![
            Box::new(TemplateRule::new("C1", &["C"], &["CC"], 0.05, Some(1.0e12)).unwrap()),
            Box::new(TemplateRule::new("C2", &["CC"], &["CCC"], 0.05, Some(1.0e12)).unwrap()),
        ]
    }

    #[test]
    fn discovers_a_linear_network() {
        let config = NetworkDiscoveryConfig {
            initial_species: vec!["C".into()],
            rules: linear_rules(),
            max_depth: 3,
            max_species: 10,
        };
        let network = discover_reaction_network(&config).unwrap();
        assert_eq!(network.species.len(), 3);
        assert_eq!(network.edges.len(), 2);
        assert!(network.edges.iter().all(|edge| edge.barrier_ev == 0.05));
    }

    #[test]
    fn canonicalisation_deduplicates_species() {
        let rule = TemplateRule::new("iso", &["CCO"], &["OCC"], 0.0, None).unwrap();
        assert_eq!(rule.reactants, rule.products);
    }

    #[test]
    fn network_converts_to_kmc_and_conserves_population() {
        let config = NetworkDiscoveryConfig {
            initial_species: vec!["C".into()],
            rules: linear_rules(),
            max_depth: 3,
            max_species: 10,
        };
        let network = discover_reaction_network(&config).unwrap();
        let kmc = network
            .to_kmc_config(
                &[100, 0, 0],
                ThermodynamicState {
                    temperature_k: 500.0,
                    pressure_bar: 1.0,
                },
                1_000,
                10.0,
                7,
            )
            .unwrap();
        let result = run_kinetic_monte_carlo(&kmc).unwrap();
        assert!(!result.events.is_empty());
        assert_eq!(result.final_counts.iter().sum::<u64>(), 100);
    }

    #[test]
    fn invalid_rule_product_is_rejected() {
        assert!(TemplateRule::new("bad", &["C"], &["Q"], 0.0, None).is_err());
    }
}
