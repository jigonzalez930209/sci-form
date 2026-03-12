use crate::graph::{Atom, Bond, BondOrder, BondStereo, ChiralType, Hybridization, Molecule};
use nalgebra::Vector3;
use petgraph::graph::NodeIndex;
use std::collections::HashMap;

pub struct SmilesParser<'a> {
    input: &'a [u8],
    pos: usize,
    rings: HashMap<u32, (NodeIndex, BondOrder, BondStereo)>,
    mol: &'a mut Molecule,
}

impl<'a> SmilesParser<'a> {
    pub fn new(smiles: &'a str, mol: &'a mut Molecule) -> Self {
        SmilesParser {
            input: smiles.as_bytes(),
            pos: 0,
            rings: HashMap::new(),
            mol,
        }
    }

    pub fn parse(&mut self) -> Result<(), String> {
        self.parse_chain(None)?;
        // After parsing, add implicit hydrogens
        self.add_implicit_hydrogens();
        Ok(())
    }

    fn peek(&self) -> Option<u8> {
        self.input.get(self.pos).copied()
    }

    fn advance(&mut self) {
        self.pos += 1;
    }

    fn parse_chain(&mut self, mut prev_node: Option<NodeIndex>) -> Result<(), String> {
        let mut prev_bond_order = BondOrder::Single;
        let mut prev_bond_stereo = BondStereo::None;

        let get_bond_order = |m: &Molecule, prev: NodeIndex, curr: NodeIndex, order: BondOrder| -> BondOrder {
            if order == BondOrder::Single 
                && m.graph[prev].hybridization == Hybridization::SP2 
                && m.graph[curr].hybridization == Hybridization::SP2 
            {
                BondOrder::Aromatic
            } else {
                order
            }
        };

        while self.pos < self.input.len() {
            let c = self.peek().unwrap();
            match c {
                b'(' => {
                    self.advance();
                    self.parse_chain(prev_node)?;
                    if self.peek() == Some(b')') {
                        self.advance();
                    } else {
                        return Err("Expected ')'".to_string());
                    }
                }
                b')' => {
                    break;
                }
                b'-' => { self.advance(); prev_bond_order = BondOrder::Single; }
                b'=' => { self.advance(); prev_bond_order = BondOrder::Double; }
                b'#' => { self.advance(); prev_bond_order = BondOrder::Triple; }
                b':' => { self.advance(); prev_bond_order = BondOrder::Aromatic; }
                b'/' => { self.advance(); prev_bond_order = BondOrder::Single; prev_bond_stereo = BondStereo::E; } // Simplify E/Z for now
                b'\\' => { self.advance(); prev_bond_order = BondOrder::Single; prev_bond_stereo = BondStereo::Z; }
                b'%' | b'0'..=b'9' => {
                    let mut ring_num = 0;
                    if c == b'%' {
                        self.advance();
                        if let (Some(d1), Some(d2)) = (self.peek(), self.input.get(self.pos + 1)) {
                            if d1.is_ascii_digit() && d2.is_ascii_digit() {
                                ring_num = ((d1 - b'0') * 10 + (d2 - b'0')) as u32;
                                self.advance();
                                self.advance();
                            }
                        }
                    } else {
                        ring_num = (c - b'0') as u32;
                        self.advance();
                    }
                    if ring_num > 0 {
                        if let Some(prev) = prev_node {
                            if let Some((target, order, stereo)) = self.rings.remove(&ring_num) {
                                let mut final_order = if prev_bond_order != BondOrder::Single { prev_bond_order.clone() } else { order };
                                final_order = get_bond_order(self.mol, prev, target, final_order);
                                let final_stereo = if prev_bond_stereo != BondStereo::None { prev_bond_stereo.clone() } else { stereo };
                                self.mol.add_bond(prev, target, Bond { order: final_order, stereo: final_stereo });
                            } else {
                                self.rings.insert(ring_num, (prev, prev_bond_order.clone(), prev_bond_stereo.clone()));
                            }
                            // Reset bond state after closing/opening ring
                            prev_bond_order = BondOrder::Single;
                            prev_bond_stereo = BondStereo::None;
                        }
                    }
                }
                b'[' => {
                    self.advance();
                    let node = self.parse_bracket_atom()?;
                    if let Some(prev) = prev_node {
                        let order = get_bond_order(self.mol, prev, node, prev_bond_order.clone());
                        self.mol.add_bond(prev, node, Bond { order, stereo: prev_bond_stereo.clone() });
                    }
                    prev_node = Some(node);
                    prev_bond_order = BondOrder::Single;
                    prev_bond_stereo = BondStereo::None;
                }
                _ if c.is_ascii_alphabetic() || c == b'*' => {
                    let node = self.parse_organic_subset()?;
                    if let Some(prev) = prev_node {
                        let order = get_bond_order(self.mol, prev, node, prev_bond_order.clone());
                        self.mol.add_bond(prev, node, Bond { order, stereo: prev_bond_stereo.clone() });
                    }
                    prev_node = Some(node);
                    prev_bond_order = BondOrder::Single;
                    prev_bond_stereo = BondStereo::None;
                }
                _ => return Err(format!("Unexpected character: {}", c as char)),
            }
        }
        Ok(())
    }

    fn parse_organic_subset(&mut self) -> Result<NodeIndex, String> {
        let mut elem_str = String::new();
        let c = self.peek().unwrap();
        self.advance();
        elem_str.push(c as char);
        if let Some(nc) = self.peek() {
            if nc.is_ascii_lowercase() && (c == b'C' && (nc == b'l' || nc == b'u' || nc == b'o' || nc == b'r') || c == b'B' && nc == b'r') {
                elem_str.push(nc as char);
                self.advance();
            }
        }
        
        // Handle aromatic
        let mut aromatic = false;
        if elem_str.chars().next().unwrap().is_lowercase() {
            aromatic = true;
            elem_str = elem_str.to_uppercase();
        }

        let element = match elem_str.as_str() {
            "C" => 6, "N" => 7, "O" => 8, "P" => 15, "S" => 16, "F" => 9, "CL" | "Cl" => 17, "BR" | "Br" => 35, "I" => 53, "B" => 5, "*" => 0,
            &_ => return Err(format!("Unknown element: {}", elem_str))
        };

        let atom = Atom {
            element,
            position: Vector3::zeros(),
            charge: 0.0,
            formal_charge: 0,
            hybridization: if aromatic { Hybridization::SP2 } else { Hybridization::SP3 },
            chiral_tag: ChiralType::Unspecified,
            explicit_h: 0,
        };

        Ok(self.mol.add_atom(atom))
    }

    fn parse_bracket_atom(&mut self) -> Result<NodeIndex, String> {
        let mut isotope = 0u32;
        let mut aromatic = false;
        let mut chiral = ChiralType::Unspecified;
        let mut h_count = 0;
        let mut charge = 0;

        // Isotope
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() {
                isotope = isotope * 10 + (c - b'0') as u32;
                self.advance();
            } else {
                break;
            }
        }

        // Element
        let mut elem_str = String::new();
        if let Some(c) = self.peek() {
            if c.is_ascii_alphabetic() {
                elem_str.push(c as char);
                self.advance();
                if let Some(nc) = self.peek() {
                    if nc.is_ascii_lowercase() && c.is_ascii_uppercase() {
                        elem_str.push(nc as char);
                        self.advance();
                    }
                }
            } else if c == b'*' {
                elem_str.push('*');
                self.advance();
            }
        }
        
        if elem_str.chars().next().unwrap().is_lowercase() {
            aromatic = true;
            elem_str = elem_str.to_uppercase();
        }

        let element: u8 = match elem_str.as_str() {
            "H" => 1, "C" => 6, "N" => 7, "O" => 8, "P" => 15, "S" => 16, "F" => 9, "CL" | "Cl" => 17, "BR" | "Br" => 35, "I" => 53,
            "FE" | "Fe" => 26, "ZN" | "Zn" => 30, "CO" | "Co" => 27, "CU" | "Cu" => 29, "B" => 5,
            "AG" | "Ag" => 47, "AU" | "Au" => 79, "PT" | "Pt" => 78, "PD" | "Pd" => 46,
            "NA" | "Na" => 11, "K" => 19, "LI" | "Li" => 3, "CA" | "Ca" => 20, "MG" | "Mg" => 12,
            "*" => 0,
            _ => 6 // Default fallback
        };

        // Chirality
        if let Some(c) = self.peek() {
            if c == b'@' {
                self.advance();
                chiral = ChiralType::TetrahedralCCW;
                if self.peek() == Some(b'@') {
                    self.advance();
                    chiral = ChiralType::TetrahedralCW;
                }
            }
        }

        // Hydrogens
        if let Some(c) = self.peek() {
            if c == b'H' {
                self.advance();
                h_count = 1;
                if let Some(hc) = self.peek() {
                    if hc.is_ascii_digit() {
                        h_count = (hc - b'0') as u8;
                        self.advance();
                    }
                }
            }
        }

        // Charge
        if let Some(c) = self.peek() {
            if c == b'+' {
                self.advance();
                charge = 1;
                if let Some(cc) = self.peek() {
                    if cc.is_ascii_digit() {
                        charge = (cc - b'0') as i8;
                        self.advance();
                    } else if cc == b'+' {
                        charge = 2;
                        self.advance();
                        if self.peek() == Some(b'+') { charge = 3; self.advance(); }
                    }
                }
            } else if c == b'-' {
                self.advance();
                charge = -1;
                if let Some(cc) = self.peek() {
                    if cc.is_ascii_digit() {
                        charge = -((cc - b'0') as i8);
                        self.advance();
                    } else if cc == b'-' {
                        charge = -2;
                        self.advance();
                        if self.peek() == Some(b'-') { charge = -3; self.advance(); }
                    }
                }
            }
        }

        if self.peek() == Some(b']') {
            self.advance();
        } else {
            return Err("Expected ']'".to_string());
        }

        let atom = Atom {
            element,
            position: Vector3::zeros(),
            charge: charge as f32,
            formal_charge: charge,
            hybridization: if aromatic { Hybridization::SP2 } else { Hybridization::SP3 }, // Simplification
            chiral_tag: chiral,
            explicit_h: h_count,
        };

        let node = self.mol.add_atom(atom);
        
        Ok(node)
    }

    fn add_implicit_hydrogens(&mut self) {
        let n = self.mol.graph.node_count();
        let mut to_add = Vec::new();

        for i in 0..n {
            let ni = NodeIndex::new(i);
            let atom = &self.mol.graph[ni];
            if atom.element == 1 || atom.element == 0 { continue; } // H or dummy

            let mut current_valents = 0;
            let mut double_bonds = 0;
            let mut triple_bonds = 0;
            let mut aromatic_bonds = 0;

            for edge in self.mol.graph.edges(ni) {
                match edge.weight().order {
                    BondOrder::Single => { current_valents += 1; },
                    BondOrder::Double => { current_valents += 2; double_bonds += 1; },
                    BondOrder::Triple => { current_valents += 3; triple_bonds += 1; },
                    BondOrder::Aromatic => { current_valents += 1; aromatic_bonds += 1; },
                    _ => { current_valents += 1; },
                };
            }
            
            // Fix hybridization for non-aromatics based on exact bond types
            if aromatic_bonds == 0 {
                let actual_hybridization = if triple_bonds > 0 || double_bonds > 1 {
                    Hybridization::SP
                } else if double_bonds == 1 {
                    Hybridization::SP2
                } else {
                    // Check for conjugation: N or O adjacent to double/aromatic bond → SP2
                    // Matches RDKit: norbs==4, degree<=3, atomHasConjugatedBond → SP2
                    let elem = atom.element;
                    let is_conjugated = if (elem == 7 || elem == 8) && double_bonds == 0 {
                        // Check if any neighbor has a double or aromatic bond (conjugation)
                        self.mol.graph.neighbors(ni).any(|nb| {
                            self.mol.graph.edges(nb).any(|e| {
                                matches!(e.weight().order, BondOrder::Double | BondOrder::Aromatic)
                            })
                        })
                    } else {
                        false
                    };
                    if is_conjugated { Hybridization::SP2 } else { Hybridization::SP3 }
                };
                
                if let Some(atom_mut) = self.mol.graph.node_weight_mut(ni) {
                    atom_mut.hybridization = actual_hybridization;
                }
            } else {
                if let Some(atom_mut) = self.mol.graph.node_weight_mut(ni) {
                    atom_mut.hybridization = Hybridization::SP2;
                }
            }
            
            let atom = &self.mol.graph[ni];
            // If it's aromatic, an uncharged Carbon expects 3 bonds normally (2 from ring, 1 external or 1 H).
            let target_val = match atom.element {
                6 => 4, // C
                7 => if atom.formal_charge == 1 { 4 } else { 3 }, // N
                8 => if atom.formal_charge == 1 { 3 } else { 2 }, // O
                9 | 17 | 35 | 53 => 1, // Halogens
                16 => 2, // S
                15 => 3, // P
                _ => 0,
            };
            
            let mut h_needed = target_val - current_valents - atom.formal_charge.abs() as i32;
            
            if atom.hybridization == Hybridization::SP2 && atom.element == 6 {
                // simple aromatic C checks
                let aromatic_bonds = self.mol.graph.edges(ni).filter(|e| e.weight().order == BondOrder::Aromatic).count();
                if aromatic_bonds == 2 {
                    h_needed = 1 - (current_valents - 2); // 1 H max on a standard aromatic carbon in a ring
                } else if aromatic_bonds == 3 {
                    h_needed = 0; // Bridgehead aromatic carbon
                }
            }

            // Bare aromatic nitrogen (lowercase 'n', no bracket/explicit H) is pyridine-type:
            // lone pair NOT in pi system, only 2 bonds needed → 0 implicit H.
            // [nH] form has explicit_h > 0 and keeps its H (pyrrole-type).
            if aromatic_bonds > 0 && atom.element == 7 && atom.explicit_h == 0 {
                h_needed = 0;
            }

            // Calculate total hydrogens to add for this atom in THIS iteration
            // to maintain exact topological index order!
            let total_h_for_this_atom = atom.explicit_h as i32 + if h_needed - (atom.explicit_h as i32) > 0 { h_needed - (atom.explicit_h as i32) } else { 0 };
            
            if total_h_for_this_atom > 0 {
                to_add.push((ni, total_h_for_this_atom as u8));
            }
        }

        // Add all hydrogens sequentially
        for (parent, count) in to_add {
            for _ in 0..count {
                let h_node = self.mol.add_atom(Atom {
                    element: 1,
                    position: Vector3::zeros(),
                    charge: 0.0,
                    formal_charge: 0,
                    hybridization: Hybridization::Unknown,
                    chiral_tag: ChiralType::Unspecified,
                    explicit_h: 0,
                });
                self.mol.add_bond(parent, h_node, Bond { order: BondOrder::Single, stereo: BondStereo::None });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_benzene_aliphatic() {
        let mol = Molecule::from_smiles("C1=CC=CC=C1").unwrap();
        assert_eq!(mol.graph.node_count(), 12); // 6 C, 6 implicit H
        assert_eq!(mol.graph.edge_count(), 12); // 6 ring bonds, 6 C-H bonds
    }

    #[test]
    fn test_parse_benzene_aromatic() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.graph.node_count(), 12); // 6 C, 6 implicit H
        assert_eq!(mol.graph.edge_count(), 12); // 6 ring bonds, 6 C-H bonds
    }

    #[test]
    fn test_parse_brackets_and_branches() {
        let mol = Molecule::from_smiles("C(C)(C)C").unwrap(); // Isobutane
        assert_eq!(mol.graph.node_count(), 14); // 4 C, 10 H
    }
}
