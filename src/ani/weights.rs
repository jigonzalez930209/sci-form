//! Binary weight file loader for pre-trained ANI models.
//!
//! Simple binary format:
//! - u32: number of element types
//! - For each element type:
//!   - u8: atomic number
//!   - u32: number of layers
//!   - For each layer:
//!     - u32: rows (output_dim)
//!     - u32: cols (input_dim)
//!     - f64 × (rows × cols): weights (row-major)
//!     - f64 × rows: bias
//!     - u8: activation (0=None, 1=Gelu, 2=Celu)

use super::nn::{Activation, DenseLayer, FeedForwardNet};
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;
use std::io::{Read, Write};

/// Load ANI model weights from a binary reader.
pub fn load_weights<R: Read>(reader: &mut R) -> Result<HashMap<u8, FeedForwardNet>, String> {
    let mut buf4 = [0u8; 4];
    let mut buf1 = [0u8; 1];
    let mut buf8 = [0u8; 8];

    reader.read_exact(&mut buf4).map_err(|e| e.to_string())?;
    let n_elements = u32::from_le_bytes(buf4) as usize;

    let mut models = HashMap::new();

    for _ in 0..n_elements {
        reader.read_exact(&mut buf1).map_err(|e| e.to_string())?;
        let element = buf1[0];

        reader.read_exact(&mut buf4).map_err(|e| e.to_string())?;
        let n_layers = u32::from_le_bytes(buf4) as usize;

        let mut layers = Vec::with_capacity(n_layers);

        for _ in 0..n_layers {
            reader.read_exact(&mut buf4).map_err(|e| e.to_string())?;
            let rows = u32::from_le_bytes(buf4) as usize;

            reader.read_exact(&mut buf4).map_err(|e| e.to_string())?;
            let cols = u32::from_le_bytes(buf4) as usize;

            let mut weights = vec![0.0f64; rows * cols];
            for w in &mut weights {
                reader.read_exact(&mut buf8).map_err(|e| e.to_string())?;
                *w = f64::from_le_bytes(buf8);
            }

            let mut bias = vec![0.0f64; rows];
            for b in &mut bias {
                reader.read_exact(&mut buf8).map_err(|e| e.to_string())?;
                *b = f64::from_le_bytes(buf8);
            }

            reader.read_exact(&mut buf1).map_err(|e| e.to_string())?;
            let activation = match buf1[0] {
                1 => Activation::Gelu,
                2 => Activation::Celu,
                _ => Activation::None,
            };

            layers.push(DenseLayer {
                weights: DMatrix::from_row_slice(rows, cols, &weights),
                bias: DVector::from_vec(bias),
                activation,
            });
        }

        models.insert(element, FeedForwardNet::new(layers));
    }

    Ok(models)
}

/// Save ANI model weights to a binary writer.
pub fn save_weights<W: Write>(
    writer: &mut W,
    models: &HashMap<u8, FeedForwardNet>,
) -> Result<(), String> {
    writer
        .write_all(&(models.len() as u32).to_le_bytes())
        .map_err(|e| e.to_string())?;

    for (&element, net) in models {
        writer.write_all(&[element]).map_err(|e| e.to_string())?;
        writer
            .write_all(&(net.layers.len() as u32).to_le_bytes())
            .map_err(|e| e.to_string())?;

        for layer in &net.layers {
            let rows = layer.weights.nrows();
            let cols = layer.weights.ncols();
            writer
                .write_all(&(rows as u32).to_le_bytes())
                .map_err(|e| e.to_string())?;
            writer
                .write_all(&(cols as u32).to_le_bytes())
                .map_err(|e| e.to_string())?;

            // Weights row-major
            for r in 0..rows {
                for c in 0..cols {
                    writer
                        .write_all(&layer.weights[(r, c)].to_le_bytes())
                        .map_err(|e| e.to_string())?;
                }
            }
            // Bias
            for b in layer.bias.iter() {
                writer
                    .write_all(&b.to_le_bytes())
                    .map_err(|e| e.to_string())?;
            }

            let act_byte: u8 = match layer.activation {
                Activation::Gelu => 1,
                Activation::Celu => 2,
                Activation::None => 0,
            };
            writer.write_all(&[act_byte]).map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

/// Create a tiny test model for one element (for unit tests).
pub fn make_test_model(input_dim: usize) -> FeedForwardNet {
    let l1 = DenseLayer {
        weights: DMatrix::from_fn(16, input_dim, |r, c| {
            ((r * input_dim + c) as f64 * 0.01).sin() * 0.1
        }),
        bias: DVector::from_element(16, 0.01),
        activation: Activation::Gelu,
    };
    let l2 = DenseLayer {
        weights: DMatrix::from_fn(1, 16, |_, c| (c as f64 * 0.1).cos() * 0.05),
        bias: DVector::from_element(1, 0.0),
        activation: Activation::None,
    };
    FeedForwardNet::new(vec![l1, l2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_roundtrip() {
        let mut models = HashMap::new();
        models.insert(1u8, make_test_model(8));
        models.insert(6u8, make_test_model(8));

        let mut buf = Vec::new();
        save_weights(&mut buf, &models).unwrap();

        let mut cursor = Cursor::new(buf);
        let loaded = load_weights(&mut cursor).unwrap();

        assert_eq!(loaded.len(), 2);
        assert!(loaded.contains_key(&1));
        assert!(loaded.contains_key(&6));

        // Check a forward pass matches
        let input = DVector::from_element(8, 0.5);
        let orig = models[&1].forward(&input);
        let load = loaded[&1].forward(&input);
        assert!((orig - load).abs() < 1e-12);
    }
}
