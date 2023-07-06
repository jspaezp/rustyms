use crate::uom::num_traits::Zero;
use crate::{r, Mass, Ratio};
use std::fmt::Write;

include!("shared/formula.rs");

impl MolecularFormula {
    pub fn monoisotopic_mass(&self) -> Option<Mass> {
        let mut mass = Mass::zero();
        for (e, i, n) in self.elements {
            mass += e.mass(i)? * Ratio::new::<r>(n as f64)
        }
        Some(mass)
    }
    pub fn average_weight(&self) -> Option<Mass> {
        let mut mass = Mass::zero();
        for (e, i, n) in self.elements {
            mass += e.average_weight(i)? * Ratio::new::<r>(n as f64)
        }
        Some(mass)
    }

    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from this collections of elements
    pub fn hill_notation(&self) -> String {
        // TODO: Handle isotopes
        let mut output = String::new();
        if let Some(carbon) = self.elements.iter().find(|e| e.0 == Element::C) {
            write!(output, "C{}", carbon.1).unwrap();
            if let Some(hydrogen) = self.elements.iter().find(|e| e.0 == Element::H) {
                write!(output, "H{}", hydrogen.1).unwrap();
            }
            for element in self
                .elements
                .iter()
                .filter(|e| e.0 != Element::H && e.0 != Element::C)
            {
                write!(output, "{}{}", element.0, element.1).unwrap();
            }
        } else {
            for element in &self.elements {
                write!(output, "{}{}", element.0, element.1).unwrap();
            }
        }
        output
    }
}
