use std::{borrow::Borrow, fmt::Display};

use thiserror::Error;

#[cfg(feature = "pyo3")]
use pyo3::{pyclass, pymethods};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Pauli {
    I,
    X,
    Y,
    Z,
}

#[derive(Copy, Clone, Debug, Error, PartialEq)]
pub enum PauliError {
    #[error("provided char is not one of \"IXYZ\"")]
    NotPauli,
}

impl TryFrom<char> for Pauli {
    type Error = PauliError;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'I' | 'i' => Ok(Pauli::I),
            'X' | 'x' => Ok(Pauli::X),
            'Y' | 'y' => Ok(Pauli::Y),
            'Z' | 'z' => Ok(Pauli::Z),
            _ => Err(PauliError::NotPauli),
        }
    }
}

impl Pauli {
    /// Try to convert a stringy type into a vector of Paulis
    pub fn try_from_str(value: impl Borrow<str>) -> Result<Vec<Pauli>, PauliError> {
        // is it possible to implement this as TryFrom?
        value.borrow().chars().map(|c| c.try_into()).collect()
    }
}

impl Display for Pauli {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Pauli::I => write!(f, "I"),
            Pauli::X => write!(f, "X"),
            Pauli::Y => write!(f, "Y"),
            Pauli::Z => write!(f, "Z"),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Gate {
    I(usize),
    X(usize),
    Y(usize),
    Z(usize),
    H(usize),
    S(usize),
    Cz(usize, usize),
    Cx(usize, usize),
    T(usize),
    Ccz(usize, usize, usize),
}

impl Display for Gate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Gate::I(q) => write!(f, "I_{}", q),
            Gate::X(q) => write!(f, "X_{}", q),
            Gate::Y(q) => write!(f, "Y_{}", q),
            Gate::Z(q) => write!(f, "Z_{}", q),
            Gate::H(q) => write!(f, "H_{}", q),
            Gate::S(q) => write!(f, "S_{}", q),
            Gate::Cz(c, t) => write!(f, "CZ_{},{}", c, t),
            Gate::Cx(c, t) => write!(f, "CX_{},{}", c, t),
            Gate::T(q) => write!(f, "T_{}", q),
            Gate::Ccz(c_1, c_2, t) => write!(f, "CCZ_{},{},{}", c_1, c_2, t),
        }
    }
}

impl Gate {
    fn get_qubits(&self) -> Vec<usize> {
        match self {
            Gate::I(q) => vec![*q],
            Gate::X(q) => vec![*q],
            Gate::Y(q) => vec![*q],
            Gate::Z(q) => vec![*q],
            Gate::H(q) => vec![*q],
            Gate::S(q) => vec![*q],
            Gate::Cz(c, t) => vec![*c, *t],
            Gate::Cx(c, t) => vec![*c, *t],
            Gate::T(q) => vec![*q],
            Gate::Ccz(c_1, c_2, t) => vec![*c_1, *c_2, *t],
        }
    }
}

#[derive(Default)]
#[cfg_attr(feature = "pyo3", pyclass)]
pub struct QuantumCircuit {
    gates: Vec<Gate>,
    num_qubits: usize,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl QuantumCircuit {
    #[new]
    pub fn new(num_qubits: usize) -> QuantumCircuit {
        QuantumCircuit {
            num_qubits,
            ..Default::default()
        }
    }
}

#[cfg(not(feature = "pyo3"))]
impl QuantumCircuit {
    pub fn new(num_qubits: usize) -> QuantumCircuit {
        QuantumCircuit {
            num_qubits,
            ..Default::default()
        }
    }
}

#[cfg_attr(feature = "pyo3", pymethods)]
impl QuantumCircuit {
    pub fn x(&mut self, q: usize) {
        self.apply_gate(Gate::X(q));
    }

    pub fn y(&mut self, q: usize) {
        self.apply_gate(Gate::Y(q));
    }

    pub fn z(&mut self, q: usize) {
        self.apply_gate(Gate::Z(q));
    }

    pub fn h(&mut self, q: usize) {
        self.apply_gate(Gate::H(q));
    }

    pub fn s(&mut self, q: usize) {
        self.apply_gate(Gate::S(q));
    }

    pub fn cz(&mut self, control: usize, target: usize) {
        self.apply_gate(Gate::Cz(control, target));
    }

    pub fn cx(&mut self, control: usize, target: usize) {
        self.apply_gate(Gate::Cx(control, target));
    }

    pub fn cnot(&mut self, control: usize, target: usize) {
        self.cx(control, target);
    }

    pub fn t(&mut self, q: usize) {
        self.apply_gate(Gate::T(q));
    }

    pub fn ccz(&mut self, control_1: usize, control_2: usize, target: usize) {
        self.apply_gate(Gate::Ccz(control_1, control_2, target));
    }
}

impl QuantumCircuit {
    pub fn apply_gate(&mut self, gate: Gate) {
        for q in gate.get_qubits() {
            assert!(q < self.num_qubits);
        }
        self.gates.push(gate);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn string_to_pauli() {
        assert_eq!(
            Pauli::try_from_str("XYZ").unwrap(),
            vec![Pauli::X, Pauli::Y, Pauli::Z]
        );
        assert_eq!(Pauli::try_from_str("ABC"), Err(PauliError::NotPauli));
    }
}
