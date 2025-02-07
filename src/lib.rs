use std::{borrow::Borrow, fmt::Display};

use thiserror::Error;

#[cfg(feature = "pyo3")]
use pyo3::{pyclass, pymethods};

#[cfg(feature = "rayon")]
use rayon::prelude::*;

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

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Gate {
    I(usize),
    X(usize),
    Y(usize),
    Z(usize),
    H(usize),
    S(usize),
    Cz(usize, usize),
    Cx(usize, usize),
    Swap(usize, usize),
    T(usize),
    Tdg(usize),
    Rx(usize, f64),
    Ry(usize, f64),
    Rz(usize, f64),
    Ccz(usize, usize, usize),
    Ccx(usize, usize, usize),
    Cswap(usize, usize, usize),
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
            Gate::Swap(a, b) => write!(f, "CX_{},{}", a, b),
            Gate::T(q) => write!(f, "T_{}", q),
            Gate::Tdg(q) => write!(f, "Tdg_{}", q),
            Gate::Rx(q, a) => write!(f, "Rx_{}({})", q, a),
            Gate::Ry(q, a) => write!(f, "Ry_{}({})", q, a),
            Gate::Rz(q, a) => write!(f, "Rz_{}({})", q, a),
            Gate::Ccz(c_1, c_2, t) => write!(f, "CCZ_{},{},{}", c_1, c_2, t),
            Gate::Ccx(c_1, c_2, t) => write!(f, "CCX_{},{},{}", c_1, c_2, t),
            Gate::Cswap(c, a, b) => write!(f, "CSWAP_{},{},{}", c, a, b),
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
            Gate::Swap(a, b) => vec![*a, *b],
            Gate::T(q) => vec![*q],
            Gate::Tdg(q) => vec![*q],
            Gate::Rx(q, _) => vec![*q],
            Gate::Ry(q, _) => vec![*q],
            Gate::Rz(q, _) => vec![*q],
            Gate::Ccz(c_1, c_2, t) => vec![*c_1, *c_2, *t],
            Gate::Ccx(c_1, c_2, t) => vec![*c_1, *c_2, *t],
            Gate::Cswap(c, a, b) => vec![*c, *a, *b],
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
    pub fn new_py(num_qubits: usize) -> QuantumCircuit {
        QuantumCircuit {
            num_qubits,
            ..Default::default()
        }
    }
}

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

    pub fn swap(&mut self, a: usize, b: usize) {
        self.apply_gate(Gate::Swap(a, b));
    }

    pub fn t(&mut self, q: usize) {
        self.apply_gate(Gate::T(q));
    }

    pub fn tdg(&mut self, q: usize) {
        self.apply_gate(Gate::Tdg(q));
    }

    pub fn rx(&mut self, q: usize, angle: f64) {
        self.apply_gate(Gate::Rx(q, angle));
    }

    pub fn ry(&mut self, q: usize, angle: f64) {
        self.apply_gate(Gate::Rx(q, angle));
    }

    pub fn rz(&mut self, q: usize, angle: f64) {
        self.apply_gate(Gate::Rx(q, angle));
    }

    pub fn ccz(&mut self, control_1: usize, control_2: usize, target: usize) {
        self.apply_gate(Gate::Ccz(control_1, control_2, target));
    }

    pub fn ccx(&mut self, control_1: usize, control_2: usize, target: usize) {
        self.apply_gate(Gate::Ccx(control_1, control_2, target));
    }

    pub fn toff(&mut self, control_1: usize, control_2: usize, target: usize) {
        self.ccx(control_1, control_2, target);
    }

    pub fn toffoli(&mut self, control_1: usize, control_2: usize, target: usize) {
        self.ccx(control_1, control_2, target);
    }

    pub fn cswap(&mut self, control: usize, a: usize, b: usize) {
        self.apply_gate(Gate::Cswap(control, a, b));
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

pub trait CliffordQuantumSimulator {
    fn x(&mut self, q: usize);
    fn y(&mut self, q: usize);
    fn z(&mut self, q: usize);
    fn s(&mut self, q: usize);
    fn h(&mut self, q: usize);
    fn cz(&mut self, control: usize, target: usize) {
        self.h(target);
        self.cx(control, target);
        self.h(target);
    }
    fn cx(&mut self, control: usize, target: usize) {
        self.h(target);
        self.cz(control, target);
        self.h(target);
    }
    fn cnot(&mut self, control: usize, target: usize) {
        self.cx(control, target);
    }
    fn swap(&mut self, a: usize, b: usize) {
        self.cx(a, b);
        self.cx(b, a);
        self.cx(a, b);
    }
}

pub trait QuantumSimulator: CliffordQuantumSimulator + Clone {
    fn run_circuit(&mut self, circ: &QuantumCircuit) {
        for gate in circ.gates.iter() {
            match gate {
                Gate::I(_) => (),
                Gate::X(q) => self.x(*q),
                Gate::Y(q) => self.y(*q),
                Gate::Z(q) => self.z(*q),
                Gate::H(q) => self.h(*q),
                Gate::S(q) => self.s(*q),
                Gate::Cz(c, t) => self.cz(*c, *t),
                Gate::Cx(c, t) => self.cx(*c, *t),
                Gate::Swap(a, b) => self.swap(*a, *b),
                Gate::T(q) => self.tdg(*q),
                Gate::Tdg(q) => self.t(*q),
                Gate::Rx(q, a) => self.rx(*q, *a),
                Gate::Ry(q, a) => self.ry(*q, *a),
                Gate::Rz(q, a) => self.rz(*q, *a),
                Gate::Ccz(c_1, c_2, t) => self.ccz(*c_1, *c_2, *t),
                Gate::Ccx(c_1, c_2, t) => self.ccx(*c_1, *c_2, *t),
                Gate::Cswap(c, a, b) => self.cswap(*c, *a, *b),
            }
        }
    }

    fn t(&mut self, q: usize) {
        self.tdg(q);
        self.s(q);
    }
    fn tdg(&mut self, q: usize) {
        self.t(q);
        self.s(q);
        self.z(q);
    }
    fn rx(&mut self, q: usize, angle: f64);
    fn ry(&mut self, q: usize, angle: f64);
    fn rz(&mut self, q: usize, angle: f64);
    fn ccz(&mut self, control_1: usize, control_2: usize, target: usize) {
        self.h(target);
        self.ccx(control_1, control_2, target);
        self.h(target);
    }
    fn ccx(&mut self, control_1: usize, control_2: usize, target: usize) {
        self.h(target);
        self.ccz(control_1, control_2, target);
        self.h(target);
    }
    fn toff(&mut self, control_1: usize, control_2: usize, target: usize) {
        self.ccx(control_1, control_2, target);
    }
    fn toffoli(&mut self, control_1: usize, control_2: usize, target: usize) {
        self.ccx(control_1, control_2, target);
    }
    fn cswap(&mut self, control: usize, a: usize, b: usize) {
        self.cx(a, b);
        self.ccx(control, b, a);
        self.cx(a, b);
    }
    // fn expectation(&mut self, ) -> f64;
    // fn measure(&mut self, q: usize);
    fn project(&mut self, q: usize) -> usize;
    fn sample(&self) -> Vec<usize> {
        let mut sim = self.clone();
        (0..self.num_qubits()).map(|i| sim.project(i)).collect()
    }
    fn num_qubits(&self) -> usize;
}

#[cfg(feature = "rayon")]
pub trait ParallelQuantumSimulator: QuantumSimulator + Sync {
    fn sample_parallel(&self, num_shots: usize) -> Vec<Vec<usize>> {
        (0..num_shots)
            .into_par_iter()
            .map(|_| self.sample())
            .collect()
    }
}

#[cfg(feature = "rayon")]
impl<T: QuantumSimulator + Sync> ParallelQuantumSimulator for T {}

#[cfg(feature = "pyo3")]
#[macro_export]
macro_rules! quantum_simulator_python {
    ($id:ty) => {
        #[pymethods]
        impl $id {
            #[pyo3(name = "x")]
            fn x_py(&mut self, q: usize) {
                self.x(q);
            }
            #[pyo3(name = "y")]
            fn y_py(&mut self, q: usize) {
                self.y(q);
            }
            #[pyo3(name = "z")]
            fn z_py(&mut self, q: usize) {
                self.z(q);
            }
            #[pyo3(name = "h")]
            fn h_py(&mut self, q: usize) {
                self.h(q);
            }
            #[pyo3(name = "s")]
            fn s_py(&mut self, q: usize) {
                self.s(q);
            }
            #[pyo3(name = "cz")]
            fn cz_py(&mut self, control: usize, target: usize) {
                self.cz(control, target);
            }
            #[pyo3(name = "cx")]
            fn cx_py(&mut self, control: usize, target: usize) {
                self.cx(control, target);
            }
            #[pyo3(name = "cnot")]
            fn cnot_py(&mut self, control: usize, target: usize) {
                self.cnot(control, target);
            }
            #[pyo3(name = "swap")]
            fn swap_py(&mut self, control: usize, target: usize) {
                self.swap(control, target);
            }
            #[pyo3(name = "t")]
            fn t_py(&mut self, q: usize) {
                self.t(q);
            }
            #[pyo3(name = "tdg")]
            fn tdg_py(&mut self, q: usize) {
                self.tdg(q);
            }
            #[pyo3(name = "rx")]
            fn rx_py(&mut self, q: usize, angle: f64) {
                self.rx(q, angle);
            }
            #[pyo3(name = "ry")]
            fn ry_py(&mut self, q: usize, angle: f64) {
                self.ry(q, angle);
            }
            #[pyo3(name = "rz")]
            fn rz_py(&mut self, q: usize, angle: f64) {
                self.rz(q, angle);
            }
            #[pyo3(name = "ccz")]
            fn ccz_py(&mut self, control_1: usize, control_2: usize, target: usize) {
                self.ccz(control_1, control_2, target);
            }
            #[pyo3(name = "ccx")]
            fn ccx_py(&mut self, control_1: usize, control_2: usize, target: usize) {
                self.ccx(control_1, control_2, target);
            }
            #[pyo3(name = "toff")]
            fn toff_py(&mut self, control_1: usize, control_2: usize, target: usize) {
                self.toff(control_1, control_2, target);
            }
            #[pyo3(name = "toffoli")]
            fn toffoli_py(&mut self, control_1: usize, control_2: usize, target: usize) {
                self.toffoli(control_1, control_2, target);
            }
            #[pyo3(name = "cswap")]
            fn cswap_py(&mut self, control: usize, a: usize, b: usize) {
                self.cswap(control, a, b);
            }

            #[pyo3(name = "project")]
            fn project_py(&mut self, q: usize) -> usize {
                self.project(q)
            }

            #[pyo3(name = "sample")]
            fn sample_py(&self) -> Vec<usize> {
                self.sample()
            }
            #[pyo3(name = "sample_parallel")]
            // TODO: why does this not work??
            // #[cfg(feature = "rayon")]
            fn sample_parallel_py(&self, num_shots: usize) -> Vec<Vec<usize>> {
                quantypes::ParallelQuantumSimulator::sample_parallel(self, num_shots)
            }
        }
    };
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
