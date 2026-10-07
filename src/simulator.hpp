#ifndef SIMULATOR
#define SIMULATOR

#include "types.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <stdexcept>
#include <vector>
#include <cassert>

constexpr Complex ZERO{0.0, 0.0};
constexpr Complex ONE{1.0, 0.0};
constexpr Complex MINUS_ONE{-1.0, 0.0};
constexpr Complex I(0, 1);
const Complex INV_SQRT_2{1.0 / sqrt((2.0)), 0};

// A matrix for our cause will be a 2x2 matrix.
using Matrix = std::array<Complex, 4>;

// ket-bra of 1
const Matrix UNITARY_ACTIVE = {ZERO, ZERO, ZERO, ONE};
// ket-bra of 0
const Matrix UNITARY_INACTIVE = {ONE, ZERO, ZERO, ZERO};

// Retrieve unitary matrix for given gate id. Parameterized gates take a
// rotation degree parameter.
static const Matrix getUnitary(GateID id, double rotation_degree) {
  switch (id) {
  case X:
    return {ZERO, ONE, ONE, ZERO};

  case Y:
    return {ZERO, -I, I, ZERO};

  case Z:
    return {ONE, ZERO, ZERO, MINUS_ONE};

  case H:
    return {INV_SQRT_2, INV_SQRT_2, INV_SQRT_2, MINUS_ONE * INV_SQRT_2};

  case S:
    return {ONE, ZERO, ZERO, I};

  case T:
    return {ONE, ZERO, ZERO, std::polar(1.0, M_PI / 4.0)};

  case AT:
    return {ONE, ZERO, ZERO, std::polar(1.0, -M_PI / 4.0)};

  case RX: {
    Complex main_diagonal(cos(0.5 * rotation_degree), 0);
    Complex anti_diagonal(0, -sin(0.5 * rotation_degree));
    return {main_diagonal, anti_diagonal, anti_diagonal, main_diagonal};
  }

  case RY: {
    Complex main_diagonal(cos(0.5 * rotation_degree), 0);
    Complex anti_diagonal(sin(0.5 * rotation_degree), 0);
    return {main_diagonal, MINUS_ONE * anti_diagonal, anti_diagonal,
            main_diagonal};
  }

  case RZ:
    return {Complex(cos(0.5 * rotation_degree), -sin(0.5 * rotation_degree)),
            ZERO, ZERO,
            Complex(cos(0.5 * rotation_degree), sin(0.5 * rotation_degree))};

  default:
    throw std::invalid_argument("Gate not supported!");
  }
}

/**
 Applies a single qubit gate to the given statevector on the given target qubit.

 @param unitary Array of length 4 representing a 2x2 gate unitary matrix.
 @param target_qubit Index of target qubit in statevector. Least significant bit
 first.
 @param number_qubits Total number of qubits in statevector.
 @param statevector Statevector to be used for gate application. Result is
 computed in place.
*/
void applyGate(const Matrix unitary, const unsigned int target_qubit,
               const unsigned int number_qubits,
               std::vector<Complex> &statevector);

/**
 Applies a controlled two qubit gate to a target qubit based on the state of a
 given control qubit.

 @param unitary Array of length 4 representing a 2x2 gate unitary matrix.
 @param control_qubit Index of control qubit in statevector. Least significant
 bit first.
 @param target_qubit Index of target qubit in statevector. Least significant bit
 first.
 @param number_qubits Total number of qubits in statevector.
 @param statevector Statevector to be used for gate application. Result is
 computed in place.
*/
void applyControlledGate(const Matrix unitary, const unsigned int control_qubit,
                         const unsigned int target_qubit,
                         const unsigned int number_qubits,
                         std::vector<Complex> &statevector);

/**
 Simulates a quantum circuit given as a list of gates using a full statevector
 representation. Applies native gates to all single and two qubit gates; uses
 decomposition for CCX gate.

 @param gates A list of QASMGate structs representing the quantum circuit.
 @param statevector Statevector for quantum simulation.
 @param number_qubits Total number of qubits in statevector.
 @param number_processed_gates Number of already simulated gates.Grows every
 iteration of the loop.
 @param number_total_gates Number of total gates in circuit. Can grow because of
 decompositions.
 @param print_progress Bool that when true prints a progress bar. True by
 default.
*/
void simulate(std::vector<QASMGate> gates, std::vector<Complex> &statevector,
              const unsigned int number_qubits,
              unsigned int &number_processed_gates,
              unsigned int &number_total_gates,
              const bool print_progress = true);

/**
 Measures a single target qubit in statevector.

 @param statevector Statevector to measure in. Collapses statevector in place.
 @param target_qubit Index of target qubit in statevector. Least significant bit
 first.
 @param number_qubits Total number of qubits in statevector.
*/
void measure(std::vector<Complex> &statevector, const unsigned int target_qubit,
             const unsigned int number_qubits);

#endif