#ifndef SIMULATOR
#define SIMULATOR

#include "qasm_parser.hpp"
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

constexpr Complex ZERO{0.0, 0.0};
constexpr Complex ONE{1.0, 0.0};
constexpr Complex MINUS_ONE{-1.0, 0.0};
constexpr Complex I(0, 1);
const Complex INV_SQRT_2{1.0 / sqrt((2.0)), 0};

using Matrix = std::array<Complex, 4>;
const Matrix UNITARY_ACTIVE = {ZERO, ZERO, ZERO, ONE};
const Matrix UNITARY_INACTIVE = {ONE, ZERO, ZERO, ZERO};

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

void applyGate(const Matrix unitary, const unsigned int target_qubit,
               const unsigned int number_qubits,
               std::vector<Complex> &statevector);

void applyControlledGate(const Matrix unitary, const unsigned int control_qubit,
                         const unsigned int target_qubit,
                         const unsigned int number_qubits,
                         std::vector<Complex> &statevector);

void simulate(std::vector<QASMGate> gates, std::vector<Complex> &statevector,
              const unsigned int number_qubits, unsigned int &processed_gates,
              unsigned int &number_total_gates,
              const bool print_progress = true);

void measure(std::vector<Complex> &statevector, const unsigned int target_qubit,
             const unsigned int number_qubits);

#endif