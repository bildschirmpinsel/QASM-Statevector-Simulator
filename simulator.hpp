#ifndef SIMULATOR
#define SIMULATOR

#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using Complex = std::complex<double>;

enum GateID { X, Y, Z, H, S, T, AT, RX, RY, RZ, CX, CZ, CCX };

static const std::unordered_map<std::string, GateID> gateMap{
    {"x", GateID::X},   {"y", GateID::Y},   {"z", GateID::Z},
    {"h", GateID::H},   {"s", GateID::S},   {"t", GateID::T},
    {"rx", GateID::RX}, {"ry", GateID::RY}, {"rz", GateID::RZ},
    {"cx", GateID::CX}, {"cz", GateID::CZ}, {"ccx", GateID::CCX}};

typedef struct {
  double rotation_degree;
  int qubit_1_global_index;
  int qubit_2_global_index;
  int qubit_3_global_index;
  GateID gate;
} QASMGate;

constexpr Complex ZERO{0.0, 0.0};
constexpr Complex ONE{1.0, 0.0};
constexpr Complex MINUS_ONE{-1.0, 0.0};
constexpr Complex I(0, 1);
const Complex INV_SQRT_2{1.0 / sqrt((2.0)), 0};

using Matrix = std::array<Complex, 4>;
const Matrix UNITARY_ACTIVE = {ZERO, ZERO, ZERO, ONE};
const Matrix UNITARY_INACTIVE = {ONE, ZERO, ZERO, ZERO};

static int number_qubits = 0;

static const Matrix getUnitary(GateID id, double *rotation_degree) {
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
    Complex main_diagonal(cos(0.5 * *rotation_degree), 0);
    Complex anti_diagonal(0, -sin(0.5 * *rotation_degree));
    return {main_diagonal, anti_diagonal, anti_diagonal, main_diagonal};
  }

  case RY: {
    Complex main_diagonal(cos(0.5 * *rotation_degree), 0);
    Complex anti_diagonal(sin(0.5 * *rotation_degree), 0);
    return {main_diagonal, MINUS_ONE * anti_diagonal, anti_diagonal,
            main_diagonal};
  }

  case RZ:
    return {Complex(cos(0.5 * *rotation_degree), -sin(0.5 * *rotation_degree)),
            ZERO, ZERO,
            Complex(cos(0.5 * *rotation_degree), sin(0.5 * *rotation_degree))};

  default:
    throw std::invalid_argument("Gate not supported!");
  }
}

void applyGate(const Matrix unitary, int target_qubit,
               std::vector<Complex> &statevector);

std::vector<QASMGate>
parseQASM(const char *file, std::vector<Complex> &statevector);

void simulate(std::vector<QASMGate> gates, std::vector<Complex>& statevector);

void measure(std::vector<Complex> &statevector, int target_qubit);

#endif