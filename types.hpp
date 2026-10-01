#ifndef TYPES
#define TYPES

#include <complex>

using Complex = std::complex<double>;

enum GateID { X, Y, Z, H, S, T, AT, RX, RY, RZ, CX, CZ, CCX };

typedef struct {
  double rotation_degree;
  int qubit_1_global_index;
  int qubit_2_global_index;
  int qubit_3_global_index;
  GateID gate;
} QASMGate;

#endif