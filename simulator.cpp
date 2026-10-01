#include "simulator.hpp"
#include "types.hpp"

int main() {
  std::vector<Complex> statevector;
  // TODO add parameter for input file path
  const char *qasm_path = "grover_n5_orig.qasm";
  std::vector<QASMGate> gates;
  const int number_qubits = parseQASM(qasm_path, statevector, gates);

  statevector[0] = ONE;
#ifndef NDEBUG
  std::cout << "Initialized statevector to |";
  for (int i = 0; i < number_qubits; i++) {
    std::cout << "0";
  }
  std::cout << ">" << std::endl;
#endif

  simulate(gates, statevector, number_qubits);

  std::cout << "Final statevector:" << std::endl;
  for (auto x : statevector) {
    std::cout << "\t" << x << std::endl;
  }

  for (int qubit = 0; qubit < number_qubits; qubit++) {
    measure(statevector, qubit, number_qubits);
  }
  auto measured_qubit_index_iterator =
      std::find_if(statevector.begin(), statevector.end(),
                   [](const Complex x) { return x.real() > 0.99; });
  int measured_qubit_index =
      std::distance(statevector.begin(), measured_qubit_index_iterator);
  std::string measured_qubit_string;
  for (int i = number_qubits - 1; i >= 0; --i) {
    measured_qubit_string += ((measured_qubit_index >> i) & 1) ? '1' : '0';
  }
  std::cout << "Measured statevector: |" << measured_qubit_string << ">"
            << std::endl;

  return 0;
}

void simulate(std::vector<QASMGate> gates, std::vector<Complex> &statevector,
              const int number_qubits) {
  unsigned int processed_gates = 0;
  for (auto gate : gates) {
    auto gate_id = gate.gate;
    if (gate_id == GateID::CZ || gate_id == GateID::CX ||
        gate_id == GateID::CCX) {
      if (gate_id == GateID::CX || gate_id == GateID::CZ) {
        if (gate_id == GateID::CX) {
          // CX
          applyControlledGate(
              getUnitary(GateID::X, 0.0), gate.qubit_1_global_index,
              gate.qubit_2_global_index, number_qubits, statevector);
        } else {
          // CZ
          applyControlledGate(
              getUnitary(GateID::Z, 0.0), gate.qubit_1_global_index,
              gate.qubit_2_global_index, number_qubits, statevector);
        }
      } else if (gate_id == GateID::CCX) {
        std::vector<QASMGate> ccx_decomposition;
        // H(3)
        ccx_decomposition.emplace_back(
            QASMGate{0.0, gate.qubit_3_global_index, 0, 0, GateID::H});
        // CNOT(2,3)
        ccx_decomposition.emplace_back(QASMGate{0.0, gate.qubit_2_global_index,
                                                gate.qubit_3_global_index, 0,
                                                GateID::CX});
        // adjoint T(3)
        ccx_decomposition.emplace_back(
            QASMGate{0.0, gate.qubit_3_global_index, 0, 0, GateID::AT});
        // CNOT(1,3)
        ccx_decomposition.emplace_back(QASMGate{0.0, gate.qubit_1_global_index,
                                                gate.qubit_3_global_index, 0,
                                                GateID::CX});
        // T(3)
        ccx_decomposition.emplace_back(
            QASMGate{0.0, gate.qubit_3_global_index, 0, 0, GateID::T});
        // CNOT(2,3)
        ccx_decomposition.emplace_back(QASMGate{0.0, gate.qubit_2_global_index,
                                                gate.qubit_3_global_index, 0,
                                                GateID::CX});
        // adjoint T(3)
        ccx_decomposition.emplace_back(
            QASMGate{0.0, gate.qubit_3_global_index, 0, 0, GateID::AT});
        // CNOT(1,3)
        ccx_decomposition.emplace_back(QASMGate{0.0, gate.qubit_1_global_index,
                                                gate.qubit_3_global_index, 0,
                                                GateID::CX});
        // T(2)
        ccx_decomposition.emplace_back(
            QASMGate{0.0, gate.qubit_2_global_index, 0, 0, GateID::T});
        // T(3)
        ccx_decomposition.emplace_back(
            QASMGate{0.0, gate.qubit_3_global_index, 0, 0, GateID::T});
        // H(3)
        ccx_decomposition.emplace_back(
            QASMGate{0.0, gate.qubit_3_global_index, 0, 0, GateID::H});
        // CNOT(1,2)
        ccx_decomposition.emplace_back(QASMGate{0.0, gate.qubit_1_global_index,
                                                gate.qubit_2_global_index, 0,
                                                GateID::CX});
        // adjoint T(2)
        ccx_decomposition.emplace_back(
            QASMGate{0.0, gate.qubit_2_global_index, 0, 0, GateID::AT});
        // T(1)
        ccx_decomposition.emplace_back(
            QASMGate{0.0, gate.qubit_1_global_index, 0, 0, GateID::T});
        // CNOT(1,2)
        ccx_decomposition.emplace_back(QASMGate{0.0, gate.qubit_1_global_index,
                                                gate.qubit_2_global_index, 0,
                                                GateID::CX});
        simulate(ccx_decomposition, statevector, number_qubits);
      }
    } else {
      // single unitary
      auto unitary = getUnitary(gate_id, gate.rotation_degree);
      applyGate(unitary, gate.qubit_1_global_index, number_qubits, statevector);
    }
    processed_gates++;
  }
}

void measure(std::vector<Complex> &statevector, int target_qubit,
             const int number_qubits) {
  std::vector<Complex> statevector_copy(statevector);
  applyGate(UNITARY_INACTIVE, target_qubit, number_qubits, statevector_copy);

  double p_0 = 0.0;
  for (const auto &z : statevector_copy) {
    p_0 += std::norm(z);
  }
  const double p_1 = 1.0 - p_0;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::bernoulli_distribution binomial(p_1);
  bool coin = binomial(gen);

  if (coin) {
    // coin is 1
    applyGate(UNITARY_ACTIVE, target_qubit, number_qubits, statevector);
    const Complex scalar{1 / std::sqrt(p_1), 0};
    for (int i = 0; i < statevector.size(); i++)
      statevector[i] *= scalar;
  } else {
    // coin is 0
    applyGate(UNITARY_INACTIVE, target_qubit, number_qubits, statevector);
    const Complex scalar{1 / std::sqrt(p_0), 0};
    for (int i = 0; i < statevector.size(); i++)
      statevector[i] *= scalar;
  }
}

void applyGate(const Matrix unitary, int target_qubit, const int number_qubits,
               std::vector<Complex> &statevector) {

  const unsigned int target_index = 1u << (target_qubit);
  const unsigned int stride = 1u << (number_qubits - target_qubit - 1);

  for (unsigned int i = 0; i < target_index; i++) {
    const unsigned int subsection_start_index = i * stride * 2;
    const unsigned int dividing_index = subsection_start_index + stride;
    const unsigned int subsection_end_index = dividing_index + i * stride;

#pragma OMP PARALLEL FOR
    for (unsigned int j = 0; j < stride; j++) {
      const auto upper_value = statevector[subsection_start_index + j];
      const auto lower_value = statevector[dividing_index + j];
      statevector[subsection_start_index + j] =
          unitary[0] * upper_value + unitary[1] * lower_value;
      statevector[dividing_index + j] =
          unitary[2] * upper_value + unitary[3] * lower_value;
    }
  }
}

void applyControlledGate(const Matrix unitary, int control_qubit,
                         int target_qubit, const int number_qubits,
                         std::vector<Complex> &statevector) {
  const unsigned int control_mask = 1u << (number_qubits - 1 - control_qubit);
  const unsigned int target_mask = 1u << (number_qubits - 1 - target_qubit);

#pragma OMP PARALLEL FOR
  for (unsigned int i = 0; i < statevector.size(); i++) {
    if (i & target_mask) {
      continue;
    }

    if ((i & control_mask) == 0) {
      continue;
    }

    const std::size_t target_zero_index = i;
    const std::size_t target_one_index = i | target_mask;

    const Complex upper = statevector[target_zero_index];
    const Complex lower = statevector[target_one_index];

    statevector[target_zero_index] = unitary[0] * upper + unitary[1] * lower;

    statevector[target_one_index] = unitary[2] * upper + unitary[3] * lower;
  }
}