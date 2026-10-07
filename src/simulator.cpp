#include "simulator.hpp"

void simulate(std::vector<QASMGate> gates, std::vector<Complex> &statevector,
              const unsigned int number_qubits,
              unsigned int &number_processed_gates,
              unsigned int &number_total_gates, const bool print_progress) {
  number_total_gates += gates.size();
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
        simulate(ccx_decomposition, statevector, number_qubits,
                 number_processed_gates, number_total_gates, print_progress);
      }
    } else {
      // single unitary
      auto unitary = getUnitary(gate_id, gate.rotation_degree);
      applyGate(unitary, gate.qubit_1_global_index, number_qubits, statevector);
    }

    number_processed_gates++;

    double norm = 0.0;
    #pragma OMP PARALLEL FOR REDUCE(+:norm) 
    for (auto z : statevector)
      norm += std::norm(z);

    assert(std::abs(norm - 1.0) < 0.01 && "Norm of statevector is not close to 1!");

    if (print_progress) {
      // print progress bar
      constexpr unsigned int width = 50;

      const auto current =
          std::clamp(number_processed_gates, 0u, number_total_gates);
      const int completed =
          std::clamp((current * width) / number_total_gates, 0u, width);

      const int remaining = width - completed;

      std::cout << '\r' << '['
                << std::string(static_cast<std::size_t>(completed), '=')
                << std::string(static_cast<std::size_t>(remaining), ' ') << "] "
                << (current * 100 / number_total_gates) << "% (" << current
                << '/' << number_total_gates << ')' << std::flush;
    }
  }
}

void measure(std::vector<Complex> &statevector, unsigned int target_qubit,
             const unsigned int number_qubits) {
  double p_0 = 0.0;

  // apply |0><0| to all values in state vector that have
  // the target bit inactive
  const unsigned int target_mask = 1u << (number_qubits - 1 - target_qubit);
  for (unsigned int i = 0; i < statevector.size(); i++) {
    if ((i & target_mask) == 0) {
      p_0 += std::norm(statevector[i]);
    }
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
#pragma OMP PARALLEL FOR
    for (int i = 0; i < statevector.size(); i++)
      statevector[i] *= scalar;
  } else {
    // coin is 0
    applyGate(UNITARY_INACTIVE, target_qubit, number_qubits, statevector);
    const Complex scalar{1 / std::sqrt(p_0), 0};
#pragma OMP PARALLEL FOR
    for (int i = 0; i < statevector.size(); i++)
      statevector[i] *= scalar;
  }
}

void applyGate(const Matrix unitary, const unsigned int target_qubit,
               const unsigned int number_qubits,
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

void applyControlledGate(const Matrix unitary, const unsigned int control_qubit,
                         const unsigned int target_qubit,
                         const unsigned int number_qubits,
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

    const unsigned int target_inactive_index = i;
    const unsigned int target_active_index = i | target_mask;

    const Complex inactive = statevector[target_inactive_index];
    const Complex active = statevector[target_active_index];

    statevector[target_inactive_index] =
        unitary[0] * inactive + unitary[1] * active;

    statevector[target_active_index] =
        unitary[2] * inactive + unitary[3] * active;
  }
}