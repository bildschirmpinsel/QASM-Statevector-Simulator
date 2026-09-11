#include "simulator.hpp"
#include <vector>

int main() {
  std::vector<Complex> statevector;
  // TODO add parameter for input file path
  const char *qasm_path = "./grover_n3_orig.qasm";
  auto gates = parseQASM(qasm_path, statevector);

  statevector[0] = ONE;
#ifndef NDEBUG
  std::cout << "Initialized statevector to |";
  for (int i = 0; i < number_qubits; i++) {
    std::cout << "0";
  }
  std::cout << ">" << std::endl;
#endif

  simulate(gates, statevector);

  // TODO add flag to toggle measurement
  // for (int qubit = 0; qubit < number_qubits; qubit++) {
  //   measure(statevector, qubit);
  // }

  for (int i = 0; i < statevector.size(); i++) {
    std::cout << statevector[i] << std::endl;
  }

  return 0;
}

std::vector<QASMGate> parseQASM(const char *file,
                                std::vector<Complex> &statevector) {
  std::vector<QASMGate> parsed_gates;
  std::unordered_map<std::string, int> registerToVector;

#ifndef NDEBUG
  std::cout << "Parsing file at path " << file << std::endl;
#endif

  std::ifstream qasm_file(file);
  std::string line;
  // command either register declaration or gate application
  std::string command;
  // captures whatever comes after command
  std::string definition;
  while (std::getline(qasm_file, line)) {
    std::istringstream stream(line);
    stream >> command;

    if (command == "qreg") {
      // new quantum register
      std::string register_name;
      int size = 0;
      stream >> definition;
      enum ParserState { NAME, INDEX };
      ParserState parser_state = ParserState::NAME;

      // parse register name and size from definition
      for (char c : definition) {
        if (parser_state == ParserState::NAME && c != '[') {
          // read until opening bracket
          register_name.push_back(c);
        } else if (parser_state == ParserState::INDEX && c != ']') {
          // parse integer until closing bracket
          size = size * 10 + (c - '0');
        } else {
          if (c == '[') {
            parser_state = ParserState::INDEX;
          } else {
            // break when all information from register declaration has been
            // read
            break;
          }
        }
      }

      // map name to vector so that references can be found for operands
      registerToVector.insert({register_name, number_qubits});
#ifndef NDEBUG
      std::cout << "\tAdded register " << register_name << " starting at qubit "
                << number_qubits << std::endl;
#endif
      number_qubits += size;

    } else {
      // get rid of parameters
      command = command.substr(0, command.find('('));
      std::vector<int> operand_qubits;
      auto gate = gateMap.find(command);
      if (gate != gateMap.end()) {
#ifndef NDEBUG
        std::cout << "\tParsing gate " << command << std::endl;
#endif
        enum ParserState { PARAMETER, OPERAND, INDEX, IDLE };
        ParserState parser_state = ParserState::IDLE;
        std::string operand_name;
        int operand_index;

        stream >> definition;
        for (char c : definition) {
          if (c == '(') {
#ifndef NDEBUG
            std::cout << "\t\tStart parsing parameter..." << std::endl;
#endif
            // parameters are always after gate name and in parenthesis
            parser_state = ParserState::PARAMETER;
          } else if (c == ')') {
#ifndef NDEBUG
            std::cout << "\t\tFinished parsing parameter with value: "
                      << std::endl;
#endif
            // terminating clause for parameter parsing
            // delimiter for parameters
            parser_state = ParserState::IDLE;
          } else if (parser_state == ParserState::PARAMETER) {
            // TODO parse parameters
          } else if (c == ' ' || c == ',' || c == ';') {
            parser_state = ParserState::OPERAND;
            if (c == ',' || c == ';') {
// parsed one of multiple operands, add global index in
// statevector
#ifndef NDEBUG
              std::cout << "\t\tFinished parsing operand " << operand_name
                        << "[" << operand_index << "]" << std::endl;
#endif
              operand_qubits.emplace_back(registerToVector.at(operand_name) +
                                          operand_index);
              // clear operand name
              operand_name.clear();
            }
#ifndef NDEBUG
            else {
              std::cout << "\t\tStart parsing operand..." << std::endl;
            }
#endif
          } else if (parser_state == ParserState::IDLE &&
                     (c >= 'a' && c <= 'z' || c >= 'A' && c <= 'Z')) {
            parser_state = ParserState::OPERAND;
            operand_name.push_back(c);
          } else if (c == '[') {
            // terminating clause for operand parsing
            // indices are always in brackets
            parser_state = ParserState::INDEX;
            operand_index = 0;
          } else if (parser_state == ParserState::OPERAND) {
            // build operand string
            operand_name.push_back(c);
          } else if (c == ']') {
            // terminating clause for index parsing
            // delimiter for indices
            parser_state = ParserState::IDLE;
          } else if (parser_state == ParserState::INDEX) {
            operand_index = operand_index * 10 + (c - '0');
          }
        }

        // fill in QASM gate
        QASMGate qasm_gate;
        qasm_gate.gate = gate->second;
        switch (operand_qubits.size()) {
        case 3:
          qasm_gate.qubit_3_global_index = operand_qubits[2];
          [[fallthrough]];
        case 2:
          qasm_gate.qubit_2_global_index = operand_qubits[1];
          [[fallthrough]];
        case 1:
          qasm_gate.qubit_1_global_index = operand_qubits[0];
          break;
        }
        parsed_gates.push_back(qasm_gate);
#ifndef NDEBUG
        std::cout << "\tFinished parsing gate" << std::endl;
#endif
      }
    }
  }

#ifndef NDEBUG
  std::cout << "\tNumber of qubits for simulator: " << number_qubits
            << std::endl;
#endif

  statevector.resize(1u << number_qubits);
  return parsed_gates;
}

void simulate(std::vector<QASMGate> gates, std::vector<Complex> &statevector) {
  for (auto gate : gates) {
    auto gate_id = gate.gate;
    if (gate_id == GateID::CZ || gate_id == GateID::CX ||
        gate_id == GateID::CCX) {
      if (gate_id == GateID::CX || gate_id == GateID::CZ) {
        std::vector<Complex> statevector_copy(statevector);
        applyGate(UNITARY_ACTIVE, gate.qubit_1_global_index, statevector);
        if (gate_id == GateID::CX) {
          // CX
          applyGate(getUnitary(GateID::X, nullptr), gate.qubit_2_global_index,
                    statevector);
        } else {
          // CZ
          applyGate(getUnitary(GateID::Z, nullptr), gate.qubit_2_global_index,
                    statevector);
        }
        applyGate(UNITARY_INACTIVE, gate.qubit_1_global_index,
                  statevector_copy);
        for (int i = 0; i < statevector.size(); i++) {
          statevector[i] += statevector_copy[i];
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
        simulate(ccx_decomposition, statevector);
      }
    } else {
      // single unitary
      auto unitary = getUnitary(gate_id, nullptr);
      applyGate(unitary, gate.qubit_1_global_index, statevector);
    }
  }
}

void measure(std::vector<Complex> &statevector, int target_qubit) {
  std::vector<Complex> statevector_copy(statevector);
  applyGate(UNITARY_INACTIVE, target_qubit, statevector_copy);

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
    applyGate(UNITARY_ACTIVE, target_qubit, statevector);
    const Complex scalar{1 / std::sqrt(p_1), 0};
    for (int i = 0; i < statevector.size(); i++)
      statevector[i] *= scalar;
  } else {
    // coin is 0
    applyGate(UNITARY_INACTIVE, target_qubit, statevector);
    const Complex scalar{1 / std::sqrt(p_0), 0};
    for (int i = 0; i < statevector.size(); i++)
      statevector[i] *= scalar;
  }
}

void applyGate(const Matrix unitary, int target_qubit,
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
