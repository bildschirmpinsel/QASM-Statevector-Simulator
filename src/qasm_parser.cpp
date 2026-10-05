#include "qasm_parser.hpp"
#include "types.hpp"
#include <stdexcept>

unsigned int parseQASM(const char *file, std::vector<Complex> &statevector,
                       std::vector<QASMGate> &gates) {
  std::unordered_map<std::string, unsigned int> registerToVector;
  unsigned int number_qubits = 0;

#ifndef NDEBUG
  std::cout << "Parsing file at path " << file << std::endl;
#endif

  std::ifstream qasm_file(file);
  std::string line;
  // command either register declaration or gate application
  std::string command;
  while (std::getline(qasm_file, line)) {
    std::istringstream stream(line);
    stream >> command;

    if (command == "OPENQASM") {
      continue;
    } else if (command == "qreg") {
      // new quantum register
      unsigned int size = 0;
      std::string definition;
      stream >> definition;
      auto register_name = parseRegisterDefinition(definition, size);

      // map name to vector so that references can be found for operands
      registerToVector.insert({register_name, number_qubits});
#ifndef NDEBUG
      std::cout << "\tAdded register " << register_name << " starting at qubit "
                << number_qubits << std::endl;
#endif
      number_qubits += size;

    } else {
      std::string operand_string;
      stream >> operand_string;

      gates.push_back(parseGate(command, operand_string, registerToVector));
    }
  }

#ifndef NDEBUG
  std::cout << "\tNumber of qubits for simulator: " << number_qubits
            << std::endl;
#endif

  statevector.resize(1u << number_qubits);
  return number_qubits;
}

std::string parseRegisterDefinition(const std::string definition,
                                    unsigned int &size) {
  std::string register_name;
  size = 0;
  enum ParserState { NAME, SIZE };
  ParserState parser_state = ParserState::NAME;

  // parse register name and size from definition
  for (char c : definition) {
    switch (parser_state) {
    case NAME:
      if (c == '[') {
        parser_state = ParserState::SIZE;
      } else {
        register_name.push_back(c);
      }
      break;
    case SIZE:
      if (c == ']') {
        return register_name;
      } else {
        size = size * 10 + (c - '0');
      }
      break;
    }
  }

  throw std::invalid_argument(
      "Declaration of register does not terminate with \']\'!");
}

QASMGate
parseGate(const std::string gate_string, const std::string operand_string,
          std::unordered_map<std::string, unsigned int> &registerToVector) {
  auto gate_name = gate_string.substr(0, gate_string.find('('));
#ifndef NDEBUG
  std::cout << "\tParsing gate " << gate_name << std::endl;
#endif
  double parameter = 0.0;
  if (gate_name.compare(gate_string)) {
    // parse parameter iff there was a parameter
    // to be trimmed away to extract the gate name
    parameter = parseParameter(gate_string);
  }
  std::vector<int> operand_qubits;
  auto gate = gateMap.find(gate_name);
  if (gate == gateMap.end()) {
    throw std::invalid_argument("Gate name not supported!");
  }
  enum ParserState { OPERAND, INDEX, IDLE };
  ParserState parser_state = ParserState::IDLE;
  std::string operand_name;
  int operand_index;

  for (char c : operand_string) {
    switch (parser_state) {
    case IDLE:
      switch (c) {
      default:
        operand_name.push_back(c);
        parser_state = ParserState::OPERAND;
      case ' ':
      case ',':
#ifndef NDEBUG
        std::cout << "\t\tStart parsing operand..." << std::endl;
#endif
        parser_state = ParserState::OPERAND;
        break;
      case ';':
        break;
      }
      break;
    case OPERAND:
      switch (c) {
      case '[':
        operand_index = 0;
        parser_state = ParserState::INDEX;
        break;
      default:
        operand_name.push_back(c);
        break;
      }
      break;
    case INDEX:
      switch (c) {
      case ']':
#ifndef NDEBUG
        std::cout << "\t\tFinished parsing operand " << operand_name << "["
                  << operand_index << "]" << std::endl;
#endif
        operand_qubits.emplace_back(registerToVector.at(operand_name) +
                                    operand_index);
        // clear operand name
        operand_name.clear();
        parser_state = ParserState::IDLE;
        break;
      default:
        operand_index = operand_index * 10 + (c - '0');
        break;
      }
    }
  }
  // fill in QASM gate
  QASMGate qasm_gate;
  qasm_gate.gate = gate->second;
  qasm_gate.rotation_degree = parameter;
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
#ifndef NDEBUG
  std::cout << "\tFinished parsing gate" << std::endl;
#endif

  return qasm_gate;
}

/**
 This function can read parameters like (-pi/2.5) but nothing else
 that is not wrapped in braces and contains only a single operator
 (except a leading minus).

 Gate names are ignored in the string as long as they do not contain a p.
*/
double parseParameter(const std::string parameter_string) {
#ifndef NDEBUG
  std::cout << "\t\tStart parsing parameter..." << std::endl;
#endif

  ParameterOperation parameter_operation = ParameterOperation::NONE;
  double parameter = 0.0;
  int decimal_place = 0;
  double parameter_buffer = 0.0;

  auto execute_operation = [&]() {
    switch (parameter_operation) {
    case ADD:
      parameter += parameter_buffer;
      break;
    case SUB:
      parameter -= parameter_buffer;
      break;
    case MUL:
      parameter *= parameter_buffer;
      break;
    case DIV:
      parameter /= parameter_buffer;
      break;
    case NONE:
      break;
    }
    parameter_buffer = 0.0;
    decimal_place = 0;
    parameter_operation = ParameterOperation::NONE;
  };

  for (char c : parameter_string) {
    switch (c) {
    case ')':
      // consume last operation
      execute_operation();
      break;
    case '+':
      [[fallthrough]];
    case '-':
      [[fallthrough]];
    case '*':
      [[fallthrough]];
    case '/':
      if (parameter_operation == ParameterOperation::SUB) {
        parameter_buffer *= -1.0;
      }
      // parsing of first operand done, clear buffer
      parameter = parameter_buffer;
      parameter_buffer = 0.0;
      parameter_operation = parameterOperationMap.at(c);
      break;
    case 'p':
      parameter_buffer = M_PI;
    case 'i':
      // read pi on p, skip i
      // no imaginary unit should appear as a parameter
      break;
    case '.':
      decimal_place = 1;
      break;
    default:
      if (c >= '0' && c <= '9') {
        // convert char to number either before
        // decimal place or after
        if (decimal_place == 0) {
          parameter_buffer = parameter_buffer * 10 + (c - '0');
        } else {
          parameter_buffer += (c - '0') / (pow(10.0, decimal_place));
          decimal_place++;
        }
      }
      break;
    }
  }

#ifndef NDEBUG
  std::cout << "\t\tFinished parsing parameter with value: " << parameter
            << std::endl;
#endif
  return parameter;
}
