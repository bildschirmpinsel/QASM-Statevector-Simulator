#include "qasm_parser.hpp"
#include "simulator.hpp"

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <qasm_file_path> <output_file_path>"
              << std::endl;
    return 1;
  }

  std::vector<Complex> statevector;
  std::vector<QASMGate> gates;
  const int number_qubits = parseQASM(argv[1], statevector, gates);

  statevector[0] = ONE;
#ifndef NDEBUG
  std::cout << "Initialized statevector to |";
  for (int i = 0; i < number_qubits; i++) {
    std::cout << "0";
  }
  std::cout << ">" << std::endl;
#endif

  unsigned int number_number_processed_gates = 0;
  unsigned int number_total_gates = 0;

  simulate(gates, statevector, number_qubits, number_number_processed_gates,
           number_total_gates);

  const char *output_file_path = argv[2];
#ifndef NDEBUG
  std::cout << std::endl
            << "Putting result in file at path " << output_file_path
            << std::endl;
#endif

  std::ofstream output_file(output_file_path, std::ios::app);
  output_file << "Final statevector:" << std::endl;
  for (auto x : statevector) {
    output_file << "\t" << x << std::endl;
  }

  for (unsigned int qubit = 0; qubit < number_qubits; qubit++) {
    measure(statevector, qubit, number_qubits);
  }

  // print collapsed statevector
  auto measured_qubit_index_iterator =
      std::max_element(statevector.begin(), statevector.end(),
                       [](const Complex &a, const Complex &b) {
                         return std::norm(a) < std::norm(b);
                       });
  int measured_qubit_index =
      std::distance(statevector.begin(), measured_qubit_index_iterator);
  std::string measured_qubit_string;
  for (int i = 0; i < number_qubits; i++) {
    measured_qubit_string += ((measured_qubit_index >> i) & 1) ? '1' : '0';
  }
  std::cout << std::endl
            << "Measured statevector (LSB first): |" << measured_qubit_string
            << ">" << std::endl;

  return 0;
}