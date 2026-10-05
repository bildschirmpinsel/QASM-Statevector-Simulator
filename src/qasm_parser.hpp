#ifndef QASM_PARSER
#define QASM_PARSER

#include "types.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Map gate symbol to internal gate id.
static const std::unordered_map<std::string, GateID> gateMap{
    {"x", GateID::X},   {"y", GateID::Y},   {"z", GateID::Z},
    {"h", GateID::H},   {"s", GateID::S},   {"t", GateID::T},
    {"rx", GateID::RX}, {"ry", GateID::RY}, {"rz", GateID::RZ},
    {"cx", GateID::CX}, {"cz", GateID::CZ}, {"ccx", GateID::CCX}};

// Internally represent the parameter operations.
enum ParameterOperation { ADD, SUB, MUL, DIV, NONE };
// Map the corresponding operation symbol to the interal representation.
static const std::unordered_map<char, ParameterOperation> parameterOperationMap{
    {'+', ParameterOperation::ADD},
    {'-', ParameterOperation::SUB},
    {'*', ParameterOperation::MUL},
    {'/', ParameterOperation::DIV}};

/**
 Parse an openQASM 3 circuit.

 @param file The file path as a char array.
 @param statevector The statevector to be resized based on the given circuit's
 qubits.
 @param gates A list of QASMGate structs to be populated by the gates of this
 circuit.

 @return The number of qubits.
*/
unsigned int parseQASM(const char *file, std::vector<Complex> &statevector,
                       std::vector<QASMGate> &gates);

/**
 Parses the register definition in openQASM 3 format.

 @param definition The string that is the register definition.
 @param size Parameter that will contain the size of the defined register.

 @return A string containing the register name.
*/
std::string parseRegisterDefinition(const std::string definition, unsigned int &size);

/**
 Parses a given gate string into a QASMGate struct, extracting the type of gate,
 parameters, and operands in the process.

 @param gate_string The string that represents the gate.
 @param operand_string The string that represents the operands of the gate.
 @param registerToVector A map mapping register names to statevector offsets.

 @return The parsed gate as a QASMGate struct.
*/
QASMGate
parseGate(const std::string gate_string, const std::string operand_string,
          std::unordered_map<std::string, unsigned int> &registerToVector);

/**
 Parse a parameter out of a string that starts with '(' and ends with ')'. Also
 must at most contain two operands. Can process 'pi' for the value of pi. Can
 also process a leading minus to negate the first parameter.

 @param parameter_string The string containing the parameter definition.

 @return The parsed parameter.
 */
double parseParameter(const std::string parameter_string);

#endif