#ifndef QASM_PARSER
#define QASM_PARSER

#include "types.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

static const std::unordered_map<std::string, GateID> gateMap{
    {"x", GateID::X},   {"y", GateID::Y},   {"z", GateID::Z},
    {"h", GateID::H},   {"s", GateID::S},   {"t", GateID::T},
    {"rx", GateID::RX}, {"ry", GateID::RY}, {"rz", GateID::RZ},
    {"cx", GateID::CX}, {"cz", GateID::CZ}, {"ccx", GateID::CCX}};


int parseQASM(const char *file, std::vector<Complex> &statevector, std::vector<QASMGate> &gates);

std::string parseRegisterDefinition(std::string definition, int &size);

QASMGate parseGate(std::string gate_string, std::string operand_string, std::unordered_map<std::string, int> &registerToVector);

double parseParameter(std::string parameter_string);

#endif