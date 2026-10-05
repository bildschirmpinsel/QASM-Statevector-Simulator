#include "../src/qasm_parser.hpp"
#include <cxxtest/TestSuite.h>

class SimulatorTestSuite : public CxxTest::TestSuite {
public:
  void testParseRegisterDefinition() {
    const std::string expected_register_name{"up"};
    const unsigned int expected_register_size = 3;

    const std::string register_definition{"up[3];"};

    unsigned int size = 0;
    const auto result = parseRegisterDefinition(register_definition, size);

    TS_ASSERT_EQUALS(size, expected_register_size);
    TS_ASSERT_EQUALS(result, expected_register_name);
  }
};