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

  void testParseParameterNaturalNumbers() {
    const double natural_number_expected = 3;
    const double negative_number_expected = -3;

    const double addition_expected = 15;
    const double subtraction_expected = 8;
    const double multiplication_expected = 21;
    const double division_expected = 3;

    const std::string natural_number = "(3)";
    const std::string negative_number = "(-3)";

    const std::string addition = "(8+7)";
    const std::string subtraction = "(15-7)";
    const std::string multiplication = "(3*7)";
    const std::string division = "(21/7)";

    TS_ASSERT_EQUALS(parseParameter(natural_number), natural_number_expected);
    TS_ASSERT_EQUALS(parseParameter(negative_number), negative_number_expected);

    TS_ASSERT_EQUALS(parseParameter(addition), addition_expected);
    TS_ASSERT_EQUALS(parseParameter(subtraction), subtraction_expected);
    TS_ASSERT_EQUALS(parseParameter(multiplication), multiplication_expected);
    TS_ASSERT_EQUALS(parseParameter(division), division_expected);
  }

  void testParseParameterRationalNumbers() {
    const double rational_number_expected = 2.0 / 5.0;
    const double negative_rational_number_expected = -2.0 / 5.0;

    const std::string rational_number = "(2/5)";
    const std::string negative_rational_number = "(-2/5)";

    TS_ASSERT_EQUALS(parseParameter(rational_number), rational_number_expected);
    TS_ASSERT_EQUALS(parseParameter(negative_rational_number),
                     negative_rational_number_expected);
  }

  void testParseParameterDecimalNumbers() {
    const double decimal_number_expected = 3.12;
    const double negative_decimal_number_expected = -3.12;

    const double addition_expected = 5.75;
    const double subtraction_expected = 3.25;
    const double multiplication_expected = 6.25;
    const double division_expected = 2.5;

    const std::string decimal_number = "(3.12)";
    const std::string negative_decimal_number = "(-3.12)";

    const std::string addition = "(3.25+2.5)";
    const std::string subtraction = "(5.75-2.5)";
    const std::string multiplication = "(2.5*2.5)";
    const std::string division = "(6.25/2.5)";

    TS_ASSERT_EQUALS(parseParameter(decimal_number), decimal_number_expected);
    TS_ASSERT_EQUALS(parseParameter(negative_decimal_number),
                     negative_decimal_number_expected);

    TS_ASSERT_EQUALS(parseParameter(addition), addition_expected);
    TS_ASSERT_EQUALS(parseParameter(subtraction), subtraction_expected);
    TS_ASSERT_EQUALS(parseParameter(multiplication), multiplication_expected);
    TS_ASSERT_EQUALS(parseParameter(division), division_expected);
  }

  void testParseParameterIrrationalNumbers() {
    const double irrational_number_expected = M_PI;

    const std::string irrational_number = "(pi)";

    TS_ASSERT_EQUALS(parseParameter(irrational_number),
                     irrational_number_expected);
  }
};