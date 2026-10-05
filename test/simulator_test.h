#include "../src/simulator.hpp"
#include <cxxtest/TestSuite.h>

class SimulatorTestSuite : public CxxTest::TestSuite {
private:
  std::vector<Complex> ket_zero;
  std::vector<Complex> ket_one;

public:
  void setUp() {
    ket_zero.push_back(ONE);
    ket_zero.push_back(ZERO);

    ket_one.push_back(ZERO);
    ket_one.push_back(ONE);
  }

  void tearDown() {
    ket_zero.clear();
    ket_one.clear();
  }

  void testGatePauliX() {
    const std::vector<Complex> expected_result_ket_zero{ZERO, ONE};
    const std::vector<Complex> expected_result_ket_one{ONE, ZERO};

    const Matrix unitary = getUnitary(GateID::X, 0.0);

    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testGatePauliY() {
    const std::vector<Complex> expected_result_ket_zero{ZERO, I};
    const std::vector<Complex> expected_result_ket_one{-I, ZERO};

    const Matrix unitary = getUnitary(GateID::Y, 0.0);

    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testGatePauliZ() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{ZERO, -ONE};

    const Matrix unitary = getUnitary(GateID::Z, 0.0);

    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testGateHadamard() {
    const std::vector<Complex> expected_result_ket_zero{INV_SQRT_2, INV_SQRT_2};
    const std::vector<Complex> expected_result_ket_one{INV_SQRT_2,
                                                       MINUS_ONE * INV_SQRT_2};

    const Matrix unitary = getUnitary(GateID::H, 0.0);

    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testGateS() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{ZERO, I};

    const Matrix unitary = getUnitary(GateID::S, 0.0);

    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testGateT() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{
        ZERO, std::polar(1.0, M_PI / 4.0)};

    const Matrix unitary = getUnitary(GateID::T, 0.0);

    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testGateAdjointT() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{
        ZERO, std::polar(1.0, -M_PI / 4.0)};

    const Matrix unitary = getUnitary(GateID::AT, 0.0);

    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testGateRotationX() {}

  void testGateRotationY() {}

  void testGateRotationZ() {}

  void testGateConditionalX() {}

  void testGateConditionalZ() {}

  void testGateDoubleConditionalX() {}

  void testMeasurementKetZero() {}

  void testMeasurementKetOne() {}
};