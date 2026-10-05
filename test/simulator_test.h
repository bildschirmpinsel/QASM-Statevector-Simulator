#include "../src/simulator.hpp"
#include <cmath>
#include <cxxtest/TestSuite.h>
#include <limits>
#include <vector>

class SimulatorTestSuite : public CxxTest::TestSuite {
private:
  std::vector<Complex> ket_zero;
  std::vector<Complex> ket_one;

  const double EPSILON = 1e-12;

  void testSingleQubitGate(const GateID gate,
                           const std::vector<Complex> expected_result_ket_zero,
                           const std::vector<Complex> expected_result_ket_one) {
    const Matrix unitary = getUnitary(gate, 0.0);
    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testSingelQubitRotationGate(
      const Matrix unitary, const std::vector<Complex> expected_result_ket_zero,
      const std::vector<Complex> expected_result_ket_one) {
    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++) {
      TS_ASSERT_DELTA(ket_zero[i].real(), expected_result_ket_zero[i].real(),
                      EPSILON);
      TS_ASSERT_DELTA(ket_zero[i].imag(), expected_result_ket_zero[i].imag(),
                      EPSILON);
    }
    for (int i = 0; i < 2; i++) {
      TS_ASSERT_DELTA(ket_one[i].real(), expected_result_ket_one[i].real(),
                      EPSILON);
      TS_ASSERT_DELTA(ket_one[i].imag(), expected_result_ket_one[i].imag(),
                      EPSILON);
    }
  }

  void testConditionalGate(const GateID gate,
                           const std::vector<Complex> expected_result_1,
                           const std::vector<Complex> expected_result_2,
                           std::vector<Complex> input_1,
                           std::vector<Complex> input_2) {
    const Matrix unitary = getUnitary(gate, 0.0);
    const int target_qubit = 0;
    const int control_qubit = 1;
    const int number_qubits = 2;

    applyControlledGate(unitary, control_qubit, target_qubit, number_qubits,
                        input_1);
    applyControlledGate(unitary, control_qubit, target_qubit, number_qubits,
                        input_2);

    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(input_1[i], expected_result_1[i]);
    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(input_2[i], expected_result_2[i]);
  }

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

    testSingleQubitGate(GateID::X, expected_result_ket_zero,
                        expected_result_ket_one);
  }

  void testGatePauliY() {
    const std::vector<Complex> expected_result_ket_zero{ZERO, I};
    const std::vector<Complex> expected_result_ket_one{-I, ZERO};

    testSingleQubitGate(GateID::Y, expected_result_ket_zero,
                        expected_result_ket_one);
  }

  void testGatePauliZ() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{ZERO, -ONE};

    testSingleQubitGate(GateID::Z, expected_result_ket_zero,
                        expected_result_ket_one);
  }

  void testGateHadamard() {
    const std::vector<Complex> expected_result_ket_zero{INV_SQRT_2, INV_SQRT_2};
    const std::vector<Complex> expected_result_ket_one{INV_SQRT_2,
                                                       MINUS_ONE * INV_SQRT_2};

    testSingleQubitGate(GateID::H, expected_result_ket_zero,
                        expected_result_ket_one);
  }

  void testGateS() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{ZERO, I};

    testSingleQubitGate(GateID::S, expected_result_ket_zero,
                        expected_result_ket_one);
  }

  void testGateT() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{
        ZERO, std::polar(1.0, M_PI / 4.0)};

    testSingleQubitGate(GateID::T, expected_result_ket_zero,
                        expected_result_ket_one);
  }

  void testGateAdjointT() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{
        ZERO, std::polar(1.0, -M_PI / 4.0)};

    testSingleQubitGate(GateID::AT, expected_result_ket_zero,
                        expected_result_ket_one);
  }

  void testGateRotationXZeroAngle() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{ZERO, ONE};

    testSingleQubitGate(GateID::RX, expected_result_ket_zero,
                        expected_result_ket_one);
  }

  void testGateRotationYZeroAngle() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{ZERO, ONE};

    testSingleQubitGate(GateID::RY, expected_result_ket_zero,
                        expected_result_ket_one);
  }

  void testGateRotationZZeroAngle() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{ZERO, ONE};

    testSingleQubitGate(GateID::RZ, expected_result_ket_zero,
                        expected_result_ket_one);
  }

  void testGateRotationXFullRotation() {
    const std::vector<Complex> expected_result_ket_zero{Complex{-1.0, -0.0},
                                                        Complex{0.0, -0.0}};
    const std::vector<Complex> expected_result_ket_one{Complex{0.0, -0.0},
                                                       Complex{-1.0, -0.0}};

    const Matrix unitary = getUnitary(GateID::RX, 2 * M_PI);

    testSingelQubitRotationGate(unitary, expected_result_ket_zero,
                                expected_result_ket_one);
  }

  void testGateRotationYFullRotation() {
    const std::vector<Complex> expected_result_ket_zero{Complex{-1.0, -0.0},
                                                        Complex{0.0, -0.0}};
    const std::vector<Complex> expected_result_ket_one{Complex{0.0, -0.0},
                                                       Complex{-1.0, -0.0}};

    const Matrix unitary = getUnitary(GateID::RY, 2 * M_PI);

    testSingelQubitRotationGate(unitary, expected_result_ket_zero,
                                expected_result_ket_one);
  }

  void testGateRotationZFullRotation() {
    const std::vector<Complex> expected_result_ket_zero{Complex{-1.0, -0.0},
                                                        Complex{0.0, -0.0}};
    const std::vector<Complex> expected_result_ket_one{Complex{0.0, -0.0},
                                                       Complex{-1.0, -0.0}};

    const Matrix unitary = getUnitary(GateID::RZ, 2 * M_PI);

    testSingelQubitRotationGate(unitary, expected_result_ket_zero,
                                expected_result_ket_one);
  }

  void testGateRotationXQuarterRotation() {
    const std::vector<Complex> expected_result_ket_zero{
        INV_SQRT_2, Complex{0.0, -1.0} * INV_SQRT_2};
    const std::vector<Complex> expected_result_ket_one{
        Complex{0.0, -1.0} * INV_SQRT_2, INV_SQRT_2};

    const Matrix unitary = getUnitary(GateID::RX, 0.5 * M_PI);

    testSingelQubitRotationGate(unitary, expected_result_ket_zero,
                                expected_result_ket_one);
  }

  void testGateRotationYQuarterRotation() {
    const std::vector<Complex> expected_result_ket_zero{INV_SQRT_2, INV_SQRT_2};
    const std::vector<Complex> expected_result_ket_one{-INV_SQRT_2, INV_SQRT_2};

    const Matrix unitary = getUnitary(GateID::RY, 0.5 * M_PI);

    testSingelQubitRotationGate(unitary, expected_result_ket_zero,
                                expected_result_ket_one);
  }

  void testGateRotationZQuarterRotation() {
    const std::vector<Complex> expected_result_ket_zero{
        Complex{cos(0.25 * M_PI), -sin(0.25 * M_PI)}, ZERO};
    const std::vector<Complex> expected_result_ket_one{
        ZERO, Complex{cos(0.25 * M_PI), sin(0.25 * M_PI)}};

    const Matrix unitary = getUnitary(GateID::RZ, 0.5 * M_PI);

    testSingelQubitRotationGate(unitary, expected_result_ket_zero,
                                expected_result_ket_one);
  }

  void testGateConditionalXInactive() {
    const std::vector<Complex> expected_result_ket_00{ONE, ZERO, ZERO, ZERO};
    const std::vector<Complex> expected_result_ket_01{ZERO, ZERO, ONE, ZERO};

    std::vector<Complex> inactive_ket_00{ONE, ZERO, ZERO, ZERO};
    std::vector<Complex> inactive_ket_01{ZERO, ZERO, ONE, ZERO};

    testConditionalGate(GateID::X, expected_result_ket_00,
                        expected_result_ket_01, inactive_ket_00,
                        inactive_ket_01);
  }

  void testGateConditionalXActive() {
    const std::vector<Complex> expected_result_ket_10{ZERO, ZERO, ZERO, ONE};
    const std::vector<Complex> expected_result_ket_11{ZERO, ONE, ZERO, ZERO};

    std::vector<Complex> active_ket_10{ZERO, ONE, ZERO, ZERO};
    std::vector<Complex> active_ket_11{ZERO, ZERO, ZERO, ONE};

    testConditionalGate(GateID::X, expected_result_ket_10,
                        expected_result_ket_11, active_ket_10, active_ket_11);
  }

  void testGateConditionalZInactive() {
    const std::vector<Complex> expected_result_ket_00{ONE, ZERO, ZERO, ZERO};
    const std::vector<Complex> expected_result_ket_01{ZERO, ZERO, ONE, ZERO};

    std::vector<Complex> inactive_ket_00{ONE, ZERO, ZERO, ZERO};
    std::vector<Complex> inactive_ket_01{ZERO, ZERO, ONE, ZERO};

    testConditionalGate(GateID::Z, expected_result_ket_00,
                        expected_result_ket_01, inactive_ket_00,
                        inactive_ket_01);
  }

  void testGateConditionalZActive() {
    const std::vector<Complex> expected_result_ket_10{ZERO, ONE, ZERO, ZERO};
    const std::vector<Complex> expected_result_ket_11{ZERO, ZERO, ZERO, -ONE};

    std::vector<Complex> active_ket_10{ZERO, ONE, ZERO, ZERO};
    std::vector<Complex> active_ket_11{ZERO, ZERO, ZERO, ONE};

    testConditionalGate(GateID::Z, expected_result_ket_10,
                        expected_result_ket_11, active_ket_10, active_ket_11);
  }

  void testGateDoubleConditionalX() {}

  void testMeasurementKetZero() {}

  void testMeasurementKetOne() {}
};