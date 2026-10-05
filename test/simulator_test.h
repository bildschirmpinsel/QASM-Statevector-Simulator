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

  void testGateRotationXZeroAngle() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{ZERO, ONE};

    const Matrix unitary = getUnitary(GateID::RX, 0.0);

    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testGateRotationYZeroAngle() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{ZERO, ONE};

    const Matrix unitary = getUnitary(GateID::RY, 0.0);

    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testGateRotationZZeroAngle() {
    const std::vector<Complex> expected_result_ket_zero{ONE, ZERO};
    const std::vector<Complex> expected_result_ket_one{ZERO, ONE};

    const Matrix unitary = getUnitary(GateID::RZ, 0.0);

    applyGate(unitary, 0, 1, ket_zero);
    applyGate(unitary, 0, 1, ket_one);

    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_zero[i], expected_result_ket_zero[i]);
    for (int i = 0; i < 2; i++)
      TS_ASSERT_EQUALS(ket_one[i], expected_result_ket_one[i]);
  }

  void testGateRotationXFullRotation() {
    const std::vector<Complex> expected_result_ket_zero{Complex{-1.0, -0.0},
                                                        Complex{0.0, -0.0}};
    const std::vector<Complex> expected_result_ket_one{Complex{0.0, -0.0},
                                                       Complex{-1.0, -0.0}};

    const Matrix unitary = getUnitary(GateID::RX, 2 * M_PI);

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

  void testGateRotationYFullRotation() {
    const std::vector<Complex> expected_result_ket_zero{Complex{-1.0, -0.0},
                                                        Complex{0.0, -0.0}};
    const std::vector<Complex> expected_result_ket_one{Complex{0.0, -0.0},
                                                       Complex{-1.0, -0.0}};

    const Matrix unitary = getUnitary(GateID::RY, 2 * M_PI);

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

  void testGateRotationZFullRotation() {
    const std::vector<Complex> expected_result_ket_zero{Complex{-1.0, -0.0},
                                                        Complex{0.0, -0.0}};
    const std::vector<Complex> expected_result_ket_one{Complex{0.0, -0.0},
                                                       Complex{-1.0, -0.0}};

    const Matrix unitary = getUnitary(GateID::RZ, 2 * M_PI);

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

  void testGateRotationXQuarterRotation() {
    const std::vector<Complex> expected_result_ket_zero{
        INV_SQRT_2, Complex{0.0, -1.0} * INV_SQRT_2};
    const std::vector<Complex> expected_result_ket_one{
        Complex{0.0, -1.0} * INV_SQRT_2, INV_SQRT_2};

    const Matrix unitary = getUnitary(GateID::RX, 0.5 * M_PI);

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

  void testGateRotationYQuarterRotation() {
    const std::vector<Complex> expected_result_ket_zero{INV_SQRT_2, INV_SQRT_2};
    const std::vector<Complex> expected_result_ket_one{-INV_SQRT_2, INV_SQRT_2};

    const Matrix unitary = getUnitary(GateID::RY, 0.5 * M_PI);

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

  void testGateRotationZQuarterRotation() {
    const std::vector<Complex> expected_result_ket_zero{
        Complex{cos(0.25 * M_PI), -sin(0.25 * M_PI)}, ZERO};
    const std::vector<Complex> expected_result_ket_one{
        ZERO, Complex{cos(0.25 * M_PI), sin(0.25 * M_PI)}};

    const Matrix unitary = getUnitary(GateID::RZ, 0.5 * M_PI);

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

  void testGateConditionalXInactive() {
    const std::vector<Complex> expected_result_ket_00{ONE, ZERO, ZERO, ZERO};
    const std::vector<Complex> expected_result_ket_01{ZERO, ZERO, ONE, ZERO};

    std::vector<Complex> inactive_ket_00{ONE, ZERO, ZERO, ZERO};
    std::vector<Complex> inactive_ket_01{ZERO, ZERO, ONE, ZERO};

    const Matrix unitary = getUnitary(GateID::X, 0.0);

    const int target_qubit = 0;
    const int control_qubit = 1;
    const int number_qubits = 2;

    applyControlledGate(unitary, control_qubit, target_qubit, number_qubits,
                        inactive_ket_00);
    applyControlledGate(unitary, control_qubit, target_qubit, number_qubits,
                        inactive_ket_01);

    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(inactive_ket_00[i], expected_result_ket_00[i]);
    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(inactive_ket_01[i], expected_result_ket_01[i]);
  }

  void testGateConditionalXActive() {
    const std::vector<Complex> expected_result_ket_10{ZERO, ZERO, ZERO, ONE};
    const std::vector<Complex> expected_result_ket_11{ZERO, ONE, ZERO, ZERO};

    std::vector<Complex> active_ket_10{ZERO, ONE, ZERO, ZERO};
    std::vector<Complex> active_ket_11{ZERO, ZERO, ZERO, ONE};

    const Matrix unitary = getUnitary(GateID::X, 0.0);

    const int target_qubit = 0;
    const int control_qubit = 1;
    const int number_qubits = 2;

    applyControlledGate(unitary, control_qubit, target_qubit, number_qubits,
                        active_ket_10);
    applyControlledGate(unitary, control_qubit, target_qubit, number_qubits,
                        active_ket_11);

    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(active_ket_10[i], expected_result_ket_10[i]);
    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(active_ket_11[i], expected_result_ket_11[i]);
  }

  void testGateConditionalZInactive() {
    const std::vector<Complex> expected_result_ket_00{ONE, ZERO, ZERO, ZERO};
    const std::vector<Complex> expected_result_ket_01{ZERO, ZERO, ONE, ZERO};

    std::vector<Complex> inactive_ket_00{ONE, ZERO, ZERO, ZERO};
    std::vector<Complex> inactive_ket_01{ZERO, ZERO, ONE, ZERO};

    const Matrix unitary = getUnitary(GateID::Z, 0.0);

    const int target_qubit = 0;
    const int control_qubit = 1;
    const int number_qubits = 2;

    applyControlledGate(unitary, control_qubit, target_qubit, number_qubits,
                        inactive_ket_00);
    applyControlledGate(unitary, control_qubit, target_qubit, number_qubits,
                        inactive_ket_01);

    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(inactive_ket_00[i], expected_result_ket_00[i]);
    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(inactive_ket_01[i], expected_result_ket_01[i]);
  }

  void testGateConditionalZActive() {
    const std::vector<Complex> expected_result_ket_10{ZERO, ONE, ZERO, ZERO};
    const std::vector<Complex> expected_result_ket_11{ZERO, ZERO, ZERO, -ONE};

    std::vector<Complex> active_ket_10{ZERO, ONE, ZERO, ZERO};
    std::vector<Complex> active_ket_11{ZERO, ZERO, ZERO, ONE};

    const Matrix unitary = getUnitary(GateID::Z, 0.0);

    const int target_qubit = 0;
    const int control_qubit = 1;
    const int number_qubits = 2;

    applyControlledGate(unitary, control_qubit, target_qubit, number_qubits,
                        active_ket_10);
    applyControlledGate(unitary, control_qubit, target_qubit, number_qubits,
                        active_ket_11);

    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(active_ket_10[i], expected_result_ket_10[i]);
    for (int i = 0; i < 4; i++)
      TS_ASSERT_EQUALS(active_ket_11[i], expected_result_ket_11[i]);
  }

  void testGateDoubleConditionalX() {}

  void testMeasurementKetZero() {}

  void testMeasurementKetOne() {}
};