#include <gtest/gtest.h>
#include "rmathWrappers.hpp"
#include <cmath>

TEST(Rmath, pf) {
  ASSERT_NEAR(rmath::pf(1, 1, 1, 1, 0), 0.5, 0.000001);
}

TEST(Rmath, qcauchyl) {
  double a1 = rmath::qcauchy(0.25, 0, 1, 0, 0);
  double b1 = rmath::qcauchyl(0.25L, 0, 1, 0, 0);
  ASSERT_NEAR(a1, b1, 0.00001);

  double a2 = rmath::qcauchy(log(0.25), 0, 1, 0, 1);
  double b2 = rmath::qcauchyl(log(0.25), 0, 1, 0, 1);
  ASSERT_NEAR(a2, b2, 0.0001);

  asm("nop");
}