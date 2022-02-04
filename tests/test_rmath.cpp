#include <gtest/gtest.h>
#include "rmathWrappers.hpp"

TEST(Rmath, pf) {
  ASSERT_NEAR(rmath::pf(1, 1, 1, 1, 0), 0.5, 0.000001);
}