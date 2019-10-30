#include <cmath>
#include <cstdio>
#include <expquad_function.h>

ExpQuadFunction::ExpQuadFunction(double a, double b, double c) {
  coefficients_[0] = c;
  coefficients_[1] = b;
  coefficients_[2] = a;
}

Probability ExpQuadFunction::Value(double x) {
  return Probability(coefficients_[2] * x * x + coefficients_[1] * x +
                     coefficients_[0]);
}

Probability ExpQuadFunction::Integral() {
  double a = coefficients_[2];
  double b = coefficients_[1];
  double c = coefficients_[0];
  return Probability(c - b * b / (4 * a) + (log(M_PI) - log(-a)) / 2);
}

ExpQuadFunction
ExpQuadFunction::ShiftMaximumPoint(double newMaximumPoint) const {
  double a = coefficients_[2];
  double b = coefficients_[1];
  double c = coefficients_[0];
  return ExpQuadFunction(a, -2 * a * newMaximumPoint,
                         c - b * b / (4 * a) +
                             a * newMaximumPoint * newMaximumPoint);
}

ExpQuadFunction &ExpQuadFunction::operator*=(const ExpQuadFunction &other) {
  for (int i = 0; i < 3; i++) {
    coefficients_[i] += other.coefficients_[i];
  }
  return *this;
}

ExpQuadFunction ExpQuadFunction::operator*(const ExpQuadFunction &other) const {
  ExpQuadFunction result;
  for (int i = 0; i < 3; i++) {
    result.coefficients_[i] = coefficients_[i] + other.coefficients_[i];
  }
  return result;
}

void ExpQuadFunction::DebugPrint() {
  printf("exp(%lf x^2 + %lf x + %lf)\n", coefficients_[0], coefficients_[1],
         coefficients_[2]);
}
