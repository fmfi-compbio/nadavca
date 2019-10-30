#ifndef NADAVCA_EXPQUAD_FUNCTION_H
#define NADAVCA_EXPQUAD_FUNCTION_H

#include <probability.h>

class ExpQuadFunction {
private:
  double coefficients_[3];

public:
  ExpQuadFunction(double a = 0, double b = 0, double c = 0);
  Probability Value(double x);
  Probability Integral();

  ExpQuadFunction ShiftMaximumPoint(double newMaximumPoint) const;
  ExpQuadFunction &operator*=(const ExpQuadFunction &other);
  ExpQuadFunction operator*(const ExpQuadFunction &other) const;

  void DebugPrint();
};

#endif
