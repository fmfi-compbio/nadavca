#include <joining_strategy.h>

Probability JoiningStrategySum::Join(const Probability &a,
                                     const Probability &b) {
  return a + b;
}

void JoiningStrategySum::JoinInto(Probability &a, const Probability &b) {
  a += b;
}

Probability JoiningStrategyMax::Join(const Probability &a,
                                     const Probability &b) {
  return a > b ? a : b;
}

void JoiningStrategyMax::JoinInto(Probability &a, const Probability &b) {
  if (b > a)
    a = b;
}
