#ifndef NADAVCA_JOINING_STRATEGY_H
#define NADAVCA_JOINING_STRATEGY_H
#include <probability.h>

class JoiningStrategySum {
public:
  static Probability Join(const Probability &a, const Probability &b);
  static void JoinInto(Probability &a, const Probability &b);
};

class JoiningStrategyMax {
public:
  static Probability Join(const Probability &a, const Probability &b);
  static void JoinInto(Probability &a, const Probability &b);
};

#endif
