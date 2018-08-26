#ifndef NADAVCA_PROBABILITY_H
#define NADAVCA_PROBABILITY_H

class Probability {
 private:
  double log_val_;

 public:
  Probability();
  Probability(double logp);
  static Probability FromP(double p);

  double Get();
  bool IsNan();
  double GetLog();

  Probability & operator*=(const Probability & other);
  Probability & operator+=(const Probability & other);
  Probability operator*(Probability other) const;
  Probability operator+(Probability other) const;
  Probability operator/(Probability other) const;
  bool operator<(Probability other) const;
  bool operator>(Probability other) const;
};
#endif