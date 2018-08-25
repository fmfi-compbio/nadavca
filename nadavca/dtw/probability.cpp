#include <probability.h>
#include <algorithm>
#include <cmath>
#include <limits>
using namespace std;

Probability::Probability() : log_val_(-numeric_limits<double>::infinity()) {}

Probability::Probability(double logp) : log_val_(logp) {}

Probability Probability::FromP(double p) { return Probability(log(p)); }

double Probability::Get() { return exp(log_val_); }

double Probability::GetLog() { return log_val_; }

Probability & Probability::operator*=(const Probability & other) {
  log_val_ += other.log_val_;
  return *this;
}

Probability & Probability::operator+=(const Probability & other) {
  log_val_ = (*this + other).log_val_;
  return *this;
}

bool Probability::IsNan() { return log_val_ != log_val_; }

Probability Probability::operator*(Probability other) const {
  return Probability(log_val_ + other.log_val_);
}

Probability Probability::operator+(Probability other) const {
  double a = log_val_, b = other.log_val_;
  if (a < b) swap(a, b);
  if (b == -numeric_limits<double>::infinity()) return Probability(a);
  return Probability(a + log(1 + exp(b - a)));
}

Probability Probability::operator/(Probability other) const {
  return Probability(log_val_ - other.log_val_);
}

bool Probability::operator<(Probability other) const { return log_val_ < other.log_val_; }
