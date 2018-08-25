#include <kmer_model.h>
#include <cmath>
#include <cstdio>
using namespace std;

KmerModel::KmerModel(int k, int central_position, int alphabet_size, vector<double> mean,
                     vector<double> sigma)
    : k_(k), central_position_(central_position), alphabet_size_(alphabet_size), mean_(mean) {
  for (double s : sigma) {
    additive_constant_.push_back(log(1 / sqrt(2 * M_PI * s * s)));
    multiplicative_constant_.push_back(1 / (2 * s * s));
  }
}

int KmerModel::GetK() const { return k_; }

int KmerModel::GetCentralPosition() const { return central_position_; }

int KmerModel::GetAlphabetSize() const { return alphabet_size_; }

int KmerModel::GetKmerId(const ExtendedSequence * sequence, int index) const {
  int result = 0;
  for (int i = index - central_position_; i < index - central_position_ + k_; i++) {
    result *= alphabet_size_;
    result += (*sequence)[i];
  }
  return result;
}

function<Probability(double)> KmerModel::GetDistribution(const ExtendedSequence * sequence,
                                                         int index) const {
  int kmer_id = GetKmerId(sequence, index);
  return [this, kmer_id](double x) {
    double diff = x - mean_[kmer_id];
    return Probability(additive_constant_[kmer_id] -
                       diff * diff * multiplicative_constant_[kmer_id]);
  };
}

function<Probability(double)> KmerModel::GetMixtureDistribution(const ExtendedSequence * sequence,
                                                                int index1, int index2) const {
  auto distribution1 = GetDistribution(sequence, index1);
  auto distribution2 = GetDistribution(sequence, index2);
  return [distribution1, distribution2](double x) {
    return (distribution1(x) + distribution2(x)) / 2;
  };
}