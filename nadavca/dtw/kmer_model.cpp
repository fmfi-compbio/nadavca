#include <cmath>
#include <cstdio>
#include <kmer_model.h>
using namespace std;

KmerModel::KmerModel(int k, int central_position, int alphabet_size,
                     vector<double> mean, vector<double> sigma)
    : k_(k), central_position_(central_position), alphabet_size_(alphabet_size),
      mean_(mean) {
  for (double s : sigma) {
    additive_constant_.push_back(log(1 / sqrt(2 * M_PI * s * s)));
    multiplicative_constant_.push_back(1 / (2 * s * s));
  }
}

int KmerModel::GetK() const { return k_; }

int KmerModel::GetCentralPosition() const { return central_position_; }

int KmerModel::GetAlphabetSize() const { return alphabet_size_; }

int KmerModel::GetKmerId(const ExtendedSequence *sequence, int index) const {
  int result = 0;
  for (int i = index - central_position_; i < index - central_position_ + k_;
       i++) {
    result *= alphabet_size_;
    result += (*sequence)[i];
  }
  return result;
}

vector<double> KmerModel::GetExpectedSignal(const vector<int> &reference,
                                            const vector<int> &context_before,
                                            const vector<int> &context_after) {
  ExtendedSequence sequence(reference, context_before, context_after);
  vector<double> result(reference.size());
  for (unsigned i = 0; i < reference.size(); i++) {
    int kmer_id = GetKmerId(&sequence, i);
    result[i] = mean_[kmer_id];
  }
  return result;
}

function<Probability(double)>
KmerModel::GetDistribution(const ExtendedSequence *sequence, int index) const {
  int kmer_id = GetKmerId(sequence, index);
  return [this, kmer_id](double x) {
    double diff = x - mean_[kmer_id];
    return Probability(additive_constant_[kmer_id] -
                       diff * diff * multiplicative_constant_[kmer_id]);
  };
}

function<Probability(double)>
KmerModel::GetMixtureDistribution(const ExtendedSequence *sequence, int index1,
                                  int index2) const {
  auto distribution1 = GetDistribution(sequence, index1);
  auto distribution2 = GetDistribution(sequence, index2);
  return [distribution1, distribution2](double x) {
    return (distribution1(x) + distribution2(x)) / 2;
  };
}

function<Probability(double)>
KmerModel::GetTransitionDistribution(const ExtendedSequence *sequence,
                                     int index1, int index2) const {
  int kmer_id1 = GetKmerId(sequence, index1);
  int kmer_id2 = GetKmerId(sequence, index2);
  double mean1 = mean_[kmer_id1];
  double mean2 = mean_[kmer_id2];

  Probability p_out = Probability::FromP(0.0);
  if (mean1 == mean2) {
    return [p_out](double x) { return p_out; };
  }
  if (mean1 > mean2)
    swap(mean1, mean2);

  Probability p_in = Probability::FromP(0.01);

  double mean = (mean1 + mean2) / 2.0;
  double sigma = (mean2 - mean1) / 3.0;
  double ac = log(1 / sqrt(2 * M_PI * sigma * sigma) * 0.1);
  double mc = 1 / (2 * sigma * sigma);
  return [mean1, mean2, p_in, p_out, mean, sigma, ac, mc](double x) {
    if (x < mean1 || x > mean2)
      return p_out;
    //    return p_in;

    double diff = x - mean;
    return Probability(ac - diff * diff * mc);
  };
}
