#ifndef NADAVCA_KMER_MODEL_H
#define NADAVCA_KMER_MODEL_H
#include <functional>
#include <probability.h>
#include <sequence.h>
#include <vector>

class KmerModel {
private:
  int k_;
  int central_position_;
  int alphabet_size_;
  std::vector<double> mean_, additive_constant_, multiplicative_constant_;

  int GetKmerId(const ExtendedSequence *sequence, int index) const;

public:
  KmerModel(int k, int central_position, int alphabet_size,
            std::vector<double> mean, std::vector<double> sigma);
  int GetK() const;
  int GetCentralPosition() const;
  int GetAlphabetSize() const;
  std::vector<double> GetExpectedSignal(const std::vector<int> &reference,
                                        const std::vector<int> &context_before,
                                        const std::vector<int> &context_after);
  std::function<Probability(double)>
  GetDistribution(const ExtendedSequence *sequence, int index) const;
  std::function<Probability(double)>
  GetMixtureDistribution(const ExtendedSequence *sequence, int index1,
                         int index2) const;
};
#endif