#include <dtw.h>
#include <node.h>
#include <sequence.h>
using namespace std;

vector<int> ComputeBandStarts(const vector<vector<int>> & approximate_alignment, int read_length,
                              int reference_length, int bandwidth) {
  vector<int> result(reference_length + 1, 0);
  for (const vector<int> & point : approximate_alignment) {
    int signal_index = point[0];
    int reference_index = point[1];
    result[reference_index] = max(0, signal_index - bandwidth);
  }
  for (int i = 1; i <= reference_length; i++) {
    result[i] = max(result[i], result[i - 1]);
  }
  return result;
}

vector<int> ComputeBandEnds(const vector<vector<int>> & approximate_alignment, int read_length,
                            int reference_length, int bandwidth) {
  vector<int> result(reference_length + 1, read_length);
  for (const vector<int> & point : approximate_alignment) {
    int signal_index = point[0];
    int reference_index = point[1];
    result[reference_index] = min(read_length, signal_index + bandwidth);
  }
  for (int i = reference_length - 1; i >= 0; i--) {
    result[i] = min(result[i], result[i + 1]);
  }
  return result;
}

vector<vector<double>> EstimateLogLikelihoods(const vector<double> & signal,
                                              const vector<int> & reference,
                                              const vector<int> & context_before,
                                              const vector<int> & context_after,
                                              const vector<vector<int>> & approximate_alignment,
                                              int bandwidth, int min_event_length,
                                              const KmerModel & kmer_model, bool model_wobbling) {
  ExtendedSequence sequence(reference, context_before, context_after);
  vector<int> band_starts =
      ComputeBandStarts(approximate_alignment, signal.size(), reference.size(), bandwidth);
  vector<int> band_ends =
      ComputeBandEnds(approximate_alignment, signal.size(), reference.size(), bandwidth);

  vector<Node> prefix_likelihoods(reference.size() + 1), suffix_likelihoods(reference.size() + 1);
  prefix_likelihoods[0] = Node(band_starts[0], band_ends[0]);
  for (unsigned i = 0; i < reference.size(); i++) {
    Node predecessor = prefix_likelihoods[i];
    if (i > 0 && model_wobbling) {
      auto distribution = kmer_model.GetMixtureDistribution(&sequence, i - 1, i);
      predecessor = Node(band_starts[i], band_ends[i], predecessor, distribution, signal, 0, false);
    }
    auto distribution = kmer_model.GetDistribution(&sequence, i);
    prefix_likelihoods[i + 1] = Node(band_starts[i + 1], band_ends[i + 1], predecessor,
                                     distribution, signal, min_event_length, false);
  }

  suffix_likelihoods[reference.size()] =
      Node(band_starts[reference.size()], band_ends[reference.size()]);
  for (unsigned i = reference.size(); i > 0; i--) {
    Node predecessor = suffix_likelihoods[i];
    if (i < reference.size() && model_wobbling) {
      auto distribution = kmer_model.GetMixtureDistribution(&sequence, i, i - 1);
      predecessor = Node(band_starts[i], band_ends[i], predecessor, distribution, signal, 0, true);
    }
    auto distribution = kmer_model.GetDistribution(&sequence, i - 1);
    suffix_likelihoods[i - 1] = Node(band_starts[i - 1], band_ends[i - 1], predecessor,
                                     distribution, signal, min_event_length, true);
  }

  Probability no_snp_likelihood =
      TotalLikelihood(prefix_likelihoods[reference.size()], suffix_likelihoods[reference.size()]);

  int alphabet_size = kmer_model.GetAlphabetSize();
  int influence_back = kmer_model.GetK() - kmer_model.GetCentralPosition() - 1;
  int influence_forward = kmer_model.GetCentralPosition();
  vector<vector<double>> result(reference.size(), vector<double>(alphabet_size));

  for (int i = 0; i < static_cast<int>(reference.size()); i++) {
    int first_influenced = max(0, i - influence_back);
    int last_influenced = min(static_cast<int>(reference.size() - 1), i + influence_forward);
    for (int base = 0; base < alphabet_size; base++) {
      if (base == sequence[i]) {
        result[i][base] = no_snp_likelihood.GetLog();
      } else {
        ModifiedSequence modified(&sequence, i, base);
        Node current = prefix_likelihoods[first_influenced];
        for (int j = first_influenced; j <= last_influenced; j++) {
          if (j > 0 && model_wobbling) {
            auto distribution = kmer_model.GetMixtureDistribution(&modified, j - 1, j);
            current = Node(band_starts[j], band_ends[j], current, distribution, signal, 0, false);
          }
          auto distribution = kmer_model.GetDistribution(&modified, j);
          current = Node(band_starts[j + 1], band_ends[j + 1], current, distribution, signal,
                         min_event_length, false);
        }
        if (last_influenced + 1 < static_cast<int>(reference.size()) && model_wobbling) {
          auto distribution =
              kmer_model.GetMixtureDistribution(&modified, last_influenced, last_influenced + 1);
          current = Node(band_starts[last_influenced], band_ends[last_influenced], current,
                         distribution, signal, 0, false);
        }
        result[i][base] =
            TotalLikelihood(current, suffix_likelihoods[last_influenced + 1]).GetLog();
      }
    }
  }
  return result;
}