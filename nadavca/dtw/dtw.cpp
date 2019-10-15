#include <dtw.h>
#include <joining_strategy.h>
#include <node.h>
#include <sequence.h>
using namespace std;

typedef Node<JoiningStrategySum> SumNode;
typedef Node<JoiningStrategyMax> MaxNode;

vector<int> ComputeBandStarts(const vector<vector<int>> &approximate_alignment,
                              int read_length, int reference_length,
                              int bandwidth) {
  vector<int> result(reference_length + 1, 0);
  for (const vector<int> &point : approximate_alignment) {
    int signal_index = point[0];
    int reference_index = point[1];
    result[reference_index] = max(0, signal_index - bandwidth);
  }
  for (int i = 1; i <= reference_length; i++) {
    result[i] = max(result[i], result[i - 1]);
  }
  return result;
}

vector<int> ComputeBandEnds(const vector<vector<int>> &approximate_alignment,
                            int read_length, int reference_length,
                            int bandwidth) {
  vector<int> result(reference_length + 1, read_length);
  for (const vector<int> &point : approximate_alignment) {
    int signal_index = point[0];
    int reference_index = point[1];
    result[reference_index] = min(read_length, signal_index + bandwidth);
  }
  for (int i = reference_length - 1; i >= 0; i--) {
    result[i] = min(result[i], result[i + 1]);
  }
  return result;
}

vector<vector<double>> EstimateLogLikelihoods(
    const vector<double> &signal, const vector<int> &reference,
    const vector<int> &context_before, const vector<int> &context_after,
    const vector<vector<int>> &approximate_alignment, int bandwidth,
    int min_event_length, const KmerModel &kmer_model, bool model_wobbling) {
  ExtendedSequence sequence(reference, context_before, context_after);
  vector<int> band_starts = ComputeBandStarts(
      approximate_alignment, signal.size(), reference.size(), bandwidth);
  vector<int> band_ends = ComputeBandEnds(approximate_alignment, signal.size(),
                                          reference.size(), bandwidth);

  vector<SumNode> prefix_likelihoods(reference.size() + 1),
      suffix_likelihoods(reference.size() + 1);
  prefix_likelihoods[0] = SumNode(band_starts[0], band_ends[0]);
  for (unsigned i = 0; i < reference.size(); i++) {
    SumNode predecessor = prefix_likelihoods[i];
    if (i > 0 && model_wobbling) {
      auto distribution =
          kmer_model.GetMixtureDistribution(&sequence, i - 1, i);
      predecessor = SumNode(band_starts[i], band_ends[i], predecessor,
                            distribution, signal, 0, false);
    }
    auto distribution = kmer_model.GetDistribution(&sequence, i);
    prefix_likelihoods[i + 1] =
        SumNode(band_starts[i + 1], band_ends[i + 1], predecessor, distribution,
                signal, min_event_length, false);
  }

  suffix_likelihoods[reference.size()] =
      SumNode(band_starts[reference.size()], band_ends[reference.size()]);
  for (unsigned i = reference.size(); i > 0; i--) {
    SumNode predecessor = suffix_likelihoods[i];
    if (i < reference.size() && model_wobbling) {
      auto distribution =
          kmer_model.GetMixtureDistribution(&sequence, i, i - 1);
      predecessor = SumNode(band_starts[i], band_ends[i], predecessor,
                            distribution, signal, 0, true);
    }
    auto distribution = kmer_model.GetDistribution(&sequence, i - 1);
    suffix_likelihoods[i - 1] =
        SumNode(band_starts[i - 1], band_ends[i - 1], predecessor, distribution,
                signal, min_event_length, true);
  }

  Probability no_snp_likelihood =
      SumNode::TotalLikelihood(prefix_likelihoods[reference.size()],
                               suffix_likelihoods[reference.size()]);

  int alphabet_size = kmer_model.GetAlphabetSize();
  int influence_back = kmer_model.GetK() - kmer_model.GetCentralPosition() - 1;
  int influence_forward = kmer_model.GetCentralPosition();
  vector<vector<double>> result(reference.size(),
                                vector<double>(alphabet_size));

  for (int i = 0; i < static_cast<int>(reference.size()); i++) {
    int first_influenced = max(0, i - influence_back);
    int last_influenced =
        min(static_cast<int>(reference.size() - 1), i + influence_forward);
    for (int base = 0; base < alphabet_size; base++) {
      if (base == sequence[i]) {
        result[i][base] = no_snp_likelihood.GetLog();
      } else {
        ModifiedSequence modified(&sequence, i, base);
        SumNode current = prefix_likelihoods[first_influenced];
        for (int j = first_influenced; j <= last_influenced; j++) {
          if (j > 0 && model_wobbling) {
            auto distribution =
                kmer_model.GetMixtureDistribution(&modified, j - 1, j);
            current = SumNode(band_starts[j], band_ends[j], current,
                              distribution, signal, 0, false);
          }
          auto distribution = kmer_model.GetDistribution(&modified, j);
          current = SumNode(band_starts[j + 1], band_ends[j + 1], current,
                            distribution, signal, min_event_length, false);
        }
        if (last_influenced + 1 < static_cast<int>(reference.size()) &&
            model_wobbling) {
          auto distribution = kmer_model.GetMixtureDistribution(
              &modified, last_influenced, last_influenced + 1);
          current =
              SumNode(band_starts[last_influenced], band_ends[last_influenced],
                      current, distribution, signal, 0, false);
        }
        result[i][base] = SumNode::TotalLikelihood(
                              current, suffix_likelihoods[last_influenced + 1])
                              .GetLog();
      }
    }
  }
  return result;
}

vector<int> RefineAlignment(const vector<double> &signal,
                            const vector<int> &reference,
                            const vector<int> &context_before,
                            const vector<int> &context_after,
                            const vector<vector<int>> &approximate_alignment,
                            int bandwidth, int min_event_length,
                            const KmerModel &kmer_model) {
  ExtendedSequence sequence(reference, context_before, context_after);
  vector<int> band_starts = ComputeBandStarts(
      approximate_alignment, signal.size(), reference.size(), bandwidth);
  vector<int> band_ends = ComputeBandEnds(approximate_alignment, signal.size(),
                                          reference.size(), bandwidth);

  vector<SumNode> prefix_likelihoods(reference.size() + 1),
      suffix_likelihoods(reference.size() + 1);

  prefix_likelihoods[0] = SumNode(band_starts[0], band_ends[0]);
  for (unsigned i = 0; i < reference.size(); i++) {
    SumNode predecessor = prefix_likelihoods[i];
    auto distribution = kmer_model.GetDistribution(&sequence, i);
    prefix_likelihoods[i + 1] =
        SumNode(band_starts[i + 1], band_ends[i + 1], predecessor, distribution,
                signal, min_event_length, false);
  }

  suffix_likelihoods[reference.size()] =
      SumNode(band_starts[reference.size()], band_ends[reference.size()]);
  for (unsigned i = reference.size(); i > 0; i--) {
    SumNode predecessor = suffix_likelihoods[i];
    auto distribution = kmer_model.GetDistribution(&sequence, i - 1);
    suffix_likelihoods[i - 1] =
        SumNode(band_starts[i - 1], band_ends[i - 1], predecessor, distribution,
                signal, min_event_length, true);
  }

  vector<vector<Probability>> all_paths_sum(reference.size() + 1);
  for (unsigned i = 0; i <= reference.size(); i++) {
    all_paths_sum[i].resize(band_ends[i] - band_starts[i] + 1);
    for (int j = band_starts[i]; j <= band_ends[i]; j++) {
      all_paths_sum[i][j - band_starts[i]] =
          prefix_likelihoods[i][j] * suffix_likelihoods[i][j];
    }
  }

  vector<vector<Probability>> dp(reference.size() + 1);
  vector<vector<int>> previous(reference.size() + 1);
  dp[0] = all_paths_sum[0];
  previous[0].resize(band_ends[0] - band_starts[0] + 1, -1);
  for (unsigned i = 1; i <= reference.size(); i++) {
    dp[i].resize(band_ends[i] - band_starts[i] + 1, Probability::FromP(0.0));
    previous[i].resize(band_ends[i] - band_starts[i] + 1, -1);
    int best_previous_index = -1;
    Probability best_previous_score = Probability::FromP(0.0);

    for (int j = band_starts[i - 1];
         j <= band_ends[i - 1] && j < band_starts[i] - min_event_length; j++) {
      Probability previous_score = dp[i - 1][j - band_starts[i - 1]];
      if (previous_score > best_previous_score) {
        best_previous_score = previous_score;
        best_previous_index = j;
      }
    }

    for (int j = band_starts[i]; j <= band_ends[i]; j++) {
      if (j - min_event_length >= band_starts[i - 1] &&
          j - min_event_length <= band_ends[i - 1]) {
        Probability previous_score =
            dp[i - 1][j - min_event_length - band_starts[i - 1]];
        if (previous_score > best_previous_score) {
          best_previous_score = previous_score;
          best_previous_index = j - min_event_length;
        }
      }
      dp[i][j - band_starts[i]] =
          best_previous_score * all_paths_sum[i][j - band_starts[i]];
      previous[i][j - band_starts[i]] = best_previous_index;
    }
  }

  vector<int> best_signal_index(reference.size() + 1);
  Probability best_score = Probability::FromP(0.0);
  int best_index = -1;
  for (int i = band_starts[reference.size()]; i <= band_ends[reference.size()];
       i++) {
    Probability score = dp[reference.size()][i - band_starts[reference.size()]];
    if (score > best_score) {
      best_score = score;
      best_index = i;
    }
  }
  for (int i = reference.size(); i >= 0; i--) {
    best_signal_index[i] = best_index;
    best_index = previous[i][best_index - band_starts[i]];
  }
  return best_signal_index;
}
