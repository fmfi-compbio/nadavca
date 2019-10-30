#include <node.h>

using namespace std;

Node::Node() = default;

Node::Node(int start_index, int end_index, double p)
    : start_index_(start_index), end_index_(end_index),
      likelihoods_(end_index_ + 1 - start_index_, Probability::FromP(p)) {}

Probability &Node::GetReference(int index) {
  return likelihoods_[index - start_index_];
}

Probability Node::operator[](int index) const {
  if (index < start_index_ || index > end_index_)
    return Probability::FromP(0.0);
  return likelihoods_[index - start_index_];
}

int Node::GetStartIndex() const { return start_index_; }

int Node::GetEndIndex() const { return end_index_; }

Node Node::operator*(const Node &other) const {
  Node result(start_index_, end_index_);
  for (int i = start_index_; i <= end_index_; i++) {
    result.GetReference(i) = (*this)[i] * other[i];
  }
  return result;
}

Probability Node::TotalLikelihood(const Node &prefix, const Node &suffix) {
  Probability result = Probability::FromP(0);
  for (int i = prefix.start_index_; i <= prefix.end_index_; i++) {
    result += prefix[i] * suffix[i];
  }
  return result;
}

PathSearchingNode::PathSearchingNode() : Node() {}

PathSearchingNode::PathSearchingNode(const Node &from)
    : Node(from), previous_(end_index_ - start_index_ + 1, -1) {}

int &PathSearchingNode::GetPrevious(int index) {
  return previous_[index - start_index_];
}

int PathSearchingNode::GetBestIndex() const {
  Probability best = Probability::FromP(0.0);
  int result = -1;
  for (int i = start_index_; i < end_index_; i++) {
    if ((*this)[i] > best) {
      best = (*this)[i];
      result = i;
    }
  }
  return result;
}

PathSearchingNode
PathSearchingNode::NextRow(const Node &scores,
                           const PathSearchingNode &predecessor,
                           int min_event_length) {
  PathSearchingNode result(scores);
  int best_previous_index = -1;
  Probability best_previous_score = Probability::FromP(0.0);

  for (int i = predecessor.start_index_;
       i <= predecessor.end_index_ &&
       i < result.start_index_ - min_event_length;
       i++) {
    if (predecessor[i] > best_previous_score) {
      best_previous_score = predecessor[i];
      best_previous_index = i;
    }
  }

  for (int i = result.start_index_; i <= result.end_index_; i++) {
    int index_from = i - min_event_length;
    if (index_from >= predecessor.start_index_ &&
        index_from <= predecessor.end_index_) {
      if (predecessor[index_from] > best_previous_score) {
        best_previous_score = predecessor[index_from];
        best_previous_index = index_from;
      }
    }
    result.GetReference(i) = best_previous_score * scores[i];
    result.GetPrevious(i) = best_previous_index;
  }
  return result;
}

Node Node::NextRowSlow(
    int start_index, int end_index, const Node &predecessor,
    const ExpQuadFunction &mean_distribution,
    const ExpQuadFunction &noise_distribution,
    vector<pair<double, Probability>> mean_distribution_samples,
    const vector<double> &signal, const vector<Probability> &short_event_priors,
    bool reverse) {
  Node result(start_index, end_index, 0);
  int short_threshold = short_event_priors.size() - 1;
  int samples_count = mean_distribution_samples.size();
  if (reverse) {
    vector<Probability> long_event_samples(samples_count,
                                           Probability::FromP(0));

    if (end_index + short_threshold < predecessor.end_index_) {
      ExpQuadFunction current_distribution;
      for (int i = end_index + 1; i <= predecessor.end_index_; i++) {
        current_distribution *=
            noise_distribution.ShiftMaximumPoint(signal[i - 1]);
        if (i - end_index > short_threshold) {
          for (int j = 0; j < samples_count; j++) {
            long_event_samples[j] +=
                current_distribution.Value(mean_distribution_samples[j].first) *
                mean_distribution_samples[j].second * predecessor[i];
          }
        }
      }
    }
    for (int i = end_index; i >= start_index; i--) {
      Probability &current = result.GetReference(i);
      for (const Probability &p : long_event_samples) {
        current += p;
      }
      ExpQuadFunction current_distribution;
      for (int j = i; j <= i + short_threshold && j <= predecessor.end_index_;
           j++) {
        if (j > i) {
          current_distribution *=
              noise_distribution.ShiftMaximumPoint(signal[j - 1]);
        }
        current += (current_distribution * mean_distribution).Integral() *
                   predecessor[j] * short_event_priors[j - i];
      }
      if (i + short_threshold <= predecessor.end_index_ &&
          i - 1 >= start_index) {
        ExpQuadFunction nextSignalpointDistribution =
            noise_distribution.ShiftMaximumPoint(signal[i - 1]);
        for (int j = 0; j < samples_count; j++) {
          long_event_samples[j] +=
              current_distribution.Value(mean_distribution_samples[j].first) *
              mean_distribution_samples[j].second *
              predecessor[i + short_threshold];
          long_event_samples[j] *= nextSignalpointDistribution.Value(
              mean_distribution_samples[j].first);
        }
      }
    }

  } else {
    vector<Probability> long_event_samples(samples_count,
                                           Probability::FromP(0));

    if (start_index - short_threshold > predecessor.start_index_) {
      ExpQuadFunction current_distribution;
      for (int i = start_index - 1; i >= predecessor.start_index_; i--) {
        current_distribution *= noise_distribution.ShiftMaximumPoint(signal[i]);
        if (start_index - i > short_threshold) {
          for (int j = 0; j < samples_count; j++) {
            long_event_samples[j] +=
                current_distribution.Value(mean_distribution_samples[j].first) *
                mean_distribution_samples[j].second * predecessor[i];
          }
        }
      }
    }
    for (int i = start_index; i <= end_index; i++) {
      Probability &current = result.GetReference(i);
      for (const Probability &p : long_event_samples) {
        current += p;
      }
      ExpQuadFunction current_distribution;
      for (int j = i; j >= i - short_threshold && j >= predecessor.start_index_;
           j--) {
        if (j < i) {
          current_distribution *=
              noise_distribution.ShiftMaximumPoint(signal[j]);
        }
        current += (current_distribution * mean_distribution).Integral() *
                   predecessor[j] * short_event_priors[i - j];
      }
      if (i - short_threshold >= predecessor.start_index_ &&
          i + 1 <= end_index) {
        ExpQuadFunction nextSignalpointDistribution =
            noise_distribution.ShiftMaximumPoint(signal[i]);
        for (int j = 0; j < samples_count; j++) {
          long_event_samples[j] +=
              current_distribution.Value(mean_distribution_samples[j].first) *
              mean_distribution_samples[j].second *
              predecessor[i - short_threshold];
          long_event_samples[j] *= nextSignalpointDistribution.Value(
              mean_distribution_samples[j].first);
        }
      }
    }
  }
  return result;
}
