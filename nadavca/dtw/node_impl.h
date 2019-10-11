#ifndef NADAVCA_NODE_IMPL_H
#define NADAVCA_NODE_IMPL_H
#include <algorithm>
#include <node.h>
using namespace std;

template <class ScoreJoiningStrategy>
Node<ScoreJoiningStrategy>::Node() = default;

template <class ScoreJoiningStrategy>
Node<ScoreJoiningStrategy>::Node(int start_index, int end_index)
    : start_index_(start_index), end_index_(end_index),
      log_likelihoods_(end_index_ + 1 - start_index_, Probability::FromP(1.0)) {
}

template <class ScoreJoiningStrategy>
Node<ScoreJoiningStrategy>::Node(
    int start_index, int end_index,
    const Node<ScoreJoiningStrategy> &predecessor,
    std::function<Probability(double)> distribution,
    const std::vector<double> &signal, int min_event_length, bool reverse)
    : start_index_(start_index), end_index_(end_index),
      log_likelihoods_(end_index_ + 1 - start_index_, Probability::FromP(0.0)) {
  if (reverse) {
    if (end_index_ + min_event_length <= static_cast<int>(signal.size())) {
      Probability p = Probability::FromP(1);
      for (int i = end_index_; i <= predecessor.end_index_; i++) {
        if (i > end_index_) {
          p *= distribution(signal[i - 1]);
        }
        if (i - end_index_ >= min_event_length) {
          ScoreJoiningStrategy::JoinInto(GetReference(end_index_),
                                         p * predecessor[i]);
        }
      }
    }
    for (int i = min(end_index_ - 1,
                     static_cast<int>(signal.size() - min_event_length));
         i >= start_index_; i--) {
      Probability p = Probability::FromP(1);
      for (int j = i; j < i + min_event_length; j++) {
        p *= distribution(signal[j]);
      }
      GetReference(i) =
          ScoreJoiningStrategy::Join(p * predecessor[i + min_event_length],
                                     distribution(signal[i]) * (*this)[i + 1]);
    }
  } else {
    if (start_index_ >= min_event_length) {
      Probability p = Probability::FromP(1);
      for (int i = start_index_; i >= predecessor.start_index_; i--) {
        if (i < start_index_) {
          p *= distribution(signal[i]);
        }
        if (start_index_ - i >= min_event_length) {
          ScoreJoiningStrategy::JoinInto(GetReference(start_index_),
                                         p * predecessor[i]);
        }
      }
    }
    for (int i = max(start_index_ + 1, min_event_length); i <= end_index_;
         i++) {
      Probability p = Probability::FromP(1);
      for (int j = i - 1; j >= i - min_event_length; j--) {
        p *= distribution(signal[j]);
      }
      GetReference(i) = ScoreJoiningStrategy::Join(
          p * predecessor[i - min_event_length],
          distribution(signal[i - 1]) * (*this)[i - 1]);
    }
  }
}

template <class ScoreJoiningStrategy>
Probability &Node<ScoreJoiningStrategy>::GetReference(int index) {
  return log_likelihoods_[index - start_index_];
}

template <class ScoreJoiningStrategy>
Probability Node<ScoreJoiningStrategy>::operator[](int index) const {
  if (index < start_index_ || index > end_index_)
    return Probability::FromP(0.0);
  return log_likelihoods_[index - start_index_];
}

template <class ScoreJoiningStrategy>
Probability Node<ScoreJoiningStrategy>::TotalLikelihood(
    const Node<ScoreJoiningStrategy> &prefix,
    const Node<ScoreJoiningStrategy> &suffix) {
  Probability result = Probability::FromP(0);
  for (int i = prefix.start_index_; i <= prefix.end_index_; i++) {
    result += prefix[i] * suffix[i];
  }
  return result;
}

template <class ScoreJoiningStrategy>
int Node<ScoreJoiningStrategy>::MostLikelyIndex(
    const Node<ScoreJoiningStrategy> &prefix,
    const Node<ScoreJoiningStrategy> &suffix) {
  Probability max_likelihood = Probability::FromP(0);
  int best_index = -1;
  for (int i = prefix.start_index_; i <= prefix.end_index_; i++) {
    Probability likelihood = prefix[i] * suffix[i];
    if (likelihood > max_likelihood) {
      max_likelihood = likelihood;
      best_index = i;
    }
  }
  return best_index;
}
#endif
