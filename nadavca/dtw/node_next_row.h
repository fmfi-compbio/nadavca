#ifndef NADAVCA_NODE_NEXT_ROW_H
#define NADAVCA_NODE_NEXT_ROW_H
#include <algorithm>
#include <node.h>

template <class ScoreJoiningStrategy>
Node Node::NextRow(int start_index, int end_index, const Node &predecessor,
                   std::function<Probability(double)> distribution,
                   const std::vector<double> &signal, int min_event_length,
                   bool reverse) {
  Node result(start_index, end_index, 0.0);
  if (reverse) {
    if (end_index + min_event_length <= static_cast<int>(signal.size())) {
      Probability p = Probability::FromP(1);
      for (int i = end_index; i <= predecessor.end_index_; i++) {
        if (i > end_index) {
          p *= distribution(signal[i - 1]);
        }
        if (i - end_index >= min_event_length) {
          ScoreJoiningStrategy::JoinInto(result.GetReference(end_index),
                                         p * predecessor[i]);
        }
      }
    }
    for (int i = std::min(end_index - 1,
                          static_cast<int>(signal.size() - min_event_length));
         i >= start_index; i--) {
      Probability p = Probability::FromP(1);
      for (int j = i; j < i + min_event_length; j++) {
        p *= distribution(signal[j]);
      }
      result.GetReference(i) =
          ScoreJoiningStrategy::Join(p * predecessor[i + min_event_length],
                                     distribution(signal[i]) * result[i + 1]);
    }
  } else {
    if (start_index >= min_event_length) {
      Probability p = Probability::FromP(1);
      for (int i = start_index; i >= predecessor.start_index_; i--) {
        if (i < start_index) {
          p *= distribution(signal[i]);
        }
        if (start_index - i >= min_event_length) {
          ScoreJoiningStrategy::JoinInto(result.GetReference(start_index),
                                         p * predecessor[i]);
        }
      }
    }
    for (int i = std::max(start_index + 1, min_event_length); i <= end_index;
         i++) {
      Probability p = Probability::FromP(1);
      for (int j = i - 1; j >= i - min_event_length; j--) {
        p *= distribution(signal[j]);
      }
      result.GetReference(i) = ScoreJoiningStrategy::Join(
          p * predecessor[i - min_event_length],
          distribution(signal[i - 1]) * result[i - 1]);
    }
  }
  return result;
}

#endif
