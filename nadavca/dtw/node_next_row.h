#ifndef NADAVCA_NODE_NEXT_ROW_H
#define NADAVCA_NODE_NEXT_ROW_H
#include <algorithm>
#include <node.h>

template <class ScoreJoiningStrategy>
Node Node::NextRow(int start_index, int end_index, const Node &predecessor,
                   std::function<Probability(double)> distribution,
                   const std::vector<double> &signal,
                   const std::vector<Probability> &short_event_acceptability,
                   bool reverse) {
  Node result(start_index, end_index, 0.0);
  int short_threshold = short_event_acceptability.size() - 1;

  if (reverse) {
    Probability long_events = Probability::FromP(0);
    if (end_index + short_threshold < predecessor.end_index_) {
      Probability p = Probability::FromP(1);
      for (int i = end_index + 1; i <= predecessor.end_index_; i++) {
        p *= distribution(signal[i - 1]);
        if (i - end_index > short_threshold) {
          ScoreJoiningStrategy::JoinInto(long_events, p * predecessor[i]);
        }
      }
    }
    for (int i = end_index; i >= start_index; i--) {
      if (i < end_index) {
        long_events *= distribution(signal[i]);
      }
      Probability &current = result.GetReference(i);
      current = long_events;
      Probability p = Probability::FromP(1);
      for (int j = i; j <= i + short_threshold && j <= predecessor.end_index_;
           j++) {
        if (j > i) {
          p *= distribution(signal[j - 1]);
        }
        ScoreJoiningStrategy::JoinInto(
            current, p * predecessor[j] * short_event_acceptability[j - i]);
      }
      if (i + short_threshold <= predecessor.end_index_) {
        ScoreJoiningStrategy::JoinInto(long_events,
                                       p * predecessor[i + short_threshold]);
      }
    }
  } else {
    Probability long_events = Probability::FromP(0);
    if (start_index - short_threshold > predecessor.start_index_) {
      Probability p = Probability::FromP(1);
      for (int i = start_index - 1; i >= predecessor.start_index_; i--) {
        p *= distribution(signal[i]);
        if (start_index - i > short_threshold) {
          ScoreJoiningStrategy::JoinInto(long_events, p * predecessor[i]);
        }
      }
    }
    for (int i = start_index; i <= end_index; i++) {
      if (i > start_index) {
        long_events *= distribution(signal[i - 1]);
      }
      Probability &current = result.GetReference(i);
      current = long_events;
      Probability p = Probability::FromP(1);
      for (int j = i; j >= i - short_threshold && j >= predecessor.start_index_;
           j--) {
        if (j < i) {
          p *= distribution(signal[j]);
        }
        ScoreJoiningStrategy::JoinInto(
            current, p * predecessor[j] * short_event_acceptability[i - j]);
      }
      if (i - short_threshold >= predecessor.start_index_) {
        ScoreJoiningStrategy::JoinInto(long_events,
                                       p * predecessor[i - short_threshold]);
      }
    }
  }
  return result;
}

#endif
