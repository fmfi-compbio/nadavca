#ifndef NADAVCA_NODE_H
#define NADAVCA_NODE_H
#include <functional>
#include <probability.h>
#include <vector>

template <class ScoreJoiningStrategy> class Node {
private:
  int start_index_, end_index_;
  std::vector<Probability> log_likelihoods_;
  Probability &GetReference(int index);

public:
  Node();
  Node(int start_index, int end_index);
  Node(int start_index, int end_index,
       const Node<ScoreJoiningStrategy> &predecessor,
       std::function<Probability(double)> distribution,
       const std::vector<double> &signal, int min_event_length = 0,
       bool reverse = false);
  Probability operator[](int index) const;

  static Probability TotalLikelihood(const Node<ScoreJoiningStrategy> &prefix,
                                     const Node<ScoreJoiningStrategy> &suffix);
};

#include <node_impl.h>
#endif
