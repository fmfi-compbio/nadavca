#ifndef NADAVCA_NODE_H
#define NADAVCA_NODE_H
#include <probability.h>
#include <functional>
#include <vector>

class Node {
 private:
  int start_index_, end_index_;
  std::vector<Probability> log_likelihoods_;
  Probability & GetReference(int index);

 public:
  Node();
  Node(int start_index, int end_index);
  Node(int start_index, int end_index, const Node & predecessor,
       std::function<Probability(double)> distribution, const std::vector<double> & signal,
       int min_event_length = 0, bool reverse = false);
  Probability operator[](int index) const;
  friend Probability TotalLikelihood(const Node & prefix, const Node & suffix);
  friend int MostLikelyIndex(const Node & prefix, const Node & suffix);
};

Probability TotalLikelihood(const Node & prefix, const Node & suffix);
int MostLikelyIndex(const Node & prefix, const Node & suffix);
#endif