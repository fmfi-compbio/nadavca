#ifndef NADAVCA_NODE_H
#define NADAVCA_NODE_H
#include <functional>
#include <probability.h>
#include <vector>

class Node {
protected:
  int start_index_, end_index_;
  std::vector<Probability> likelihoods_;

public:
  Node();
  Node(int start_index, int end_index, double p = 1.0);
  Probability &GetReference(int index);
  Probability operator[](int index) const;

  int GetStartIndex() const;
  int GetEndIndex() const;

  Node operator*(const Node &other) const;

  static Probability TotalLikelihood(const Node &prefix, const Node &suffix);

  template <class ScoreJoiningStrategy>
  static Node NextRow(int start_index, int end_index, const Node &predecessor,
                      std::function<Probability(double)> distribution,
                      const std::vector<double> &signal,
                      const std::vector<Probability> &short_event_acceptability,
                      bool reverse = false);
  // short_event_acceptability must contain at least one element
};

class PathSearchingNode : public Node {
private:
  std::vector<int> previous_;

public:
  PathSearchingNode();
  PathSearchingNode(const Node &from);

  int &GetPrevious(int index);

  int GetBestIndex() const;
  static PathSearchingNode NextRow(const Node &values,
                                   const PathSearchingNode &predecessor,
                                   int min_event_length);
};

#include <node_next_row.h>
#endif
