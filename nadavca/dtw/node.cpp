#include <node.h>

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
  for (int i = start_index_; i <= end_index_; i++) {
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
