#ifndef NADAVCA_SEQUENCE_H
#define NADAVCA_SEQUENCE_H
#include <vector>

class ExtendedSequence {
private:
  std::vector<int> values;
  int offset;

public:
  ExtendedSequence();
  ExtendedSequence(const std::vector<int> &main_sequence,
                   const std::vector<int> &context_before,
                   const std::vector<int> &context_after);
  virtual int operator[](int index) const;
};

class ModifiedSequence : public ExtendedSequence {
private:
  const ExtendedSequence *wrapped_;
  int changed_index_;
  int changed_value_;

public:
  ModifiedSequence(const ExtendedSequence *original, int index, int new_value);
  int operator[](int index) const override;
};

#endif