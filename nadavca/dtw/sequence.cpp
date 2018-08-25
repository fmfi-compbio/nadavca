#include <sequence.h>
using namespace std;

ExtendedSequence::ExtendedSequence() = default;

ExtendedSequence::ExtendedSequence(const vector<int> & main_sequence,
                                   const vector<int> & context_before,
                                   const vector<int> & context_after) {
  offset = context_before.size();
  values.reserve(main_sequence.size() + context_before.size() + context_after.size());
  for (int base : context_before) {
    values.push_back(base);
  }
  for (int base : main_sequence) {
    values.push_back(base);
  }
  for (int base : context_after) {
    values.push_back(base);
  }
}

int ExtendedSequence::operator[](int index) const {
  index += offset;
  if (index < 0 || index >= static_cast<int>(values.size())) return 0;
  return values[index];
}

ModifiedSequence::ModifiedSequence(const ExtendedSequence * original, int index, int new_value)
    : wrapped_(original), changed_index_(index), changed_value_(new_value) {}

int ModifiedSequence::operator[](int index) const {
  if (index == changed_index_) return changed_value_;
  return (*wrapped_)[index];
}
