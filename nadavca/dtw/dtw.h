#ifndef NADAVCA_DTW_H
#define NADAVCA_DTW_H
#include <kmer_model.h>
#include <vector>

std::vector<std::vector<double>> EstimateLogLikelihoods(
    const std::vector<double> &signal, const std::vector<int> &reference,
    const std::vector<int> &context_before,
    const std::vector<int> &context_after,
    const std::vector<std::vector<int>> &approximate_alignment, int bandwidth,
    const std::vector<double> &short_event_acceptability,
    const KmerModel &kmer_model, bool model_wobbling);

std::vector<std::vector<int>> RefineAlignment(
    const std::vector<double> &signal, const std::vector<int> &reference,
    const std::vector<int> &context_before,
    const std::vector<int> &context_after,
    const std::vector<std::vector<int>> &approximate_alignment, int bandwidth,
    const std::vector<double> &short_event_acceptability, int min_event_length,
    const KmerModel &kmer_model, bool model_transitions);
#endif
