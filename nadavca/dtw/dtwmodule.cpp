#include <cstdio>
#include <dtw.h>
#include <kmer_model.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(dtw, m) {
  m.doc() = "A C++ extension for fast DTW computation";
  py::class_<KmerModel>(m, "KmerModel")
      .def(py::init<int, int, int, vector<double>, vector<double>>())
      .def("get_k", &KmerModel::GetK)
      .def("get_central_position", &KmerModel::GetCentralPosition)
      .def("get_expected_signal", &KmerModel::GetExpectedSignal, "",
           py::arg("reference"), py::arg("context_before"),
           py::arg("context_after"));
  m.def("estimate_log_likelihoods", &EstimateLogLikelihoods, "",
        py::arg("signal"), py::arg("reference"), py::arg("context_before"),
        py::arg("context_after"), py::arg("approximate_alignment"),
        py::arg("bandwidth"), py::arg("min_event_length"),
        py::arg("kmer_model"), py::arg("model_wobbling"));
  m.def("refine_alignment", &RefineAlignment, "", py::arg("signal"),
        py::arg("reference"), py::arg("context_before"),
        py::arg("context_after"), py::arg("approximate_alignment"),
        py::arg("bandwidth"), py::arg("min_event_length"),
        py::arg("kmer_model"));
}