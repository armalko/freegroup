#include "magnus.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_magnus, m) {
    m.doc() = "Magnus normal‐closure bindings";

    m.def("magnus_is_from_normal_closure_cpp",
          // two‐arg wrapper
          [](const std::vector<int>& w,
             const std::vector<int>& r) {
              return MagnusIsFromNormalClosure(w, r);
          },
          py::arg("word"), py::arg("relator"),
          "Decide if `word` lies in the normal closure of `relator`");
}
