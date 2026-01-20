PYTHON ?= .venv/bin/python
CXX ?= c++

EXT_SUFFIX := $(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '')")
PYBIND11_INCLUDES := $(shell $(PYTHON) -m pybind11 --includes 2>/dev/null)

CXXFLAGS ?= -O3 -std=c++17 -fPIC
LDFLAGS ?= -undefined dynamic_lookup

CORE_SOURCES := cpp/core/entanglement.cpp cpp/core/geometry2d.cpp cpp/core/geometry3d.cpp cpp/core/surface_builder.cpp cpp/core/pseudoknot_decomposition.cpp cpp/bindings/pybind_module.cpp
CORE_TARGET := python/rnaknotdetector_core$(EXT_SUFFIX)

.PHONY: all clean

all: $(CORE_TARGET)

$(CORE_TARGET): $(CORE_SOURCES)
	$(CXX) $(CXXFLAGS) $(PYBIND11_INCLUDES) -Icpp/core $^ -shared -o $@ $(LDFLAGS)

clean:
	rm -f $(CORE_TARGET)
