PYTHON ?= python3
CXX ?= c++

EXT_SUFFIX := $(shell $(PYTHON) - <<'PY'\nimport sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '')\nPY)
PYBIND11_INCLUDES := $(shell $(PYTHON) -m pybind11 --includes)

CXXFLAGS ?= -O3 -std=c++17 -fPIC
LDFLAGS ?=

CORE_SOURCES := cpp/core/entanglement_core.cpp cpp/bindings/pybind_module.cpp
CORE_TARGET := python/rnaknotdetector_core$(EXT_SUFFIX)

.PHONY: all clean

all: $(CORE_TARGET)

$(CORE_TARGET): $(CORE_SOURCES)
	$(CXX) $(CXXFLAGS) $(PYBIND11_INCLUDES) -Icpp/core $^ -shared -o $@ $(LDFLAGS)

clean:
	rm -f $(CORE_TARGET)
