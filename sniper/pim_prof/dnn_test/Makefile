# See LICENSE.txt for license details.
TARGET = gapbs
#include ../shared/Makefile.shared

DNNLROOT = /usr/local
SNIPER_ROOT = /media/zmac/SpinDisk/Research/PIM/full_sys_simulator/test/sniper_pim_prof/sniper

CXX = g++
CXX_FLAGS = -std=c++11 -O3 -Wall -fopenmp -g
#DNNL_LDFLAGS = -I${DNNLROOT}/include -L${DNNLROOT}/lib -ldnnl
SNIPER_LDFLAGS:= -L${SNIPER_ROOT}/lib -pthread -L${DNNLROOT}/lib -ldnnl
SNIPER_CFLAGS:=-mno-sse4 -mno-sse4.1 -mno-sse4.2 -mno-sse4a -mno-avx -mno-avx2 -I${SNIPER_ROOT}/include -I${DNNLROOT}/include 
#CXX_FLAGS = -std=c++11 -O3 -Wall -fopenmp -g -march=skylake-avx512 -I$(PIMPROF_ROOT)/LLVMAnalysis
#CXX_FLAGS_2 = -std=c++11 -O3 -Wall -g -I$(PIMPROF_ROOT)/LLVMAnalysis
#CXX_INJ_FLAGS = $(CXX_FLAGS) -Xclang -load -Xclang $(PIMPROF_ROOT)/build/LLVMAnalysis/libAnnotationInjection.so
#CXX_OFF_FLAGS += $(CXX_FLAGS) -Xclang -load -Xclang $(PIMPROF_ROOT)/build/LLVMAnalysis/libOffloaderInjection.so
#LD_FLAGS = -L$(PIMPROF_ROOT)/build/LLVMAnalysis -Wl,-rpath=$(PIMPROF_ROOT)/build/LLVMAnalysis -lPIMProfAnnotation -lzsimhooks
#LD_FLAGS = -L$(PIMPROF_ROOT)/build/LLVMAnalysis -Wl,-rpath=$(PIMPROF_ROOT)/build/LLVMAnalysis -lPIMProfAnnotation -lzsimhooks

KERNELS = test alexnet
SUITE = $(KERNELS)
#SUITE_ZSIM = $(SUITE:=.zsim)
#SUITE_FINAL = $(SUITE:=.final)

.PHONY: all
all: $(SUITE)

% : src/%.cpp src/*.h
	$(CXX) $(CXX_FLAGS) $< -o $@ $(LD_FLAGS) $(DNNL_LDFLAGS) -lm $(SNIPER_LDFLAGS) $(SNIPER_CFLAGS)

test : src/test.cpp
	$(CXX) $(CXX_FLAGS) src/example_utils.hpp src/test.cpp -o test $(LD_FLAGS) $(DNNL_LDFLAGS) -lm $(SNIPER_LDFLAGS) $(SNIPER_CFLAGS)

alexnet : src/alexnet.cpp
	$(CXX) $(CXX_FLAGS) src/example_utils.hpp src/alexnet.cpp -o alexnet $(LD_FLAGS) $(DNNL_LDFLAGS) -lm $(SNIPER_LDFLAGS) $(SNIPER_CFLAGS)

.PHONY: clean
clean:
	rm -f $(SUITE) test/out/* sim.* mem.trace

clean_kernel:
	rm -f $(KERNELS)
