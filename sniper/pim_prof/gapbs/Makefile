# See LICENSE.txt for license details.
TARGET = gapbs
include ../shared/Makefile.shared

#CXX = /usr/bin/clang++-7
#OPT = /usr/bin/opt-7
#PIMPROF_ROOT = /home/warsier/Documents/PIMProf/

CXX_FLAGS = -std=c++11 -O3 -Wall -fopenmp -g
#CXX_FLAGS = -std=c++11 -O3 -Wall -fopenmp -g -march=skylake-avx512 -I$(PIMPROF_ROOT)/LLVMAnalysis
#CXX_FLAGS_2 = -std=c++11 -O3 -Wall -g -I$(PIMPROF_ROOT)/LLVMAnalysis
#CXX_INJ_FLAGS = $(CXX_FLAGS) -Xclang -load -Xclang $(PIMPROF_ROOT)/build/LLVMAnalysis/libAnnotationInjection.so
#CXX_OFF_FLAGS += $(CXX_FLAGS) -Xclang -load -Xclang $(PIMPROF_ROOT)/build/LLVMAnalysis/libOffloaderInjection.so
#LD_FLAGS = -L$(PIMPROF_ROOT)/build/LLVMAnalysis -Wl,-rpath=$(PIMPROF_ROOT)/build/LLVMAnalysis -lPIMProfAnnotation -lzsimhooks
#LD_FLAGS = -L$(PIMPROF_ROOT)/build/LLVMAnalysis -Wl,-rpath=$(PIMPROF_ROOT)/build/LLVMAnalysis -lPIMProfAnnotation -lzsimhooks

KERNELS = bc bfs bfs_cpu bfs_pim cc cc_sv pr pr_pim sssp tc
SUITE = $(KERNELS) converter
#SUITE_ZSIM = $(SUITE:=.zsim)
#SUITE_FINAL = $(SUITE:=.final)

.PHONY: all
all: $(SUITE)

% : src/%.cc src/*.h
	$(CXX) $(CXX_FLAGS) $< -o $@ $(LD_FLAGS) -lm $(SNIPER_LDFLAGS) $(SNIPER_CFLAGS)

.PHONY: clean
clean:
	rm -f $(SUITE) test/out/*

clean_kernel:
	rm -f $(KERNELS)
