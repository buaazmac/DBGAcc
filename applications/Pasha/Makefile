ZLIB_DIR = zlib
COMMON_DIR = Common
PREGRAPH_DIR = PreGraph
#GRAPH_DIR = Graph

BIN_DIR = bin
OBJS_DIR = objs
LIB_DIR = lib

#directory for Intel Threading Building Blocks (TBB)
#can be downloaded from http://threadingbuildingblocks.org/download.php
export TBB_DIR = ~/install
export TBB_INC = -I $(TBB_DIR)/include/tbb
export TBB_LIB = -L $(TBB_DIR)/lib -ltbb -ltbbmalloc

#directory for Message Passing Interface (MPI) library
#please give the correct values for these macros
export MPI_CXX = mpicxx
export MPI_DIR = /usr/lib/openmpi
export MPI_INC = -I $(MPI_DIR)/include
export MPI_LIB = -L $(MPI_DIR)/lib -lmpi

#common options for general c++ compiler
export MY_CXX = g++
export OMP_FLAG = -fopenmp
export MYLIBS = -L $(LIB_DIR) -lz

#executable binaries
#export KMERGEN = pasha-kmergen
#export PREGRAPH = pasha-pregraph
#export GRAPH= pasha-graph
export STEP1 = pasha-step1

#all: dir zlib kmergen_pregraph graph
all: dir zlib step1

zlib:dir
	cd $(ZLIB_DIR); ./configure; make; rm minigzip.o; rm example.o
	cp -rf $(ZLIB_DIR)/libz.a $(LIB_DIR)

#kmergen_pregraph: dir
#	make -C $(PREGRAPH_DIR)
#	cp -rf $(PREGRAPH_DIR)/$(KMERGEN) $(BIN_DIR)
#	cp -rf $(PREGRAPH_DIR)/$(PREGRAPH) $(BIN_DIR)

step1: dir
	make -C $(PREGRAPH_DIR)
	cp -rf $(PREGRAPH_DIR)/$(STEP1) $(BIN_DIR)

#graph: dir
#	make -C $(GRAPH_DIR) STAGE=ALSHA_GRAPH
#	cp -rf $(GRAPH_DIR)/$(GRAPH) $(BIN_DIR)

dir:
	mkdir -p $(OBJS_DIR)
	mkdir -p $(BIN_DIR)
	mkdir -p $(LIB_DIR)

clean:
	make -C $(ZLIB_DIR) clean
	-rm -rf $(BIN_DIR)
	-rm -rf $(OBJS_DIR)
	-rm -rf $(LIB_DIR)
