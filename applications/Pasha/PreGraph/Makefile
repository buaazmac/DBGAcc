CXX_SOURCES = Alsha.cpp AlshaMessage.cpp AlshaComms.cpp  AlshaThread.cpp AlshaHash.cpp
#KMERGEN_SOURCES = KmerGen.cpp
#PREGRAPH_SOURCES = PreGraph.cpp
STEP1_SOURCES = pasha-step1.cpp

COMMON_DIR = ../Common
COMMON_INC = -I $(COMMON_DIR)
COMMON_LIB = -L ../lib -lcommon

# export CXX = $(MPI_CXX)
export CXXFLAGS = -O3 -m64 -funroll-loops -Wno-deprecated $(MPI_INC) $(TBB_INC) $(COMMON_INC) -DHAVE_GOOGLE_SPARSE_HASH_SET -I ..

OBJS_DIR = ../objs
SOURCE_DIR = .
#KMERGEN_OBJS += $(patsubst %.cpp, $(OBJS_DIR)/%.cpp.kmer.o, $(KMERGEN_SOURCES))
#KMERGEN_OBJS += $(patsubst %.cpp, $(OBJS_DIR)/%.cpp.kmer.o, $(CXX_SOURCES))
#PREGRAPH_OBJS += $(patsubst %.cpp, $(OBJS_DIR)/%.cpp.pregraph.o, $(PREGRAPH_SOURCES))
#PREGRAPH_OBJS += $(patsubst %.cpp, $(OBJS_DIR)/%.cpp.pregraph.o, $(CXX_SOURCES))
STEP1_OBJS += $(patsubst %.cpp, $(OBJS_DIR)/%.cpp.step1.o, $(STEP1_SOURCES))
STEP1_OBJS += $(patsubst %.cpp, $(OBJS_DIR)/%.cpp.step1.o, $(CXX_SOURCES))

#all:kmergen pregraph
all: step1

#kmergen: kmergen_clean $(KMERGEN_OBJS)
#	$(MPI_CXX) $(CXXFLAGS) -o $(KMERGEN) $(SHARED_OBJS) $(KMERGEN_OBJS) $(COMMON_LIB) $(TBB_LIB) $(LIBS) $(MYLIBS)
#	strip $(KMERGEN)

step1: step1_clean $(STEP1_OBJS)
	$(info ***** step1: step1_clean finished****)
	$(MPI_CXX) $(CXXFLAGS) -o $(STEP1) $(SHARED_OBJS) $(STEP1_OBJS) $(COMMON_LIB) $(TBB_LIB) $(LIBS) $(MYLIBS) -lpthread
	strip $(STEP1)


#kmergen_clean:
#	rm -rf $(OBJS_DIR)/*
#	make -C $(COMMON_DIR) STAGE=ALSHA_KMERGEN USE_MPI=1
#	make -C $(COMMON_DIR) STAGE=ALSHA_PREGRAPH USE_MPI=1


step1_clean:
	rm -rf $(OBJS_DIR)/*
	make -C $(COMMON_DIR) STAGE=ALSHA_KMERGEN USE_MPI=1
	$(info ***** step1_clean make1****)
	make -C $(COMMON_DIR) STAGE=ALSHA_PREGRAPH USE_MPI=1
	$(info ***** step1_clean make2****)

#pregraph: pregraph_clean $(PREGRAPH_OBJS)
#	$(MPI_CXX) $(CXXFLAGS) -o $(PREGRAPH) $(SHARED_OBJS) $(PREGRAPH_OBJS) $(COMMON_LIB) $(TBB_LIB) $(LIBS) $(MYLIBS)
#	strip $(PREGRAPH)

#$(OBJS_DIR)/%.cpp.kmer.o: $(SOURCE_DIR)/%.cpp
#	$(MPI_CXX) $(CXXFLAGS) -o $@ -c $<

#$(OBJS_DIR)/%.cpp.pregraph.o: $(SOURCE_DIR)/%.cpp
#	$(MPI_CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJS_DIR)/%.cpp.step1.o: $(SOURCE_DIR)/%.cpp
	$(MPI_CXX) $(CXXFLAGS) -o $@ -c $<


#pregraph_clean:
#	rm -rf $(OBJS_DIR)/*
#	make -C $(COMMON_DIR) STAGE=ALSHA_PREGRAPH USE_MPI=1

