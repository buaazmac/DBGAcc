CXX_SOURCES = AlshaUtils.cpp AlshaKmer.cpp AlshaParams.cpp AlshaFileParser.cpp \
		AlshaFastaFileParser.cpp AlshaFastqFileParser.cpp AlshaHistogram.cpp

STAGE ?=
USE_MPI ?= 0

ifeq ($(USE_MPI), 1)
	CXXFLAGS += -DALSHA_MPI
endif
CXXFLAGS += -D$(STAGE)

OBJS_DIR = ../objs
SOURCE_DIR = .
OBJS += $(patsubst %.cpp, $(OBJS_DIR)/%.cpp.o, $(CXX_SOURCES))

COMMON_LIB = libcommon.a

all: $(OBJS)
	ar -r $(COMMON_LIB) $(OBJS)
	mv $(COMMON_LIB) ../lib

$(OBJS_DIR)/%.cpp.o: $(SOURCE_DIR)/%.cpp
	$(MPI_CXX) $(CXXFLAGS) -o $@ -c $<

