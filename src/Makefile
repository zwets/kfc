CXXFLAGS += -std=c++14 -O3 -DNDEBUG -Wall -Wextra -pedantic -mtune=native

OBJS = kfc.o kmercounter.o kmerencoder.o seqreader.o utils.o 

LIBS =

HDRS = implpicker.h kmercounter.h tallyman.h kmerencoder.h kmercodec.h basecodec.h bitfiddle.h seqreader.h utils.h

TARGET = kfc

ifeq (,$(wildcard /usr/include/boost/iostreams/filter/gzip.hpp))
  $(warning "NOTE: kfc will be built without gzip support, so won't be able to open compressed fasta.")
  $(warning "      To build with decompression support, install libboost-iostream-dev and libz-dev,")
  $(warning "      then run 'make clean; make' again.")
  CXXFLAGS += -DNO_ZLIB
else
  LIBS += -Wl,-Bstatic -lboost_iostreams -lz -Wl,-Bdynamic
endif

ifdef NO_THREADS
  CXXFLAGS += -DNO_THREADS
else
  CXXFLAGS += -pthread
  LIBS += -lpthread
endif

$(TARGET): $(OBJS) $(HDRS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all: $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
	$(MAKE) -C unit-test clean

test:
	$(MAKE) -C unit-test test

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $<

.cpp.s:
	$(CXX) $(CXXFLAGS) -g -S $<

