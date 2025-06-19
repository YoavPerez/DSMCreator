CXX = g++
CXXFLAGS = -std=c++17 -O2 -fopenmp
INCLUDES = 
LDFLAGS = -ltiff -fopenmp
LIBLAS_LIB_PATH = /usr/local/lib

SRC = las_parser.cpp
OUT = main

all: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS)

run: $(OUT)
	LD_LIBRARY_PATH=$(LIBLAS_LIB_PATH):$$LD_LIBRARY_PATH ./$(OUT)

clean:
	rm -f $(OUT)