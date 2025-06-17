CXX = g++
CXXFLAGS = -std=c++17 -O2
LDFLAGS = $(shell pkg-config --libs pdal)
INCLUDES = $(shell pkg-config --cflags pdal)

SRC = las_info.cpp
OUT = main

all: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(OUT)