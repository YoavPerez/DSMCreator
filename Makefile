CXX = g++
CXXFLAGS = -std=c++17 -O2 -fopenmp
INCLUDES = -Iinclude $(shell gdal-config --cflags)
LDFLAGS = -ltiff -fopenmp $(shell gdal-config --libs)

SRC_DIR = src
OBJ_DIR = build
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
OUT = las_dsm

all: $(OUT)

$(OUT): $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -rf $(OBJ_DIR) $(OUT)
