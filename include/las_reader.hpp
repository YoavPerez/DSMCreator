#include "las_header.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstring>
#include <tuple>


/* 
Point Data Format 3 – Structure (34 bytes)
Offset (bytes)	Size (bytes)	Field name	Type	Description
0	4	X	int32_t	X coordinate (to be scaled and offset)
4	4	Y	int32_t	Y coordinate
8	4	Z	int32_t	Z coordinate
12	2	Intensity	uint16_t	LIDAR return intensity
14	1	Return Info	uint8_t	Bit-packed: Return Number, Num Returns, Scan Dir Flag, Edge of Flight Line
15	1	Classification	uint8_t	ASPRS classification code
16	1	Scan Angle Rank	int8_t	Angle from -90 to +90
17	1	User Data	uint8_t	User-defined (optional)
18	2	Point Source ID	uint16_t	Typically flightline ID
20	8	GPS Time	float64	Time of return
28	2	Red	uint16_t	Red channel of color (0–65535)
30	2	Green	uint16_t	Green channel
32	2	Blue	uint16_t	Blue channel
*/

struct LASPointData{
    uint32_t x; // X coordinate (scaled and offset)
    uint32_t y; // Y coordinate
    uint32_t z; // Z coordinate
    uint16_t intensity; // LIDAR return intensity
    uint8_t returnInfo; // Bit-packed: Return Number, Num Returns, Scan Dir Flag, Edge of Flight Line
    uint8_t classification; // ASPRS classification code
    int8_t scanAngleRank; // Angle from -90 to +90
    uint8_t userData; // User-defined (optional)
    uint16_t pointSourceID; // Typically flightline ID
    double gpsTime; // Time of return
    uint16_t red; // Red channel of color (0–65535)
    uint16_t green; // Green channel
    uint16_t blue; // Blue channel
};

struct LASReader
{
    public:
        std::ifstream file;
        LASHeader header;
        uint32_t minX;
        uint32_t minY;
        uint32_t maxX;
        uint32_t maxY;
        double resolution;


    LASReader(const std::string& filename, double res = 0.3) {
        file.open(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Could not open LAS file: " + filename);
        }
        if (!header.read(file)) {
            throw std::runtime_error("Failed to read LAS header.");
        }
        file.seekg(header.dataOffset, std::ios::beg);
        resolution = res;
        minX = std::numeric_limits<uint32_t>::max();
        minY = std::numeric_limits<uint32_t>::max();
        maxX = std::numeric_limits<uint32_t>::lowest();
        maxY = std::numeric_limits<uint32_t>::lowest();
    }

    ~LASReader() {
        if (file.is_open()) {
            file.close();
        }
    }
    
    bool is_open() const {
        return file.is_open();
    }

    const LASHeader& get_header() const {
        return header;
    }

    std::vector<LASPointData> read_points();
    std::vector<uint32_t> create_DSM(const std::vector<LASPointData>& points, double minX, double minY, double maxX, double maxY, double res);
    void refine_DSM(std::vector<uint32_t>& dsm, int width, int height);
};
