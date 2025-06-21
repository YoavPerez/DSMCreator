#include "las_reader.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./las_dsm input.las [resolution]" << std::endl;
        return 1;
    }

    std::string las = argv[1];

    double resolution = (argc >= 3) ? std::stod(argv[2]) : 0.3;


    auto startRead = std::chrono::high_resolution_clock::now();

    LASReader reader(las, resolution);

    std::cout << "Point count: " << reader.header.pointCount << ", record size: " << reader.header.pointRecordLength << std::endl;
    std::cout << "Scale: (" << reader.header.scaleX << ", " << reader.header.scaleY << ", " << reader.header.scaleZ << ")" << std::endl;
    std::cout << "Offset: (" << reader.header.offsetX << ", " << reader.header.offsetY << ", " << reader.header.offsetZ << ")" << std::endl;
    std::cout << "Data offset: " << reader.header.dataOffset << std::endl; 
    std::cout << "Point data format ID: " << static_cast<int>(reader.header.pointRecordLength) << std::endl;
    
    int minX, minY, maxX, maxY;
    auto points = reader.read_points(&minX, &minY, &maxX, &maxY);
    double minXCoord = reader.header.offsetX + reader.header.scaleX * minX;
    double maxXCoord = reader.header.offsetX + reader.header.scaleX * maxX;
    double minYCoord = reader.header.offsetY + reader.header.scaleY * minY;
    double maxYCoord = reader.header.offsetY + reader.header.scaleY * maxY;
    auto dsm = reader.create_DSM(points, minXCoord, minYCoord, maxXCoord, maxYCoord, resolution);
    reader.refine_DSM(dsm, std::ceil((maxX - minX) / resolution), std::ceil((maxY - minY) / resolution));
    auto endRead = std::chrono::high_resolution_clock::now();
    std::cout << "\nReading complete in " << std::chrono::duration<double>(endRead - startRead).count() << " seconds.\n";


    
    return 0;
};
