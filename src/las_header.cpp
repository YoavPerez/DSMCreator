#include "las_header.hpp"
#include <iostream>

bool LASHeader::read(std::ifstream &file)
{
    file.seekg(24, std::ios::beg);
    uint8_t version_major = 0, version_minor = 0;
    file.read(reinterpret_cast<char *>(&version_major), 1);
    file.read(reinterpret_cast<char *>(&version_minor), 1);
    std::cout << "LAS version: " << (int)version_major << "." << (int)version_minor << std::endl;

    if (version_major != 1 || version_minor > 2)
    {
        std::cerr << "Only LAS 1.0 to 1.2 supported." << std::endl;
        return false;
    }

    file.seekg(96, std::ios::beg);
    file.read(reinterpret_cast<char *>(&dataOffset), 4);

    file.seekg(104, std::ios::beg);
    uint8_t format_id;
    file.read(reinterpret_cast<char *>(&format_id), 1);
    std::cout << "Point Data Format ID: " << (int)format_id << std::endl;
    if (format_id != 3)
    {
        std::cerr << "Only format 3 supported in this example." << std::endl;
        return false;
    }

    file.seekg(105);
    file.read(reinterpret_cast<char *>(&pointRecordLength), 2);

    file.seekg(107);
    file.read(reinterpret_cast<char *>(&pointCount), 4);

    file.seekg(131);
    file.read(reinterpret_cast<char *>(&scaleX), 8);
    file.read(reinterpret_cast<char *>(&scaleY), 8);
    file.read(reinterpret_cast<char *>(&scaleZ), 8);
    file.read(reinterpret_cast<char *>(&offsetX), 8);
    file.read(reinterpret_cast<char *>(&offsetY), 8);
    file.read(reinterpret_cast<char *>(&offsetZ), 8);

    return true;
}
