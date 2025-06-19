struct LASHeader {
    uint32_t dataOffset;
    uint32_t pointCount;
    uint16_t pointRecordLength;
    double scaleX, scaleY, scaleZ;
    double offsetX, offsetY, offsetZ;
    bool read(std::ifstream& file);
};