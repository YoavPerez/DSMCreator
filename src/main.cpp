#include "las_reader.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <gdal_priv.h>
#include <cpl_conv.h>           // for CPLMalloc/CPLFree
#include <ogr_spatialref.h>     // for OGRSpatialReference



int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: ./las_dsm <input.las> <output.tif> [resolution]" << std::endl;
        return 1;
    }

    std::string las = argv[1];

    std::string output = argv[2];
    double resolution = (argc >= 4) ? std::stod(argv[3]) : 0.3;


    auto startRead = std::chrono::high_resolution_clock::now();

    LASReader reader(las, resolution);

    std::cout << "Point count: " << reader.header.pointCount << ", record size: " << reader.header.pointRecordLength << std::endl;
    std::cout << "Scale: (" << reader.header.scaleX << ", " << reader.header.scaleY << ", " << reader.header.scaleZ << ")" << std::endl;
    std::cout << "Offset: (" << reader.header.offsetX << ", " << reader.header.offsetY << ", " << reader.header.offsetZ << ")" << std::endl;
    std::cout << "Data offset: " << reader.header.dataOffset << std::endl; 
    std::cout << "Point data format ID: " << static_cast<int>(reader.header.pointRecordLength) << std::endl;

    auto points = reader.read_points();
    std::cout << "points read: " << points.size() << std::endl;
    auto dsm = reader.create_DSM(points);

    auto endRead = std::chrono::high_resolution_clock::now();
    std::cout << "\nReading complete in " << std::chrono::duration<double>(endRead - startRead).count() << " seconds.\n";
    
    int width = reader.width;
    int height = reader.height;
    double pixleSize = reader.resolution;
    double origin_x = reader.minXCoord;
    double origin_y = reader.maxYCoord; // Y is inverted in GDAL, so we use maxY for top-left corner

     GDALAllRegister();

    const char *pszFormat = "GTiff";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if (poDriver == nullptr) {
        throw std::runtime_error("GTiff driver not available.");
    }

    GDALDataset *poDstDS = poDriver->Create(output.c_str(), width, height, 1, GDT_Float64, nullptr);
    if (!poDstDS) {
        throw std::runtime_error("Failed to create output dataset.");
    }


    double adfGeoTransform[6] = {
        origin_x, resolution, 0,
        origin_y, 0, -resolution
    };
    poDstDS->SetGeoTransform(adfGeoTransform);

    // Set projection to UTM Zone 36N (WGS84)
    OGRSpatialReference oSRS;
    oSRS.SetUTM(36, TRUE); // TRUE for Northern Hemisphere
    oSRS.SetWellKnownGeogCS("WGS84");

    char *pszSRS_WKT = nullptr;
    oSRS.exportToWkt(&pszSRS_WKT);
    poDstDS->SetProjection(pszSRS_WKT);
    CPLFree(pszSRS_WKT);

    // Write the data and set NoData value
    GDALRasterBand *poBand = poDstDS->GetRasterBand(1);

    // Set NoData value to match your "no data" marker
    double noDataValue = std::numeric_limits<double>::min();
    poBand->SetNoDataValue(noDataValue);

    CPLErr err = poBand->RasterIO(GF_Write, 0, 0, width, height,
                                  (void *)dsm.data(), width, height, GDT_Float64,
                                  0, 0);
    if (err != CE_None) {
        std::cerr << "RasterIO write failed!" << std::endl;
    }

    GDALClose(poDstDS);

    return 0;
}
