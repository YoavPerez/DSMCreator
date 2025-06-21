#include <tiffio.h>
#include <tiff.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <chrono>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <omp.h>

#define MIN_FLOAT -3.402823466e+38F
#define BATCH_SIZE 1000000

struct LASHeader {
    uint32_t dataOffset;
    uint32_t pointCount;
    uint16_t pointRecordLength;
    double scaleX, scaleY, scaleZ;
    double offsetX, offsetY, offsetZ;
};

bool readHeader(std::ifstream& file, LASHeader& header) {
    file.seekg(24, std::ios::beg);
    uint8_t version_major = 0, version_minor = 0;
    file.read(reinterpret_cast<char*>(&version_major), 1);
    file.read(reinterpret_cast<char*>(&version_minor), 1);
    std::cout << "LAS version: " << (int)version_major << "." << (int)version_minor << std::endl;

    if (version_major != 1 || version_minor > 2) {
        std::cerr << "Only LAS 1.0 to 1.2 supported." << std::endl;
        return false;
    }

    file.seekg(96, std::ios::beg);
    file.read(reinterpret_cast<char*>(&header.dataOffset), 4);

    file.seekg(104, std::ios::beg);
    uint8_t format_id;
    file.read(reinterpret_cast<char*>(&format_id), 1);
    std::cout << "Point Data Format ID: " << (int)format_id << std::endl;
    if (format_id != 3) {
        std::cerr << "Only format 3 supported in this example." << std::endl;
        return false;
    }

    file.seekg(105);
    file.read(reinterpret_cast<char*>(&header.pointRecordLength), 2);

    file.seekg(107);
    file.read(reinterpret_cast<char*>(&header.pointCount), 4);

    file.seekg(131);
    file.read(reinterpret_cast<char*>(&header.scaleX), 8);
    file.read(reinterpret_cast<char*>(&header.scaleY), 8);
    file.read(reinterpret_cast<char*>(&header.scaleZ), 8);
    file.read(reinterpret_cast<char*>(&header.offsetX), 8);
    file.read(reinterpret_cast<char*>(&header.offsetY), 8);
    file.read(reinterpret_cast<char*>(&header.offsetZ), 8);

    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./las_dsm input.las [resolution]" << std::endl;
        return 1;
    }

    std::ifstream las(argv[1], std::ios::binary);
    if (!las) {
        std::cerr << "Could not open LAS file." << std::endl;
        return 1;
    }

    double resolution = (argc >= 3) ? std::stod(argv[2]) : 0.3;
    auto startRead = std::chrono::high_resolution_clock::now();

    LASHeader header;
    if (!readHeader(las, header)) return 1;

    std::cout << "Point count: " << header.pointCount << ", record size: " << header.pointRecordLength << std::endl;
    las.seekg(header.dataOffset, std::ios::beg);

    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();

    std::vector<uint8_t> buffer(BATCH_SIZE * header.pointRecordLength);

    std::vector<std::tuple<float, float, float>> xyz;
    xyz.reserve(header.pointCount);


    for (uint32_t done = 0; done < header.pointCount;) {
        size_t count = std::min<size_t>(BATCH_SIZE, header.pointCount - done);
        las.read(reinterpret_cast<char*>(buffer.data()), count * header.pointRecordLength);

        #pragma omp parallel for
        for (size_t i = 0; i < count; ++i) {
            size_t base = i * header.pointRecordLength;
            int32_t ix, iy, iz;
            std::memcpy(&ix, &buffer[base + 0], 4);
            std::memcpy(&iy, &buffer[base + 4], 4);
            std::memcpy(&iz, &buffer[base + 8], 4);

            float x = static_cast<float>(ix * header.scaleX + header.offsetX);
            float y = static_cast<float>(iy * header.scaleY + header.offsetY);
            float z = static_cast<float>(iz * header.scaleZ + header.offsetZ);

            #pragma omp critical
            {
                xyz.emplace_back(x, y, z);
                minX = std::min(minX, (double)x);
                maxX = std::max(maxX, (double)x);
                minY = std::min(minY, (double)y);
                maxY = std::max(maxY, (double)y);
            }
        }

        done += count;
        std::cout << "\rRead " << done << "/" << header.pointCount << std::flush;
    }

    auto endRead = std::chrono::high_resolution_clock::now();
    std::cout << "\nReading complete in " << std::chrono::duration<double>(endRead - startRead).count() << " seconds.\n";

    int width = static_cast<int>(std::ceil((maxX - minX) / resolution));
    int height = static_cast<int>(std::ceil((maxY - minY) / resolution));

    std::vector<float> dsm(width * height, MIN_FLOAT);

    auto startDSM = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (size_t i = 0; i < xyz.size(); ++i) {
        float x = std::get<0>(xyz[i]);
        float y = std::get<1>(xyz[i]);
        float z = std::get<2>(xyz[i]);

        int col = static_cast<int>((x - minX) / resolution);
        int row = static_cast<int>((y - minY) / resolution);

        if (col >= 0 && col < width && row >= 0 && row < height) {
            size_t idx = row * width + col;
            #pragma omp critical
            dsm[idx] = std::max(dsm[idx], z);
        }
    }

    // Parallelize both row and col loops using collapse(2)
    #pragma omp parallel for collapse(2)
    for (int row = 1; row < height - 1; ++row) {
        for (int col = 1; col < width - 1; ++col) {
            size_t idx = row * width + col;
            if (dsm[idx] == MIN_FLOAT) {
                float sum = 0;
                int count = 0;
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        size_t n_idx = (row + dy) * width + (col + dx);
                        float val = dsm[n_idx];
                        if (val != MIN_FLOAT) {
                            sum += val;
                            count++;
                        }
                    }
                }
                if (count > 0) {
                    dsm[idx] = sum / count;  // Mean interpolation
                }
            }
        }
    }


    auto endDSM = std::chrono::high_resolution_clock::now();
    std::cout << "DSM computation took " << std::chrono::duration<double>(endDSM - startDSM).count() << " seconds.\n";

    // output DSM statistics
    float minDSM = *std::min_element(dsm.begin(), dsm.end());
    float maxDSM = *std::max_element(dsm.begin(), dsm.end());
    std::cout << "DSM size: " << width << " x " << height << " pixels" << std::endl;
    std::cout << "DSM min: " << minDSM << ", max: " << maxDSM << std::endl;
    
    
    // Export DSM as a 32-bit float GeoTIFF using libtiff

    const char* tiff_filename = "dsm.tif";
    TIFF* tif = TIFFOpen(tiff_filename, "w");
    if (!tif) {
        std::cerr << "Could not open " << tiff_filename << " for writing." << std::endl;
        return 1;
    }

    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

    // Optionally set GeoTIFF tags here (requires libgeotiff)

    // Write row by row (TIFF expects rows from top to bottom)
    for (int row = 0; row < height; ++row) {
        float* buf = &dsm[row * width];
        if (TIFFWriteScanline(tif, buf, row, 0) < 0) {
            std::cerr << "Error writing TIFF scanline." << std::endl;
            TIFFClose(tif);
            return 1;
        }
    }
    TIFFClose(tif);
    std::cout << "DSM written to " << tiff_filename << " (size: " << width << " x " << height << ")" << std::endl;

    return 0;
}
