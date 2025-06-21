#include "las_reader.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <tuple>
#include <cmath>

#define BATCH_SIZE 1'000'000

std::vector<LASPointData> LASReader::read_points(int *minX, int *minY, int *maxX, int *maxY)
{
    std::vector<LASPointData> points(header.pointCount);

    std::vector<uint8_t> buffer(BATCH_SIZE * header.pointRecordLength);

    uint32_t total_points = header.pointCount;
    uint16_t record_len = header.pointRecordLength;

    // Initialize min/max with extreme values
    *minX = std::numeric_limits<int>::max();
    *minY = std::numeric_limits<int>::max();
    *maxX = std::numeric_limits<int>::min();
    *maxY = std::numeric_limits<int>::min();

    for (uint32_t done = 0; done < total_points; done += BATCH_SIZE)
    {
        size_t count = std::min<size_t>(BATCH_SIZE, total_points - done);
        std::cout << "\rRead " << std::min<uint32_t>(static_cast<uint32_t>(done + count), total_points) << "/" << total_points << std::flush;
        file.read(reinterpret_cast<char *>(buffer.data()), count * record_len);

        // #pragma omp parallel for
        for (size_t i = 0; i < count; ++i)
        {
            size_t base = i * record_len;
            LASPointData pt;

            pt.x =  static_cast<int32_t>(
                static_cast<uint32_t>(buffer[base + 0]) |
                (static_cast<uint32_t>(buffer[base + 1]) << 8) |
                (static_cast<uint32_t>(buffer[base + 2]) << 16) |
                (static_cast<uint32_t>(buffer[base + 3]) << 24)
            );
            pt.y =  static_cast<int32_t>(
                static_cast<uint32_t>(buffer[base + 4]) |
                (static_cast<uint32_t>(buffer[base + 5]) << 8) |
                (static_cast<uint32_t>(buffer[base + 6]) << 16) |
                (static_cast<uint32_t>(buffer[base + 7]) << 24)
            );
            pt.z =  static_cast<int32_t>(
                static_cast<uint32_t>(buffer[base + 8]) |
                (static_cast<uint32_t>(buffer[base + 9]) << 8) |
                (static_cast<uint32_t>(buffer[base + 10]) << 16) |
                (static_cast<uint32_t>(buffer[base + 11]) << 24)
            );
            pt.intensity = static_cast<uint16_t>(
                static_cast<uint16_t>(buffer[base + 12]) |
                (static_cast<uint16_t>(buffer[base + 13]) << 8)
            );
            pt.returnInfo = buffer[base + 14];
            pt.classification = buffer[base + 15];
            pt.scanAngleRank = buffer[base + 16];
            pt.userData = buffer[base + 17];
            pt.pointSourceID = static_cast<uint16_t>(
                static_cast<uint16_t>(buffer[base + 18]) |
                (static_cast<uint16_t>(buffer[base + 19]) << 8)
            );
            // gpsTime is 8 bytes, little-endian double
            uint64_t gpsTimeBits =
                static_cast<uint64_t>(buffer[base + 20]) |
                (static_cast<uint64_t>(buffer[base + 21]) << 8) |
                (static_cast<uint64_t>(buffer[base + 22]) << 16) |
                (static_cast<uint64_t>(buffer[base + 23]) << 24) |
                (static_cast<uint64_t>(buffer[base + 24]) << 32) |
                (static_cast<uint64_t>(buffer[base + 25]) << 40) |
                (static_cast<uint64_t>(buffer[base + 26]) << 48) |
                (static_cast<uint64_t>(buffer[base + 27]) << 56);
            std::memcpy(&pt.gpsTime, &gpsTimeBits, sizeof(double)); // This is allowed, as it's not memcpy from buffer

            pt.red = static_cast<uint16_t>(
                static_cast<uint16_t>(buffer[base + 28]) |
                (static_cast<uint16_t>(buffer[base + 29]) << 8)
            );
            pt.green = static_cast<uint16_t>(
                static_cast<uint16_t>(buffer[base + 30]) |
                (static_cast<uint16_t>(buffer[base + 31]) << 8)
            );
            pt.blue = static_cast<uint16_t>(
                static_cast<uint16_t>(buffer[base + 32]) |
                (static_cast<uint16_t>(buffer[base + 33]) << 8)
            );
            points[done + i] = pt;

            // Update min/max coordinates (atomic for thread safety)
            // #pragma omp critical
            // {
                int x = static_cast<int>(pt.x);
                int y = static_cast<int>(pt.y);
                if (x < *minX) *minX = x;
                if (x > *maxX) *maxX = x;
                if (y < *minY) *minY = y;
                if (y > *maxY) *maxY = y;
            // }
        }
    }

    std::cout << "\nReading complete. Total points read: " << points.size() << std::endl;
    double minXCoord = header.offsetX + header.scaleX * (*minX);
    double maxXCoord = header.offsetX + header.scaleX * (*maxX);
    double minYCoord = header.offsetY + header.scaleY * (*minY);
    double maxYCoord = header.offsetY + header.scaleY * (*maxY);

    std::cout << "Coordinate range: X [" << minXCoord << ", " << maxXCoord << "], Y [" << minYCoord << ", " << maxYCoord << "]" << std::endl;
    std::cout << "Resolution: " << resolution << std::endl;
    return points;
}

std::vector<uint32_t> LASReader::create_DSM(const std::vector<LASPointData>& points, double minX, double minY, double maxX, double maxY, double res)
{
    // create DSM based on the points read
    int width = static_cast<int>(std::ceil((maxX - minX) / res));
    int height = static_cast<int>(std::ceil((maxY - minY) / res));
    std::cout << "DSM dimensions: " << width << " x " << height << " pixels" << std::endl;
    std::cout << "Min X: " << minX << ", Min Y: " << minY << ", Max X: " << maxX << ", Max Y: " << maxY << std::endl;
    if (width <= 0 || height <= 0) {
        throw std::runtime_error("Invalid DSM dimensions: width and height must be positive.");
    }
    if (static_cast<int64_t>(width) * static_cast<int64_t>(height) > 1000000000) {
        throw std::runtime_error("DSM dimensions too large, possible memory allocation failure.");
    }
    std::vector<uint32_t> dsm(static_cast<size_t>(width) * static_cast<size_t>(height), std::numeric_limits<uint32_t>::min());
    std::cout << "Creating DSM with resolution: " << res << " meters" << std::endl;
    std::cout << "DSM size: " << width << " x " << height << " pixels" << std::endl;
    int i = 0;
    for (const auto& pt : points)
    {
        if (i % 100000 == 0) // Show progress every 100k points
            // Show progress bar
            std::cout << "\rProcessed " << (&pt - &points[0] + 1) << "/" << points.size() << " points" << std::flush;
        int col = static_cast<int>((pt.x - minX) / res);
        int row = static_cast<int>((pt.y - minY) / res);

        if (col >= 0 && col < width && row >= 0 && row < height)
        {
            size_t idx = row * width + col;
            dsm[idx] = std::max(dsm[idx], pt.z);
        }
        i++;
        // show progress bar
    }
    std::cout << "\rProcessed " << points.size() << "/" << points.size() << " points" << std::flush;
    return dsm;

}

void LASReader::refine_DSM(std::vector<uint32_t>& dsm, int width, int height)
{
    // Refine DSM by interpolating missing values
    auto startDSM = std::chrono::high_resolution_clock::now();

    // #pragma omp parallel for collapse(2)
    for (int row = 1; row < height - 1; ++row) {
        for (int col = 1; col < width - 1; ++col) {
            size_t idx = row * width + col;
            if (dsm[idx] == std::numeric_limits<uint32_t>::min()) {
                float sum = 0;
                int count = 0;
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        // size_t n_idx = (row + dy) * width + (col + dx);
                        // float val = dsm[n_idx];
                        // if (val != std::numeric_limits<uint32_t>::min()) {
                        //     sum += val;
                        //     count++;
                        // }
                    }
                }
                if (count > 0) {
                    dsm[idx] = sum / count;  // Mean interpolation
                }
            }
        }
    }

    auto endDSM = std::chrono::high_resolution_clock::now();
    std::cout << "DSM refinement took " << std::chrono::duration<double>(endDSM - startDSM).count() << " seconds.\n";
}