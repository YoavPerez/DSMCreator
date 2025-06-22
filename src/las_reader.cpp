#include "las_reader.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <tuple>
#include <cmath>

#define BATCH_SIZE 1'000'000

std::vector<LASPointData> LASReader::read_points()
{
    std::vector<LASPointData> points(header.pointCount);

    std::vector<uint8_t> buffer(BATCH_SIZE * header.pointRecordLength);

    uint32_t total_points = header.pointCount;
    uint16_t record_len = header.pointRecordLength;


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

            // Read X, Y, Z as int32_t (little-endian)
            std::memcpy(&pt.x, &buffer[base + 0], sizeof(int32_t));
            std::memcpy(&pt.y, &buffer[base + 4], sizeof(int32_t));
            std::memcpy(&pt.z, &buffer[base + 8], sizeof(int32_t));

            // Intensity (uint16_t)
            std::memcpy(&pt.intensity, &buffer[base + 12], sizeof(uint16_t));

            pt.returnInfo      = buffer[base + 14];
            pt.classification  = buffer[base + 15];
            pt.scanAngleRank   = buffer[base + 16];
            pt.userData        = buffer[base + 17];

            // Point Source ID (uint16_t)
            std::memcpy(&pt.pointSourceID, &buffer[base + 18], sizeof(uint16_t));

            // gpsTime (double, 8 bytes, little-endian)
            std::memcpy(&pt.gpsTime, &buffer[base + 20], sizeof(double));

            // Color (uint16_t, little-endian)
            std::memcpy(&pt.red,   &buffer[base + 28], sizeof(uint16_t));
            std::memcpy(&pt.green, &buffer[base + 30], sizeof(uint16_t));
            std::memcpy(&pt.blue,  &buffer[base + 32], sizeof(uint16_t));
            points[done + i] = pt;

            // Update min/max coordinates (atomic for thread safety)
            // #pragma omp critical
            // {
                int x = static_cast<int>(pt.x);
                int y = static_cast<int>(pt.y);
                if (x < minX) minX = x;
                if (x > maxX) maxX = x;
                if (y < minY) minY = y;
                if (y > maxY) maxY = y;
            }
        }
    }

    std::cout << "Reading complete. Total points read: " << points.size() << std::endl;
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
    std::cout << "\rProcessed " << points.size() << "/" << points.size() << " points" << std::endl;
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