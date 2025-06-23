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

            pt.returnInfo = buffer[base + 14];
            pt.classification = buffer[base + 15];
            pt.scanAngleRank = buffer[base + 16];
            pt.userData = buffer[base + 17];

            // Point Source ID (uint16_t)
            std::memcpy(&pt.pointSourceID, &buffer[base + 18], sizeof(uint16_t));

            points[done + i] = pt;

            int x = static_cast<int>(pt.x);
            int y = static_cast<int>(pt.y);
            int z = static_cast<int>(pt.z);
            if (x < minX)
                minX = x;
            if (x > maxX)
                maxX = x;
            if (y < minY)
                minY = y;
            if (y > maxY)
                maxY = y;
            if (z < minZ)
                minZ = z;
            if (z > maxZ)
                maxZ = z;
        }
    }
    minXCoord = static_cast<double>(minX) * header.scaleX + header.offsetX;
    minYCoord = static_cast<double>(minY) * header.scaleY + header.offsetY;
    maxXCoord = static_cast<double>(maxX) * header.scaleX + header.offsetX;
    maxYCoord = static_cast<double>(maxY) * header.scaleY + header.offsetY;
    minZCoord = static_cast<double>(minZ) * header.scaleZ + header.offsetZ;
    maxZCoord = static_cast<double>(maxZ) * header.scaleZ + header.offsetZ;

    width = static_cast<int>(std::ceil((maxXCoord - minXCoord) / resolution));
    height = static_cast<int>(std::ceil((maxYCoord - minYCoord) / resolution));
    if (width <= 0 || height <= 0)
    {
        throw std::runtime_error("Invalid dimensions: width and height must be positive.");
    }
    if (static_cast<int64_t>(width) * static_cast<int64_t>(height) > 1000000000)
    {
        throw std::runtime_error("Dimensions too large, possible memory allocation failure.");
    }

    std::cout << "\nReading complete. Total points read: " << points.size() << std::endl;
    std::cout << "X: [" << minXCoord << ", " << maxXCoord << "]  -> Total: " << (maxXCoord - minXCoord) << " Meters" << std::endl;
    std::cout << "Y: [" << minYCoord << ", " << maxYCoord << "]  -> Total: " << (maxYCoord - minYCoord) << " Meters" << std::endl;
    std::cout << "Z: [" << minZCoord << ", " << maxZCoord << "]  -> Total: " << (maxZCoord - minZCoord) << " Meters" << std::endl;

    return points;
}

std::vector<uint32_t> LASReader::create_DSM(const std::vector<LASPointData> &points)
{
    // create DSM based on the points read
    std::vector<uint32_t> dsm(static_cast<size_t>(width) * static_cast<size_t>(height), std::numeric_limits<uint32_t>::min());
    std::cout << "Creating DSM with resolution: " << resolution << " meters" << std::endl;
    std::cout << "DSM size: " << width << " x " << height << " pixels" << std::endl;
    int i = 0;
    for (const auto &pt : points)
    {
        if (i % 100000 == 0) // Show progress every 100k points
            // Show progress bar
            std::cout << "\rProcessed " << (&pt - &points[0] + 1) << "/" << points.size() << " points" << std::flush;
        int col = static_cast<int>((pt.x - minX) / resolution);
        int row = static_cast<int>((pt.y - minY) / resolution);

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

void LASReader::refine_DSM(std::vector<uint32_t> &dsm)
{
    // Refine DSM by interpolating missing values
    auto startDSM = std::chrono::high_resolution_clock::now();

    // #pragma omp parallel for collapse(2)
    for (int row = 1; row < height - 1; ++row)
    {
        for (int col = 1; col < width - 1; ++col)
        {
            size_t idx = row * width + col;
            if (dsm[idx] == std::numeric_limits<uint32_t>::min())
            {
                float sum = 0;
                int count = 0;
                for (int dy = -1; dy <= 1; ++dy)
                {
                    for (int dx = -1; dx <= 1; ++dx)
                    {
                        // size_t n_idx = (row + dy) * width + (col + dx);
                        // float val = dsm[n_idx];
                        // if (val != std::numeric_limits<uint32_t>::min()) {
                        //     sum += val;
                        //     count++;
                        // }
                    }
                }
                if (count > 0)
                {
                    dsm[idx] = sum / count; // Mean interpolation
                }
            }
        }
    }

    auto endDSM = std::chrono::high_resolution_clock::now();
    std::cout << "DSM refinement took " << std::chrono::duration<double>(endDSM - startDSM).count() << " seconds.\n";
}