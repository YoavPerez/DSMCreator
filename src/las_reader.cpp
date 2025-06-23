#include "las_reader.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <tuple>
#include <cmath>
#include <iomanip>
#define BATCH_SIZE 1'000'000
#define REFINMENT_TIMES 5

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

std::vector<double> LASReader::create_DSM(const std::vector<LASPointData> &points)
{
    // create DSM based on the points read
    std::vector<double> dsm(static_cast<size_t>(width) * static_cast<size_t>(height), std::numeric_limits<double>::min());
    std::vector<LASPointData> selected_points(static_cast<size_t>(width) * static_cast<size_t>(height), LASPointData{0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0, 0, 0});
    std::cout << "Creating DSM with resolution: " << resolution << " meters" << std::endl;
    std::cout << "DSM size: " << width << " x " << height << " pixels" << std::endl;
    int i = 0;
    for (const auto &pt : points)
    {
        if (i % 1'000'000 == 0) // Show progress every 1M points
            // Show progress bar
            std::cout << "\rProcessed " << (&pt - &points[0] + 1) << "/" << points.size() << " points" << std::flush;
        
        int col = static_cast<int>(std::round((pt.x * header.scaleX + header.offsetX - minXCoord) / resolution));
        int row = static_cast<int>(std::round((maxYCoord - (pt.y * header.scaleY + header.offsetY)) / resolution));

        if (col >= 0 && col < width && row >= 0 && row < height)
        {
            size_t idx = row * width + col;
            double z = static_cast<double>(pt.z) * header.scaleZ + header.offsetZ;
            dsm[idx] = std::max(dsm[idx], z);
            selected_points[idx] = pt; // Store the point data for potential further use
        }
        i++;
        // show progress bar
    }
    std::cout << "\rProcessed " << points.size() << "/" << points.size() << " points" << std::endl;
    std::cout << "DSM sample values: " << std::endl;

    for(int i=0; i < REFINMENT_TIMES; i++){
        std::cout << "Refining DSM, iteration: " << i + 1 << std::endl;
        refine_DSM(dsm, selected_points);
    }
    return dsm;
}

std::tuple<int, double> count_non_null_in_range(const std::vector<double> &dsm, int row, int col, int range, int width)
{
    int count = 0;
    double sum = 0.0;
    int height = static_cast<int>(dsm.size() / width);
    for (int r = -range; r <= range; ++r)
    {
        for (int c = -range; c <= range; ++c)
        {
            if (r == 0 && c == 0) continue; // Skip the center
            int nr = row + r;
            int nc = col + c;
            if (nr >= 0 && nr < height && nc >= 0 && nc < width)
            {
                size_t nidx = nr * width + nc;
                if (dsm[nidx] != std::numeric_limits<double>::min())
                {
                    count++;
                    sum += dsm[nidx];
                }
            }
        }
    }
    return std::make_tuple(count, sum);
}

double calculate_stddev(const std::vector<double> &dsm, int row, int col, int range, int width, double mean)
{
    double sum = 0.0;
    int count = 0;
    int height = static_cast<int>(dsm.size() / width);
    for (int r = -range; r <= range; ++r)
    {
        for (int c = -range; c <= range; ++c)
        {
            if (r == 0 && c == 0) continue; // Skip the center
            int nr = row + r;
            int nc = col + c;
            if (nr >= 0 && nr < height && nc >= 0 && nc < width)
            {
                size_t nidx = nr * width + nc;
                if (dsm[nidx] != std::numeric_limits<double>::min())
                {
                    double diff = dsm[nidx] - mean;
                    sum += diff * diff;
                    count++;
                }
            }
        }
    }
    return std::sqrt(sum / count);
}

void LASReader::refine_DSM(std::vector<double> &dsm, std::vector<LASPointData> &selected_points)
{
    // Refine DSM by interpolating missing values
    auto startDSM = std::chrono::high_resolution_clock::now();
    // interpolate miising values, if the standard deviation of the surrounding points is less than 0.5, use the mean of the surrounding points
    std::cout << "Refining DSM..." << std::endl;
    int fixed_points = 0;
    int null_points = 0;
    for (int row = 0; row < height; ++row)
    {
        for (int col = 0; col < width; ++col)
        {
            size_t idx = row * width + col;
            if (dsm[idx] == std::numeric_limits<double>::min())
            {
                null_points++;
                double meter_range = 1.5; // 1.5 meters range for interpolation
                int range = static_cast<int>(std::round(meter_range / resolution)); // Range in pixels, 1.5 meters
                auto [count, sum] = count_non_null_in_range(dsm, row, col, range, width);
                // If we have enough surrounding points, calculate the mean and standard deviation
                if (count > range*range /2.0)
                {
                    double mean = sum / count;
                    double stddev = calculate_stddev(dsm, row, col, range, width, mean);
                    if (stddev < meter_range) // If standard deviation is less than resolution, use the mean
                    {
                        dsm[idx] = mean;
                        selected_points[idx].z = static_cast<int32_t>(std::round(mean / header.scaleZ - header.offsetZ));
                        fixed_points++;
                    }
                    else
                    {
                        // If standard deviation is too high, keep the value as is (or set to a specific value if needed)
                        dsm[idx] = std::numeric_limits<double>::min(); // or somakeme other logic
                    }
                }
            }
        }
    }
    auto endDSM = std::chrono::high_resolution_clock::now();
    std::cout << "DSM refinement took " << std::chrono::duration<double>(endDSM - startDSM).count() << " seconds.\n";
    std::cout << "Fixed " << fixed_points << " points, null points: " << null_points << std::endl;
}