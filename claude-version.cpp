#pragma once

#include <vector>
#include <string>
#include <memory>
#include <chrono>
#include <thread>
#include <future>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <cstring>

// LAS file format structures
#pragma pack(push, 1)
struct LASHeader {
    char fileSignature[4];          // "LASF"
    uint16_t fileSourceId;
    uint16_t globalEncoding;
    uint32_t guidData1;
    uint16_t guidData2;
    uint16_t guidData3;
    uint8_t guidData4[8];
    uint8_t versionMajor;
    uint8_t versionMinor;
    char systemId[32];
    char generatingSoftware[32];
    uint16_t creationDay;
    uint16_t creationYear;
    uint16_t headerSize;
    uint32_t offsetToPointData;
    uint32_t numberOfVariableLengthRecords;
    uint8_t pointDataRecordFormat;
    uint16_t pointDataRecordLength;
    uint32_t numberOfPointRecords;
    uint32_t numberOfPointsByReturn[5];
    double xScaleFactor;
    double yScaleFactor;
    double zScaleFactor;
    double xOffset;
    double yOffset;
    double zOffset;
    double maxX;
    double minX;
    double maxY;
    double minY;
    double maxZ;
    double minZ;
};

struct LASPointFormat0 {
    uint32_t x;
    uint32_t y;
    uint32_t z;
    uint16_t intensity;
    uint8_t returnInfo;
    uint8_t classification;
    int8_t scanAngle;
    uint8_t userData;
    uint16_t pointSourceId;
};

struct LASPointFormat1 : LASPointFormat0 {
    double gpsTime;
};

struct LASPointFormat2 : LASPointFormat0 {
    uint16_t red;
    uint16_t green;
    uint16_t blue;
};

struct LASPointFormat3 : LASPointFormat1 {
    uint16_t red;
    uint16_t green;
    uint16_t blue;
};
#pragma pack(pop)

// Forward declarations
struct Point3D;
struct GridBounds;
class Grid2D;

//=============================================================================
// INTERFACES (Dependency Inversion Principle)
//=============================================================================

class IPointCloudReader {
public:
    virtual ~IPointCloudReader() = default;
    virtual bool open(const std::string& filename) = 0;
    virtual std::vector<Point3D> readAllPoints() = 0;
    virtual GridBounds getBounds() const = 0;
    virtual size_t getPointCount() const = 0;
    virtual void close() = 0;
    virtual std::string getFileInfo() const = 0;
};

class IRasterWriter {
public:
    virtual ~IRasterWriter() = default;
    virtual bool writeRaster(const std::string& filename, 
                           const Grid2D& grid, 
                           const GridBounds& bounds,
                           double resolution) = 0;
};

class IGridProcessor {
public:
    virtual ~IGridProcessor() = default;
    virtual void processPoints(const std::vector<Point3D>& points,
                             Grid2D& grid,
                             const GridBounds& bounds,
                             double resolution) = 0;
};

class IProgressReporter {
public:
    virtual ~IProgressReporter() = default;
    virtual void reportProgress(const std::string& stage, double percentage) = 0;
    virtual void reportTime(const std::string& operation, double seconds) = 0;
};

//=============================================================================
// DATA STRUCTURES
//=============================================================================

struct Point3D {
    double x, y, z;
    
    Point3D() : x(0), y(0), z(0) {}
    Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

struct GridBounds {
    double minX, maxX, minY, maxY;
    
    GridBounds() : minX(0), maxX(0), minY(0), maxY(0) {}
    GridBounds(double minX_, double maxX_, double minY_, double maxY_)
        : minX(minX_), maxX(maxX_), minY(minY_), maxY(maxY_) {}
    
    double width() const { return maxX - minX; }
    double height() const { return maxY - minY; }
};

//=============================================================================
// GRID CLASS (Single Responsibility Principle)
//=============================================================================

class Grid2D {
private:
    std::vector<std::vector<double>> data_;
    size_t width_, height_;
    double noDataValue_;

public:
    Grid2D(size_t width, size_t height, double noDataValue = -9999.0)
        : width_(width), height_(height), noDataValue_(noDataValue) {
        data_.resize(height_, std::vector<double>(width_, noDataValue_));
    }

    void setValue(size_t row, size_t col, double value) {
        if (row < height_ && col < width_) {
            data_[row][col] = value;
        }
    }

    void setMaxValue(size_t row, size_t col, double value) {
        if (row < height_ && col < width_) {
            if (data_[row][col] == noDataValue_ || value > data_[row][col]) {
                data_[row][col] = value;
            }
        }
    }

    double getValue(size_t row, size_t col) const {
        if (row < height_ && col < width_) {
            return data_[row][col];
        }
        return noDataValue_;
    }

    size_t getWidth() const { return width_; }
    size_t getHeight() const { return height_; }
    double getNoDataValue() const { return noDataValue_; }
    
    const std::vector<std::vector<double>>& getData() const { return data_; }
};

//=============================================================================
// CONCRETE IMPLEMENTATIONS
//=============================================================================

// Real LAS File Reader
class LASReader : public IPointCloudReader {
private:
    LASHeader header_;
    std::ifstream file_;
    std::vector<Point3D> points_;
    GridBounds bounds_;
    bool isOpen_;
    std::string filename_;

public:
    LASReader() : isOpen_(false) {}

    bool open(const std::string& filename) override {
        filename_ = filename;
        file_.open(filename, std::ios::binary);
        if (!file_.is_open()) {
            std::cerr << "Failed to open LAS file: " << filename << std::endl;
            return false;
        }

        // Read LAS header
        file_.read(reinterpret_cast<char*>(&header_), sizeof(LASHeader));
        
        // Validate LAS file
        if (std::strncmp(header_.fileSignature, "LASF", 4) != 0) {
            std::cerr << "Not a valid LAS file (missing LASF signature)" << std::endl;
            file_.close();
            return false;
        }

        // Extract bounds from header
        bounds_ = GridBounds(header_.minX, header_.maxX, header_.minY, header_.maxY);
        
        std::cout << "LAS File Info:" << std::endl;
        std::cout << "  Version: " << static_cast<int>(header_.versionMajor) 
                  << "." << static_cast<int>(header_.versionMinor) << std::endl;
        std::cout << "  Point Format: " << static_cast<int>(header_.pointDataRecordFormat) << std::endl;
        std::cout << "  Point Count: " << header_.numberOfPointRecords << std::endl;
        std::cout << "  Point Record Length: " << header_.pointDataRecordLength << " bytes" << std::endl;
        std::cout << "  Bounds: X[" << header_.minX << ", " << header_.maxX << "] "
                  << "Y[" << header_.minY << ", " << header_.maxY << "] "
                  << "Z[" << header_.minZ << ", " << header_.maxZ << "]" << std::endl;
        
        double area = (header_.maxX - header_.minX) * (header_.maxY - header_.minY);
        double density = header_.numberOfPointRecords / area;
        std::cout << "  Point Density: " << std::fixed << std::setprecision(2) 
                  << density << " points/m²" << std::endl;

        isOpen_ = true;
        return true;
    }

    std::vector<Point3D> readAllPoints() override {
        if (!isOpen_) {
            throw std::runtime_error("LAS file not open");
        }

        points_.clear();
        points_.reserve(header_.numberOfPointRecords);

        // Seek to point data
        file_.seekg(header_.offsetToPointData, std::ios::beg);

        // Read points based on format
        switch (header_.pointDataRecordFormat) {
            case 0:
                readPointsFormat0();
                break;
            case 1:
                readPointsFormat1();
                break;
            case 2:
                readPointsFormat2();
                break;
            case 3:
                readPointsFormat3();
                break;
            default:
                throw std::runtime_error("Unsupported point format: " + 
                    std::to_string(header_.pointDataRecordFormat));
        }

        std::cout << "Successfully read " << points_.size() << " points" << std::endl;
        return points_;
    }

    GridBounds getBounds() const override {
        return bounds_;
    }

    size_t getPointCount() const override {
        return header_.numberOfPointRecords;
    }

    std::string getFileInfo() const override {
        if (!isOpen_) return "File not open";
        
        std::ostringstream info;
        info << "LAS " << static_cast<int>(header_.versionMajor) 
             << "." << static_cast<int>(header_.versionMinor)
             << ", Format " << static_cast<int>(header_.pointDataRecordFormat)
             << ", " << header_.numberOfPointRecords << " points";
        return info.str();
    }

    void close() override {
        if (file_.is_open()) {
            file_.close();
        }
        points_.clear();
        isOpen_ = false;
    }

private:
    Point3D convertLASPoint(uint32_t x, uint32_t y, uint32_t z) const {
        double realX = (x * header_.xScaleFactor) + header_.xOffset;
        double realY = (y * header_.yScaleFactor) + header_.yOffset;
        double realZ = (z * header_.zScaleFactor) + header_.zOffset;
        return Point3D(realX, realY, realZ);
    }

    void readPointsFormat0() {
        LASPointFormat0 point;
        for (uint32_t i = 0; i < header_.numberOfPointRecords; ++i) {
            file_.read(reinterpret_cast<char*>(&point), sizeof(LASPointFormat0));
            points_.push_back(convertLASPoint(point.x, point.y, point.z));
            
            // Progress reporting for large files
            if (i % 1000000 == 0 && i > 0) {
                double progress = (static_cast<double>(i) / header_.numberOfPointRecords) * 100.0;
                std::cout << "\rReading points: " << std::fixed << std::setprecision(1) 
                          << progress << "%" << std::flush;
            }
        }
        std::cout << std::endl;
    }

    void readPointsFormat1() {
        LASPointFormat1 point;
        for (uint32_t i = 0; i < header_.numberOfPointRecords; ++i) {
            file_.read(reinterpret_cast<char*>(&point), sizeof(LASPointFormat1));
            points_.push_back(convertLASPoint(point.x, point.y, point.z));
            
            if (i % 1000000 == 0 && i > 0) {
                double progress = (static_cast<double>(i) / header_.numberOfPointRecords) * 100.0;
                std::cout << "\rReading points: " << std::fixed << std::setprecision(1) 
                          << progress << "%" << std::flush;
            }
        }
        std::cout << std::endl;
    }

    void readPointsFormat2() {
        LASPointFormat2 point;
        for (uint32_t i = 0; i < header_.numberOfPointRecords; ++i) {
            file_.read(reinterpret_cast<char*>(&point), sizeof(LASPointFormat2));
            points_.push_back(convertLASPoint(point.x, point.y, point.z));
            
            if (i % 1000000 == 0 && i > 0) {
                double progress = (static_cast<double>(i) / header_.numberOfPointRecords) * 100.0;
                std::cout << "\rReading points: " << std::fixed << std::setprecision(1) 
                          << progress << "%" << std::flush;
            }
        }
        std::cout << std::endl;
    }

    void readPointsFormat3() {
        LASPointFormat3 point;
        for (uint32_t i = 0; i < header_.numberOfPointRecords; ++i) {
            file_.read(reinterpret_cast<char*>(&point), sizeof(LASPointFormat3));
            points_.push_back(convertLASPoint(point.x, point.y, point.z));
            
            if (i % 1000000 == 0 && i > 0) {
                double progress = (static_cast<double>(i) / header_.numberOfPointRecords) * 100.0;
                std::cout << "\rReading points: " << std::fixed << std::setprecision(1) 
                          << progress << "%" << std::flush;
            }
        }
        std::cout << std::endl;
    }
};

// Parallel Grid Processor (Open/Closed Principle - easily extensible)
class ParallelMaxGridProcessor : public IGridProcessor {
private:
    size_t numThreads_;

public:
    explicit ParallelMaxGridProcessor(size_t numThreads = std::thread::hardware_concurrency())
        : numThreads_(numThreads) {}

    void processPoints(const std::vector<Point3D>& points,
                      Grid2D& grid,
                      const GridBounds& bounds,
                      double resolution) override {
        
        const size_t pointsPerThread = points.size() / numThreads_;
        std::vector<std::future<void>> futures;

        // Create thread-local grids to avoid contention
        std::vector<std::unique_ptr<Grid2D>> threadGrids;
        for (size_t i = 0; i < numThreads_; ++i) {
            threadGrids.push_back(std::make_unique<Grid2D>(
                grid.getWidth(), grid.getHeight(), grid.getNoDataValue()));
        }

        // Process points in parallel
        for (size_t threadId = 0; threadId < numThreads_; ++threadId) {
            size_t startIdx = threadId * pointsPerThread;
            size_t endIdx = (threadId == numThreads_ - 1) ? 
                           points.size() : (threadId + 1) * pointsPerThread;

            futures.push_back(std::async(std::launch::async, [&, threadId, startIdx, endIdx]() {
                processPointRange(points, *threadGrids[threadId], bounds, 
                                resolution, startIdx, endIdx);
            }));
        }

        // Wait for all threads to complete
        for (auto& future : futures) {
            future.wait();
        }

        // Merge thread-local grids
        mergeGrids(grid, threadGrids);
    }

private:
    void processPointRange(const std::vector<Point3D>& points,
                          Grid2D& grid,
                          const GridBounds& bounds,
                          double resolution,
                          size_t startIdx,
                          size_t endIdx) {
        
        for (size_t i = startIdx; i < endIdx; ++i) {
            const auto& point = points[i];
            
            // Calculate grid coordinates
            int col = static_cast<int>((point.x - bounds.minX) / resolution);
            int row = static_cast<int>((bounds.maxY - point.y) / resolution);
            
            // Set maximum elevation value
            grid.setMaxValue(static_cast<size_t>(row), static_cast<size_t>(col), point.z);
        }
    }

    void mergeGrids(Grid2D& mainGrid, const std::vector<std::unique_ptr<Grid2D>>& threadGrids) {
        for (size_t row = 0; row < mainGrid.getHeight(); ++row) {
            for (size_t col = 0; col < mainGrid.getWidth(); ++col) {
                double maxValue = mainGrid.getNoDataValue();
                
                for (const auto& threadGrid : threadGrids) {
                    double value = threadGrid->getValue(row, col);
                    if (value != threadGrid->getNoDataValue()) {
                        if (maxValue == mainGrid.getNoDataValue() || value > maxValue) {
                            maxValue = value;
                        }
                    }
                }
                
                if (maxValue != mainGrid.getNoDataValue()) {
                    mainGrid.setValue(row, col, maxValue);
                }
            }
        }
    }
};

// Enhanced GeoTIFF Writer with proper georeferencing
class GeoTIFFWriter : public IRasterWriter {
public:
    bool writeRaster(const std::string& filename,
                    const Grid2D& grid,
                    const GridBounds& bounds,
                    double resolution) override {
        
        // Write as ASCII Grid format (easily readable, can be converted to GeoTIFF)
        std::ofstream file(filename + ".asc");
        if (!file.is_open()) {
            std::cerr << "Failed to create output file: " << filename << ".asc" << std::endl;
            return false;
        }

        // Write ESRI ASCII Grid header
        file << "ncols " << grid.getWidth() << "\n";
        file << "nrows " << grid.getHeight() << "\n";
        file << "xllcorner " << std::fixed << std::setprecision(6) << bounds.minX << "\n";
        file << "yllcorner " << std::fixed << std::setprecision(6) << bounds.minY << "\n";
        file << "cellsize " << std::fixed << std::setprecision(6) << resolution << "\n";
        file << "NODATA_value " << grid.getNoDataValue() << "\n";

        // Write grid data
        const auto& data = grid.getData();
        for (size_t row = 0; row < grid.getHeight(); ++row) {
            for (size_t col = 0; col < grid.getWidth(); ++col) {
                file << std::fixed << std::setprecision(3) << data[row][col];
                if (col < grid.getWidth() - 1) file << " ";
            }
            file << "\n";
        }

        file.close();

        // Also write a world file for georeferencing
        std::ofstream worldFile(filename + ".wld");
        if (worldFile.is_open()) {
            worldFile << std::fixed << std::setprecision(6) << resolution << "\n";  // pixel width
            worldFile << "0.0\n";                                                   // rotation
            worldFile << "0.0\n";                                                   // rotation
            worldFile << std::fixed << std::setprecision(6) << -resolution << "\n"; // pixel height (negative)
            worldFile << std::fixed << std::setprecision(6) << bounds.minX + (resolution / 2.0) << "\n"; // x of upper left
            worldFile << std::fixed << std::setprecision(6) << bounds.maxY - (resolution / 2.0) << "\n"; // y of upper left
            worldFile.close();
        }

        std::cout << "DSM written to: " << filename << ".asc" << std::endl;
        std::cout << "World file written to: " << filename << ".wld" << std::endl;
        std::cout << "Grid dimensions: " << grid.getWidth() << " x " << grid.getHeight() << std::endl;
        std::cout << "Resolution: " << resolution << " meters" << std::endl;
        
        return true;
    }
};

// Console Progress Reporter
class ConsoleProgressReporter : public IProgressReporter {
public:
    void reportProgress(const std::string& stage, double percentage) override {
        std::cout << stage << ": " << std::fixed << std::setprecision(1) 
                  << percentage << "%" << std::endl;
    }

    void reportTime(const std::string& operation, double seconds) override {
        std::cout << operation << " completed in " << std::fixed 
                  << std::setprecision(2) << seconds << " seconds" << std::endl;
    }
};

//=============================================================================
// MAIN DSM GENERATOR CLASS (Single Responsibility + Dependency Injection)
//=============================================================================

class DSMGenerator {
private:
    std::unique_ptr<IPointCloudReader> reader_;
    std::unique_ptr<IRasterWriter> writer_;
    std::unique_ptr<IGridProcessor> processor_;
    std::unique_ptr<IProgressReporter> reporter_;

public:
    DSMGenerator(std::unique_ptr<IPointCloudReader> reader,
                std::unique_ptr<IRasterWriter> writer,
                std::unique_ptr<IGridProcessor> processor,
                std::unique_ptr<IProgressReporter> reporter)
        : reader_(std::move(reader))
        , writer_(std::move(writer))
        , processor_(std::move(processor))
        , reporter_(std::move(reporter)) {}

    bool generateDSM(const std::string& inputFile,
                    const std::string& outputFile,
                    double resolution = 1.0) {
        
        auto startTime = std::chrono::high_resolution_clock::now();

        try {
            // Step 1: Open and read point cloud
            reporter_->reportProgress("Opening file", 0.0);
            if (!reader_->open(inputFile)) {
                std::cerr << "Failed to open input file: " << inputFile << std::endl;
                return false;
            }

            // Get file information
            std::cout << "\n" << reader_->getFileInfo() << std::endl;
            
            // Calculate optimal resolution if not specified
            if (resolution <= 0.0) {
                resolution = calculateOptimalResolution();
                std::cout << "Auto-calculated resolution: " << resolution << " meters" << std::endl;
            }

            auto readStart = std::chrono::high_resolution_clock::now();
            auto points = reader_->readAllPoints();
            auto readEnd = std::chrono::high_resolution_clock::now();
            
            auto readDuration = std::chrono::duration<double>(readEnd - readStart).count();
            reporter_->reportTime("Point cloud reading", readDuration);
            reporter_->reportProgress("Reading points", 25.0);

            if (points.empty()) {
                std::cerr << "No points found in the file" << std::endl;
                return false;
            }

            // Step 2: Create grid
            GridBounds bounds = reader_->getBounds();
            size_t gridWidth = static_cast<size_t>(std::ceil(bounds.width() / resolution)) + 1;
            size_t gridHeight = static_cast<size_t>(std::ceil(bounds.height() / resolution)) + 1;
            
            Grid2D grid(gridWidth, gridHeight);
            reporter_->reportProgress("Grid initialized", 35.0);

            // Step 3: Process points into grid
            auto processStart = std::chrono::high_resolution_clock::now();
            processor_->processPoints(points, grid, bounds, resolution);
            auto processEnd = std::chrono::high_resolution_clock::now();
            
            auto processDuration = std::chrono::duration<double>(processEnd - processStart).count();
            reporter_->reportTime("Grid processing", processDuration);
            reporter_->reportProgress("Processing complete", 80.0);

            // Step 4: Write output
            auto writeStart = std::chrono::high_resolution_clock::now();
            if (!writer_->writeRaster(outputFile, grid, bounds, resolution)) {
                std::cerr << "Failed to write output file: " << outputFile << std::endl;
                return false;
            }
            auto writeEnd = std::chrono::high_resolution_clock::now();
            
            auto writeDuration = std::chrono::duration<double>(writeEnd - writeStart).count();
            reporter_->reportTime("File writing", writeDuration);
            reporter_->reportProgress("DSM generation complete", 100.0);

            // Cleanup
            reader_->close();

            auto endTime = std::chrono::high_resolution_clock::now();
            auto totalDuration = std::chrono::duration<double>(endTime - startTime).count();
            reporter_->reportTime("Total processing", totalDuration);

            return true;

        } catch (const std::exception& e) {
            std::cerr << "Error during DSM generation: " << e.what() << std::endl;
            return false;
        }
    }

private:
    double calculateOptimalResolution() {
        // Calculate optimal resolution based on point density
        GridBounds bounds = reader_->getBounds();
        size_t pointCount = reader_->getPointCount();
        
        double area = bounds.width() * bounds.height();
        double density = static_cast<double>(pointCount) / area;
        
        // Aim for roughly 1-4 points per cell
        double optimalResolution = std::sqrt(4.0 / density);
        
        // Round to sensible values
        if (optimalResolution < 0.1) return 0.1;
        else if (optimalResolution < 0.25) return 0.25;
        else if (optimalResolution < 0.5) return 0.5;
        else if (optimalResolution < 1.0) return 1.0;
        else if (optimalResolution < 2.0) return 2.0;
        else if (optimalResolution < 5.0) return 5.0;
        else return 10.0;
    }
};

//=============================================================================
// FACTORY CLASS (Factory Pattern)
//=============================================================================

class DSMGeneratorFactory {
public:
    static std::unique_ptr<DSMGenerator> createFastDSMGenerator() {
        auto reader = std::make_unique<LASReader>();
        auto writer = std::make_unique<GeoTIFFWriter>();
        auto processor = std::make_unique<ParallelMaxGridProcessor>();
        auto reporter = std::make_unique<ConsoleProgressReporter>();

        return std::make_unique<DSMGenerator>(
            std::move(reader),
            std::move(writer), 
            std::move(processor),
            std::move(reporter)
        );
    }
    
    // Alternative factory for different configurations
    static std::unique_ptr<DSMGenerator> createCustomDSMGenerator(
        std::unique_ptr<IPointCloudReader> reader,
        std::unique_ptr<IRasterWriter> writer,
        size_t numThreads = 0) {
        
        if (numThreads == 0) {
            numThreads = std::thread::hardware_concurrency();
        }
        
        auto processor = std::make_unique<ParallelMaxGridProcessor>(numThreads);
        auto reporter = std::make_unique<ConsoleProgressReporter>();

        return std::make_unique<DSMGenerator>(
            std::move(reader),
            std::move(writer), 
            std::move(processor),
            std::move(reporter)
        );
    }
};

//=============================================================================
// USAGE EXAMPLE
//=============================================================================

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input.las> [output_name] [resolution]" << std::endl;
        std::cout << "  input.las:    Path to input LAS file" << std::endl;
        std::cout << "  output_name:  Output filename (optional, default: dsm_output)" << std::endl;
        std::cout << "  resolution:   Grid resolution in meters (optional, auto-calculated if not specified)" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = (argc > 2) ? argv[2] : "dsm_output";
    double resolution = (argc > 3) ? std::atof(argv[3]) : -1.0; // -1.0 means auto-calculate
    
    std::cout << "=== Fast DSM Generator ===" << std::endl;
    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;
    if (resolution > 0) {
        std::cout << "Resolution: " << resolution << " meters" << std::endl;
    } else {
        std::cout << "Resolution: Auto-calculated" << std::endl;
    }
    std::cout << "Threads available: " << std::thread::hardware_concurrency() << std::endl;
    std::cout << std::endl;

    try {
        // Create DSM generator using factory
        auto dsmGenerator = DSMGeneratorFactory::createFastDSMGenerator();
        
        // Generate DSM
        bool success = dsmGenerator->generateDSM(inputFile, outputFile, resolution);
        
        if (success) {
            std::cout << "\n✓ DSM generation completed successfully!" << std::endl;
        } else {
            std::cout << "\n✗ DSM generation failed!" << std::endl;
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}