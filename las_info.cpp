// las_info.cpp
// Simple LAS file parser: prints header info and first 10 points
// Requires PDAL

#include <iostream>
#include <pdal/StageFactory.hpp>
#include <pdal/Options.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/Dimension.hpp>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.las>\n";
        return 1;
    }
    std::string las_file = argv[1];
    pdal::StageFactory factory;
    pdal::Options las_opts;
    las_opts.add("filename", las_file);
    pdal::Stage* reader = factory.createStage("readers.las");
    if (!reader) {
        std::cerr << "PDAL LAS reader not found!\n";
        return 2;
    }

    auto start = std::chrono::high_resolution_clock::now();

    reader->setOptions(las_opts);
    pdal::PointTable table;
    reader->prepare(table);
    pdal::PointViewSet viewSet = reader->execute(table);

    // Combine all PointViews into one for faster sequential access
    pdal::PointViewPtr mergedView = std::make_shared<pdal::PointView>(table);
    for (auto& v : viewSet)
        mergedView->append(*v);

    // Map to store the point with largest Z for each (X, Y)
    std::map<std::pair<double, double>, pdal::PointId> maxZPoints;
    for (pdal::PointId i = 0; i < mergedView->size(); ++i) {
        double x = mergedView->getFieldAs<double>(pdal::Dimension::Id::X, i);
        double y = mergedView->getFieldAs<double>(pdal::Dimension::Id::Y, i);
        double z = mergedView->getFieldAs<double>(pdal::Dimension::Id::Z, i);
        auto key = std::make_pair(x, y);
        auto it = maxZPoints.find(key);
        if (it == maxZPoints.end() || z > mergedView->getFieldAs<double>(pdal::Dimension::Id::Z, it->second)) {
            maxZPoints[key] = i;
        }
    }

    // Create a new PointView with only the max Z points
    pdal::PointViewPtr view = std::make_shared<pdal::PointView>(table);
    for (const auto& kv : maxZPoints) {
        view->appendPoint(*mergedView, kv.second);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Processing time: " << elapsed.count() << " seconds\n";

    std::cout << "LAS file: " << las_file << "\n";
    std::cout << "Number of points: " << view->size() << "\n";
    std::cout << "First 10 points (X, Y, Z):\n";
    for (pdal::PointId i = 0; i < std::min<pdal::PointId>(10, view->size()); ++i) {
        double x = view->getFieldAs<double>(pdal::Dimension::Id::X, i);
        double y = view->getFieldAs<double>(pdal::Dimension::Id::Y, i);
        double z = view->getFieldAs<double>(pdal::Dimension::Id::Z, i);
        std::cout << i << ": " << x << ", " << y << ", " << z << "\n";
    }
    return 0;
}
