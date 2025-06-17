#include "hypergraph_ordering.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <stdexcept>

using namespace hypergraph_ordering;

int main(int argc, char *argv[])
{
    try
    {
        std::string matrix_file;
        // Parse command line arguments
        if (argc < 2)
        {
            std::cout << "No matrix file provided" << std::endl;
            return 0;
        }
        else
        {
            matrix_file = argv[1];
        }

        // Load matrix
        std::cout << "\nLoading matrix: " << matrix_file << std::endl;
        auto matrix = SparseMatrix::loadFromMatrixMarket(matrix_file);

        // Print matrix information
        std::cout << "\n=== Matrix Information ===" << std::endl;
        matrix.printInfo();

        // TODO: check if matrix is symmetric

        // Get matrix statistics
        Index min_deg, max_deg;
        double avg_deg, std_deg;
        matrix.getStatistics(min_deg, max_deg, avg_deg, std_deg);
        std::cout << "Degree statistics: min=" << min_deg << ", max=" << max_deg
                  << ", avg=" << std::fixed << std::setprecision(1) << avg_deg
                  << ", std=" << std_deg << std::endl;

        // Configure hypergraph ordering
        OrderingConfig config;
        config.verbose = false;
        config.clique_type = OrderingConfig::C2_CLIQUES;
        // config.max_recursion_depth = 10;  // Allow deeper recursion
        // config.min_subproblem_size = 200; // Continue partitioning longer
        // // config.min_nodes_for_partitioning = 20; // Partition smaller subproblems
        // config.use_minimum_degree = true; // Use AMD for base cases

        config.imbalance = 0.001; // Tight balance constraint (1%)
        // Create ordering object
        std::cout << "\n=== Running Hypergraph Ordering ===" << std::endl;
        HypergraphOrdering ordering(config);

        // Run the algorithm
        auto start_time = std::chrono::high_resolution_clock::now();
        auto result = ordering.orderMatrix(matrix);
        auto end_time = std::chrono::high_resolution_clock::now();

        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "\nTotal execution time: " << total_duration.count() << " ms" << std::endl;

        // TODO: Validate ordering result?

        // Print first and last few vertices in ordering
        std::cout << "\n=== Ordering Result ===" << std::endl;
        std::cout << "First 10 vertices: ";
        for (size_t i = 0; i < std::min(size_t(10), result.size()); ++i)
        {
            std::cout << result[i] << " ";
        }
        std::cout << std::endl;

        if (result.size() > 10)
        {
            std::cout << "Last 10 vertices: ";
            for (size_t i = std::max(size_t(0), result.size() - 10); i < result.size(); ++i)
            {
                std::cout << result[i] << " ";
            }
            std::cout << std::endl;
        }

        // Save ordering to file
        std::string output_file = matrix_file + ".ordering";
        std::ofstream out(output_file);
        if (out.is_open())
        {
            for (Index v : result)
            {
                out << v << std::endl;
            }
            out.close();
            std::cout << "\nOrdering saved to: " << output_file << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}