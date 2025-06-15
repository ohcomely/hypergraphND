#include "hypergraph_ordering.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <stdexcept>

using namespace hypergraph_ordering;

// Function to create a simple test matrix
void createTestMatrix(const std::string &filename, int n = 10)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot create test file: " + filename);
    }

    // Write Matrix Market header
    file << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
    file << "% Simple test matrix for hypergraph ordering" << std::endl;
    file << n << " " << n << " ";

    // Create a simple grid-like structure
    std::vector<std::pair<int, int>> edges;

    // Add diagonal entries
    for (int i = 1; i <= n; ++i)
    {
        edges.push_back({i, i});
    }

    // Add off-diagonal entries (chain + some connections)
    for (int i = 1; i < n; ++i)
    {
        edges.push_back({i, i + 1}); // Chain
    }

    // Add some additional connections for triangles
    if (n >= 3)
    {
        edges.push_back({1, 3}); // Creates triangle (1,2,3)
    }
    if (n >= 5)
    {
        edges.push_back({3, 5}); // Creates triangle (3,4,5)
        edges.push_back({1, 5}); // Additional connection
    }
    if (n >= 7)
    {
        edges.push_back({5, 7}); // Continues pattern
    }

    file << edges.size() << std::endl;

    // Write edges
    for (const auto &edge : edges)
    {
        file << edge.first << " " << edge.second << " 1.0" << std::endl;
    }

    file.close();
    std::cout << "Created test matrix: " << filename << " (" << n << "x" << n
              << ", " << edges.size() << " entries)" << std::endl;
}

// Function to validate ordering quality
void validateOrdering(const SparseMatrix &matrix, const std::vector<Index> &ordering)
{
    const Index n = matrix.rows();

    std::cout << "\n=== Ordering Validation ===" << std::endl;

    // Check completeness
    if (ordering.size() != n)
    {
        std::cout << "❌ Invalid ordering size: " << ordering.size() << " vs " << n << std::endl;
        return;
    }

    // Check permutation
    std::vector<bool> seen(n, false);
    bool valid_permutation = true;
    for (Index v : ordering)
    {
        if (v >= n || seen[v])
        {
            valid_permutation = false;
            break;
        }
        seen[v] = true;
    }

    if (valid_permutation)
    {
        std::cout << "✅ Valid permutation" << std::endl;
    }
    else
    {
        std::cout << "❌ Invalid permutation" << std::endl;
        return;
    }

    // Simple fill-in estimation (very basic)
    std::vector<std::vector<bool>> filled(n, std::vector<bool>(n, false));

    // Mark original nonzeros
    for (Index i = 0; i < n; ++i)
    {
        for (Index ptr = matrix.rowPtr()[i]; ptr < matrix.rowPtr()[i + 1]; ++ptr)
        {
            Index j = matrix.colInd()[ptr];
            filled[i][j] = filled[j][i] = true;
        }
    }

    Index original_nnz = matrix.nnz();
    Index fill_in = 0;

    // Simulate elimination in given order
    std::vector<Index> inv_perm(n);
    for (Index i = 0; i < n; ++i)
    {
        inv_perm[ordering[i]] = i;
    }

    for (Index step = 0; step < n; ++step)
    {
        Index pivot = ordering[step];

        // Find neighbors of pivot that come later in ordering
        std::vector<Index> later_neighbors;
        for (Index ptr = matrix.rowPtr()[pivot]; ptr < matrix.rowPtr()[pivot + 1]; ++ptr)
        {
            Index neighbor = matrix.colInd()[ptr];
            if (inv_perm[neighbor] > step)
            {
                later_neighbors.push_back(neighbor);
            }
        }

        // Add fill-in edges between later neighbors
        for (size_t i = 0; i < later_neighbors.size(); ++i)
        {
            for (size_t j = i + 1; j < later_neighbors.size(); ++j)
            {
                Index u = later_neighbors[i];
                Index v = later_neighbors[j];
                if (!filled[u][v])
                {
                    filled[u][v] = filled[v][u] = true;
                    fill_in++;
                }
            }
        }
    }

    std::cout << "📊 Fill-in analysis:" << std::endl;
    std::cout << "   Original NNZ: " << original_nnz << std::endl;
    std::cout << "   Fill-in: " << fill_in << std::endl;
    std::cout << "   Total NNZ after elimination: " << (original_nnz + fill_in) << std::endl;
    std::cout << "   Fill-in ratio: " << std::fixed << std::setprecision(2)
              << (100.0 * fill_in / original_nnz) << "%" << std::endl;
}

// Function to test different configurations
void testConfigurations(const SparseMatrix &matrix)
{
    std::cout << "\n=== Testing Different Configurations ===" << std::endl;

    // Test C2 configuration
    {
        std::cout << "\n--- C2 (Edge-based) Configuration ---" << std::endl;
        OrderingConfig config;
        config.clique_type = OrderingConfig::C2_CLIQUES;
        config.verbose = false; // Reduce output for comparison

        HypergraphOrdering ordering_c2(config);

        auto start = std::chrono::high_resolution_clock::now();
        auto result_c2 = ordering_c2.orderMatrix(matrix);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "C2 time: " << duration.count() << " ms" << std::endl;

        const auto &stats = ordering_c2.getStatistics();
        std::cout << "C2 stats: " << stats.hypergraph_partitions << " partitions, "
                  << stats.minimum_degree_calls << " MD calls, depth "
                  << stats.max_recursion_depth_reached << std::endl;
    }

    // Test C3 configuration
    {
        std::cout << "\n--- C3 (Triangle-based) Configuration ---" << std::endl;
        OrderingConfig config;
        config.clique_type = OrderingConfig::C3_CLIQUES;
        config.verbose = false;

        HypergraphOrdering ordering_c3(config);

        auto start = std::chrono::high_resolution_clock::now();
        auto result_c3 = ordering_c3.orderMatrix(matrix);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "C3 time: " << duration.count() << " ms" << std::endl;

        const auto &stats = ordering_c3.getStatistics();
        std::cout << "C3 stats: " << stats.hypergraph_partitions << " partitions, "
                  << stats.minimum_degree_calls << " MD calls, depth "
                  << stats.max_recursion_depth_reached << std::endl;
    }
}

int main(int argc, char *argv[])
{
    try
    {
        std::cout << "=== Hypergraph-based Matrix Ordering Test ===" << std::endl;

        std::string matrix_file;
        bool use_test_matrix = false;

        // Parse command line arguments
        if (argc < 2)
        {
            std::cout << "No matrix file provided, creating test matrix..." << std::endl;
            matrix_file = "test_matrix.mtx";
            use_test_matrix = true;
        }
        else
        {
            matrix_file = argv[1];
        }

        // Create test matrix if needed
        if (use_test_matrix)
        {
            createTestMatrix(matrix_file, 20); // 20x20 test matrix
        }

        // Load matrix
        std::cout << "\nLoading matrix: " << matrix_file << std::endl;
        auto matrix = SparseMatrix::loadFromMatrixMarket(matrix_file);

        // Print matrix information
        std::cout << "\n=== Matrix Information ===" << std::endl;
        matrix.printInfo();

        // // Check if matrix is symmetric
        // if (matrix.isSymmetric())
        // {
        //     std::cout << "✅ Matrix is symmetric" << std::endl;
        // }
        // else
        // {
        //     std::cout << "⚠️  Matrix is not perfectly symmetric" << std::endl;
        // }

        // Get matrix statistics
        Index min_deg, max_deg;
        double avg_deg, std_deg;
        matrix.getStatistics(min_deg, max_deg, avg_deg, std_deg);
        std::cout << "Degree statistics: min=" << min_deg << ", max=" << max_deg
                  << ", avg=" << std::fixed << std::setprecision(1) << avg_deg
                  << ", std=" << std_deg << std::endl;

        // Configure hypergraph ordering
        OrderingConfig config;
        config.verbose = true;
        config.clique_type = OrderingConfig::C2_CLIQUES;
        config.min_subproblem_size = 5; // Small for demo
        config.max_recursion_depth = 8;

        // Test if KaHyPar config exists
        // if (!utils::fileExists(config.kahypar_config_path))
        // {
        //     std::cout << "\n⚠️  KaHyPar config file not found: " << config.kahypar_config_path << std::endl;
        //     std::cout << "Please provide a valid KaHyPar configuration file." << std::endl;
        //     std::cout << "You can download example configs from KaHyPar repository." << std::endl;
        //     return 1;
        // }

        // Create ordering object
        std::cout << "\n=== Running Hypergraph Ordering ===" << std::endl;
        HypergraphOrdering ordering(config);

        // Run the algorithm
        auto start_time = std::chrono::high_resolution_clock::now();
        auto result = ordering.orderMatrix(matrix);
        auto end_time = std::chrono::high_resolution_clock::now();

        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "\nTotal execution time: " << total_duration.count() << " ms" << std::endl;

        // Validate the result
        validateOrdering(matrix, result);

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

        // Test different configurations if matrix is not too large
        if (matrix.rows() <= 100)
        {
            testConfigurations(matrix);
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
            std::cout << "\n📁 Ordering saved to: " << output_file << std::endl;
        }

        std::cout << "\n✅ Test completed successfully!" << std::endl;

        // Clean up test file
        if (use_test_matrix)
        {
            std::remove(matrix_file.c_str());
            std::cout << "Test matrix file cleaned up." << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "❌ Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}