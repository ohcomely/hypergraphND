#include <hypergraph_ordering.hpp>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace hypergraph_ordering
{
    bool OrderingConfig::isValid() const
    {
        struct stat buffer;
        if (stat(kahypar_config_path.c_str(), &buffer) != 0)
        {
            std::cerr << "KaHyPar configuration file not found: " << kahypar_config_path << std::endl;
            return false;
        }

        if (max_recursion_depth <= 0 || min_subproblem_size <= 0 || min_nodes_for_partitioning <= 0 || num_blocks != 2)
        {
            std::cerr << "Invalid configuration parameters: "
                      << "max_recursion_depth=" << max_recursion_depth
                      << ", min_subproblem_size=" << min_subproblem_size
                      << ", min_nodes_for_partitioning=" << min_nodes_for_partitioning
                      << ", num_blocks=" << num_blocks << std::endl;
            return false;
        }

        return true;
    }

    void HypergraphOrdering::recursiveNestedDissection(
        const SparseMatrix &matrix,
        const std::vector<Index> &vertices,
        Index depth,
        std::vector<Index> &ordering) const
    {

        Timer timer;

        if (config_.verbose && depth == 0)
        {
            std::cout << "Starting recursive nested dissection..." << std::endl;
        }

        const Index n = vertices.size();

        // Update statistics
        stats_.total_subproblems++;
        if (depth > stats_.max_recursion_depth_reached)
        {
            stats_.max_recursion_depth_reached = depth;
        }

        if (config_.verbose)
        {
            std::cout << "Depth " << depth << ": Processing " << n
                      << " vertices" << std::endl;
        }

        // Base case 1: Small subproblem - use minimum degree
        if (n <= config_.min_subproblem_size || depth >= config_.max_recursion_depth)
        {
            if (config_.verbose)
            {
                std::cout << "Base case reached (size=" << n << ", depth=" << depth
                          << "), using minimum degree ordering" << std::endl;
            }

            auto md_ordering = minimumDegreeOrdering(matrix, vertices);
            ordering.insert(ordering.end(), md_ordering.begin(), md_ordering.end());
            return;
        }

        // // Base case 2: Too few nodes for meaningful partitioning
        // if (n < config_.min_nodes_for_partitioning)
        // {
        //     if (config_.verbose)
        //     {
        //         std::cout << "Too few vertices for partitioning, using minimum degree" << std::endl;
        //     }

        //     auto md_ordering = minimumDegreeOrdering(matrix, vertices);
        //     ordering.insert(ordering.end(), md_ordering.begin(), md_ordering.end());
        //     return;
        // }

        try
        {
            // Step 1: Extract submatrix for current vertices
            std::unordered_map<Index, Index> vertex_map;
            SparseMatrix submatrix;

            if (n == matrix.rows())
            {
                // Working with full matrix
                submatrix = matrix;
                for (Index i = 0; i < n; ++i)
                {
                    vertex_map[vertices[i]] = i;
                }
            }
            else
            {
                // Extract submatrix
                submatrix = matrix.extractSubmatrix(vertices, vertex_map);
            }

            if (config_.verbose)
            {
                std::cout << "Submatrix extracted: " << submatrix.rows()
                          << "x" << submatrix.cols() << ", nnz=" << submatrix.nnz() << std::endl;
            }

            // Step 2: Construct hypergraph for submatrix
            auto hg_data = constructCliqueNodeHypergraph(submatrix);

            if (hg_data.num_nodes == 0 || hg_data.num_nets == 0)
            {
                if (config_.verbose)
                {
                    std::cout << "Empty hypergraph, using minimum degree ordering" << std::endl;
                }
                auto md_ordering = minimumDegreeOrdering(matrix, vertices);
                ordering.insert(ordering.end(), md_ordering.begin(), md_ordering.end());
                return;
            }

            if (config_.verbose)
            {
                std::cout << "Hypergraph constructed: " << hg_data.num_nodes
                          << " nodes, " << hg_data.num_nets << " nets" << std::endl;
            }

            // Step 3: Partition the hypergraph
            auto partition = partitionHypergraph(hg_data);

            // Step 4: Get edge information for decoding (C2 case)
            std::vector<Index> edge_rows, edge_cols;
            if (config_.clique_type == OrderingConfig::C2_CLIQUES)
            {
                auto upper_tri = submatrix.upperTriangular();

                for (Index i = 0; i < submatrix.rows(); ++i)
                {
                    for (Index ptr = upper_tri.rowPtr()[i]; ptr < upper_tri.rowPtr()[i + 1]; ++ptr)
                    {
                        Index j = upper_tri.colInd()[ptr];
                        if (j > i)
                        {
                            edge_rows.push_back(i);
                            edge_cols.push_back(j);
                        }
                    }
                }
            }

            // Step 5: Decode partition to vertex separator
            auto separator_result = decodePartitionToVertexSeparator(
                submatrix, partition, edge_rows, edge_cols);

            // Step 6: Validate separator quality
            // if (!validateVertexSeparator(submatrix, separator_result))
            // {
            //     if (config_.verbose)
            //     {
            //         std::cout << "Invalid separator, falling back to minimum degree" << std::endl;
            //     }
            //     auto md_ordering = minimumDegreeOrdering(matrix, vertices);
            //     ordering.insert(ordering.end(), md_ordering.begin(), md_ordering.end());
            //     return;
            // }

            // Step 7: Map separator result back to original vertex indices
            std::vector<Index> orig_separator, orig_part1, orig_part2;

            for (Index local_v : separator_result.separator)
            {
                orig_separator.push_back(vertices[local_v]);
            }
            for (Index local_v : separator_result.part1)
            {
                orig_part1.push_back(vertices[local_v]);
            }
            for (Index local_v : separator_result.part2)
            {
                orig_part2.push_back(vertices[local_v]);
            }

            if (config_.verbose)
            {
                std::cout << "Separator found: " << orig_separator.size()
                          << " separator, " << orig_part1.size() << " + "
                          << orig_part2.size() << " parts" << std::endl;
            }

            // Step 8: Check for degenerate partitions
            if (orig_part1.empty() || orig_part2.empty())
            {
                if (config_.verbose)
                {
                    std::cout << "Degenerate partition, using minimum degree" << std::endl;
                }
                auto md_ordering = minimumDegreeOrdering(matrix, vertices);
                ordering.insert(ordering.end(), md_ordering.begin(), md_ordering.end());
                return;
            }

            // Step 9: Recursive calls on parts (nested dissection order)
            if (!orig_part1.empty())
            {
                recursiveNestedDissection(matrix, orig_part1, depth + 1, ordering);
            }

            if (!orig_part2.empty())
            {
                recursiveNestedDissection(matrix, orig_part2, depth + 1, ordering);
            }

            // Step 10: Add separator vertices last (nested dissection property)
            for (Index v : orig_separator)
            {
                ordering.push_back(v);
            }

            if (config_.verbose)
            {
                std::cout << "Depth " << depth << " completed in "
                          << timer.elapsed() << "s" << std::endl;
            }
        }
        catch (const std::exception &e)
        {
            if (config_.verbose)
            {
                std::cout << "Exception in nested dissection: " << e.what()
                          << ", falling back to minimum degree" << std::endl;
            }

            // Fallback to minimum degree on any error
            auto md_ordering = minimumDegreeOrdering(matrix, vertices);
            ordering.insert(ordering.end(), md_ordering.begin(), md_ordering.end());
        }
    }

    // Main entry point methods
    std::vector<Index> HypergraphOrdering::orderMatrix(const SparseMatrix &matrix)
    {
        Timer total_timer;

        if (config_.verbose)
        {
            std::cout << "=== Hypergraph-based Matrix Ordering ===" << std::endl;
            std::cout << "Matrix: " << matrix.rows() << "x" << matrix.cols()
                      << ", nnz=" << matrix.nnz() << std::endl;
            std::cout << "Clique type: " << (config_.clique_type == OrderingConfig::C2_CLIQUES ? "C2" : "C3") << std::endl;
        }

        // Reset statistics
        stats_ = Statistics{};

        // Validate input
        if (!matrix.isValid())
        {
            throw std::invalid_argument("Invalid sparse matrix");
        }

        if (matrix.rows() == 0)
        {
            return {};
        }

        // Create initial vertex list
        std::vector<Index> vertices;
        vertices.reserve(matrix.rows());
        for (Index i = 0; i < matrix.rows(); ++i)
        {
            vertices.push_back(i);
        }

        // Main algorithm
        std::vector<Index> ordering;
        ordering.reserve(matrix.rows());

        recursiveNestedDissection(matrix, vertices, 0, ordering);

        // Update final statistics
        stats_.total_time = total_timer.elapsed();

        if (config_.verbose)
        {
            std::cout << "=== Ordering Complete ===" << std::endl;
            std::cout << "Total time: " << stats_.total_time << "s" << std::endl;
            std::cout << "Hypergraph construction: " << stats_.hypergraph_construction_time << "s" << std::endl;
            std::cout << "Partitioning: " << stats_.partitioning_time << "s" << std::endl;
            std::cout << "Separator decoding: " << stats_.separator_decode_time << "s" << std::endl;
            std::cout << "Max recursion depth: " << stats_.max_recursion_depth_reached << std::endl;
            std::cout << "Total subproblems: " << stats_.total_subproblems << std::endl;
            std::cout << "Hypergraph partitions: " << stats_.hypergraph_partitions << std::endl;
            std::cout << "Minimum degree calls: " << stats_.minimum_degree_calls << std::endl;
        }

        // Validate result
        if (!utils::isValidOrdering(ordering, matrix.rows()))
        {
            throw std::runtime_error("Generated invalid ordering");
        }

        return ordering;
    }

    std::vector<Index> HypergraphOrdering::orderMatrix(const std::string &matrix_file)
    {
        if (config_.verbose)
        {
            std::cout << "Loading matrix from: " << matrix_file << std::endl;
        }

        auto matrix = SparseMatrix::loadFromMatrixMarket(matrix_file);

        if (config_.verbose)
        {
            std::cout << "Matrix loaded successfully" << std::endl;
        }

        return orderMatrix(matrix);
    }

    // Helper method to determine if hypergraph partitioning should be used
    bool HypergraphOrdering::shouldUseHypergraphPartitioning(
        const SparseMatrix &matrix,
        const std::vector<Index> &vertices,
        Index depth) const
    {

        const Index n = vertices.size();

        // Size-based criteria
        // if (n < config_.min_nodes_for_partitioning)
        // {
        //     return false;
        // }

        if (n <= config_.min_subproblem_size)
        {
            return false;
        }

        // Depth-based criteria
        if (depth >= config_.max_recursion_depth)
        {
            return false;
        }

        // Density-based criteria (optional)
        // Could add checks for matrix density, connectivity, etc.

        return true;
    }

    bool HypergraphOrdering::validateVertexSeparator(
        const SparseMatrix &matrix,
        const VertexSeparatorResult &result) const
    {

        const Index n = matrix.rows();

        // Check that all vertices are accounted for
        std::unordered_set<Index> all_vertices;
        for (Index v : result.separator)
            all_vertices.insert(v);
        for (Index v : result.part1)
            all_vertices.insert(v);
        for (Index v : result.part2)
            all_vertices.insert(v);

        if (all_vertices.size() != n)
        {
            if (config_.verbose)
            {
                std::cout << "Warning: Separator validation failed - vertex count mismatch" << std::endl;
            }
            return false;
        }

        return true;
    }
}