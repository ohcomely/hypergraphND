#include "hypergraph_ordering.hpp"
#include <iomanip>

namespace hypergraph_ordering
{
    VertexSeparatorResult HypergraphOrdering::decodePartitionToVertexSeparator(
        const SparseMatrix &matrix,
        const std::vector<PartitionID> &partition,
        const std::vector<Index> &edge_rows,
        const std::vector<Index> &edge_cols) const
    {

        Timer timer;

        if (config_.verbose)
        {
            std::cout << "Decoding hypergraph partition to vertex separator..." << std::endl;
        }

        if (config_.clique_type == OrderingConfig::C2_CLIQUES)
        {
            auto result = decodeC2PartitionToVertexSeparator(matrix, partition, edge_rows, edge_cols);
            stats_.separator_decode_time += timer.elapsed();
            return result;
        }
        else
        {
            // For C3, we need to extract triangle information
            auto result = decodeC3PartitionToVertexSeparator(matrix, partition);
            stats_.separator_decode_time += timer.elapsed();
            return result;
        }
    }
    VertexSeparatorResult HypergraphOrdering::decodeC2PartitionToVertexSeparator(
        const SparseMatrix &matrix,
        const std::vector<PartitionID> &partition,
        const std::vector<Index> &edge_rows,
        const std::vector<Index> &edge_cols) const
    {
        Timer timer;

        const Index n = matrix.rows();
        const Index num_edges = edge_rows.size();

        if (config_.verbose)
        {
            std::cout << "Decoding C2 clique-node hypergraph partition..." << std::endl;
            std::cout << "  Matrix size: " << n << "x" << n << std::endl;
            std::cout << "  Number of edges: " << num_edges << std::endl;
        }

        // Validate input
        if (partition.size() != num_edges || edge_cols.size() != num_edges)
        {
            throw std::runtime_error("Partition size mismatch with number of edges");
        }

        // Step 1: Classify edges by partition
        std::vector<std::unordered_set<Index>> part_edges(2);
        std::vector<Index> cut_edges; // Edges that cross partitions

        for (Index edge_idx = 0; edge_idx < num_edges; ++edge_idx)
        {
            PartitionID part = partition[edge_idx];
            if (part == 0)
            {
                part_edges[0].insert(edge_idx);
            }
            else if (part == 1)
            {
                part_edges[1].insert(edge_idx);
            }
            else
            {
                // Edge is in cut (if KaHyPar uses special partition ID for cut)
                cut_edges.push_back(edge_idx);
            }
        }

        if (config_.verbose)
        {
            std::cout << "  Partition 0 edges: " << part_edges[0].size() << std::endl;
            std::cout << "  Partition 1 edges: " << part_edges[1].size() << std::endl;
            std::cout << "  Cut edges: " << cut_edges.size() << std::endl;
        }

        // Step 2: Build vertex-to-edges mapping for each partition
        std::vector<std::unordered_set<Index>> vertex_part_edges(n);
        std::vector<std::unordered_set<PartitionID>> vertex_partitions(n);

        // Track which partitions each vertex is connected to
        for (Index edge_idx = 0; edge_idx < num_edges; ++edge_idx)
        {
            Index u = edge_rows[edge_idx];
            Index v = edge_cols[edge_idx];
            PartitionID part = partition[edge_idx];

            if (part < 2)
            { // Valid partition (0 or 1)
                vertex_part_edges[u].insert(edge_idx);
                vertex_part_edges[v].insert(edge_idx);
                vertex_partitions[u].insert(part);
                vertex_partitions[v].insert(part);
            }
        }

        // Step 3: Classify vertices based on partition connectivity
        std::unordered_set<Index> separator_candidates;
        std::unordered_set<Index> part1_vertices, part2_vertices;
        std::unordered_set<Index> isolated_vertices;

        for (Index v = 0; v < n; ++v)
        {
            const auto &v_partitions = vertex_partitions[v];

            if (v_partitions.empty())
            {
                // Vertex has no edges (isolated)
                isolated_vertices.insert(v);
            }
            else if (v_partitions.size() == 1)
            {
                // Vertex connected to only one partition
                PartitionID single_part = *(v_partitions.begin());
                if (single_part == 0)
                {
                    part1_vertices.insert(v);
                }
                else
                {
                    part2_vertices.insert(v);
                }
            }
            else
            {
                // Vertex connected to multiple partitions -> separator candidate
                separator_candidates.insert(v);
            }
        }

        if (config_.verbose)
        {
            std::cout << "  Initial separator candidates: " << separator_candidates.size() << std::endl;
            std::cout << "  Initial part 1 vertices: " << part1_vertices.size() << std::endl;
            std::cout << "  Initial part 2 vertices: " << part2_vertices.size() << std::endl;
            std::cout << "  Isolated vertices: " << isolated_vertices.size() << std::endl;
        }

        // Step 4: Optimize separator by reducing its size
        // Try to move vertices from separator to parts if connectivity allows
        std::vector<Index> separator_list(separator_candidates.begin(), separator_candidates.end());

        for (Index sep_v : separator_list)
        {
            // Count connections to each part
            Index conn_to_part0 = 0, conn_to_part1 = 0;

            for (Index edge_idx : vertex_part_edges[sep_v])
            {
                PartitionID edge_part = partition[edge_idx];
                if (edge_part == 0)
                    conn_to_part0++;
                else if (edge_part == 1)
                    conn_to_part1++;
            }

            // Additional check: count direct neighbors in each part
            Index neighbors_in_part1 = 0, neighbors_in_part2 = 0;
            for (Index ptr = matrix.rowPtr()[sep_v]; ptr < matrix.rowPtr()[sep_v + 1]; ++ptr)
            {
                Index neighbor = matrix.colInd()[ptr];
                if (neighbor != sep_v)
                {
                    if (part1_vertices.count(neighbor))
                        neighbors_in_part1++;
                    if (part2_vertices.count(neighbor))
                        neighbors_in_part2++;
                }
            }

            // Decide if vertex can be moved to a specific part
            bool move_to_part1 = false, move_to_part2 = false;

            if (conn_to_part0 > 0 && conn_to_part1 == 0 && neighbors_in_part2 == 0)
            {
                move_to_part1 = true;
            }
            else if (conn_to_part1 > 0 && conn_to_part0 == 0 && neighbors_in_part1 == 0)
            {
                move_to_part2 = true;
            }
            else if (neighbors_in_part1 > 0 && neighbors_in_part2 == 0)
            {
                move_to_part1 = true;
            }
            else if (neighbors_in_part2 > 0 && neighbors_in_part1 == 0)
            {
                move_to_part2 = true;
            }

            // Move vertex if possible
            if (move_to_part1)
            {
                separator_candidates.erase(sep_v);
                part1_vertices.insert(sep_v);
            }
            else if (move_to_part2)
            {
                separator_candidates.erase(sep_v);
                part2_vertices.insert(sep_v);
            }
        }

        // Step 5: Handle isolated vertices
        // Assign them to the smaller part for balance
        for (Index iso_v : isolated_vertices)
        {
            if (part1_vertices.size() <= part2_vertices.size())
            {
                part1_vertices.insert(iso_v);
            }
            else
            {
                part2_vertices.insert(iso_v);
            }
        }

        // Step 6: Final validation and cleanup
        // Ensure no vertex appears in multiple sets
        for (Index v : separator_candidates)
        {
            part1_vertices.erase(v);
            part2_vertices.erase(v);
        }

        // Step 7: Build result
        VertexSeparatorResult result;

        result.separator.assign(separator_candidates.begin(), separator_candidates.end());
        result.part1.assign(part1_vertices.begin(), part1_vertices.end());
        result.part2.assign(part2_vertices.begin(), part2_vertices.end());

        // Sort for consistency
        std::sort(result.separator.begin(), result.separator.end());
        std::sort(result.part1.begin(), result.part1.end());
        std::sort(result.part2.begin(), result.part2.end());

        // Step 8: Validate result
        if (result.separator.size() + result.part1.size() + result.part2.size() != n)
        {
            if (config_.verbose)
            {
                std::cout << "Warning: Vertex count mismatch in separator result!" << std::endl;
                std::cout << "  Separator: " << result.separator.size() << std::endl;
                std::cout << "  Part1: " << result.part1.size() << std::endl;
                std::cout << "  Part2: " << result.part2.size() << std::endl;
                std::cout << "  Expected total: " << n << std::endl;
            }
        }

        // Check for degenerate partitions
        if (result.part1.empty() || result.part2.empty())
        {
            if (config_.verbose)
            {
                std::cout << "Warning: Degenerate partition detected!" << std::endl;
            }

            // Create a simple balanced partition as fallback
            result.separator.clear();
            result.part1.clear();
            result.part2.clear();

            for (Index v = 0; v < n; ++v)
            {
                if (v % 2 == 0)
                {
                    result.part1.push_back(v);
                }
                else
                {
                    result.part2.push_back(v);
                }
            }
        }

        stats_.separator_decode_time += timer.elapsed();

        if (config_.verbose)
        {
            std::cout << "C2 separator decoding completed:" << std::endl;
            std::cout << "  Final separator size: " << result.separator.size() << std::endl;
            std::cout << "  Final part 1 size: " << result.part1.size() << std::endl;
            std::cout << "  Final part 2 size: " << result.part2.size() << std::endl;
            std::cout << "  Separator ratio: " << std::fixed << std::setprecision(2)
                      << (100.0 * result.separator.size() / n) << "%" << std::endl;
            std::cout << "  Decoding time: " << timer.elapsed() << "s" << std::endl;
        }

        return result;
    }
    // VertexSeparatorResult HypergraphOrdering::decodeC2PartitionToVertexSeparator(
    //     const SparseMatrix &matrix,
    //     const std::vector<PartitionID> &partition,
    //     const std::vector<Index> &edge_rows,
    //     const std::vector<Index> &edge_cols) const
    // {

    //     const Index n = matrix.rows();
    //     const Index num_edges = edge_rows.size();

    //     if (partition.size() != num_edges)
    //     {
    //         throw std::runtime_error("Partition size mismatch with number of edges");
    //     }

    //     // Step 1: Build edge partition mapping
    //     std::unordered_set<Index> part0_edges, part1_edges;

    //     for (Index edge_idx = 0; edge_idx < num_edges; ++edge_idx)
    //     {
    //         if (partition[edge_idx] == 0)
    //         {
    //             part0_edges.insert(edge_idx);
    //         }
    //         else if (partition[edge_idx] == 1)
    //         {
    //             part1_edges.insert(edge_idx);
    //         }
    //     }

    //     // Step 2: Determine vertex partitions based on edge connectivity
    //     std::vector<std::unordered_set<Index>> vertex_part0_edges(n), vertex_part1_edges(n);

    //     // Collect edges incident to each vertex by partition
    //     for (Index edge_idx = 0; edge_idx < num_edges; ++edge_idx)
    //     {
    //         Index u = edge_rows[edge_idx];
    //         Index v = edge_cols[edge_idx];

    //         if (partition[edge_idx] == 0)
    //         {
    //             vertex_part0_edges[u].insert(edge_idx);
    //             vertex_part0_edges[v].insert(edge_idx);
    //         }
    //         else if (partition[edge_idx] == 1)
    //         {
    //             vertex_part1_edges[u].insert(edge_idx);
    //             vertex_part1_edges[v].insert(edge_idx);
    //         }
    //     }

    //     // Step 3: Classify vertices
    //     std::unordered_set<Index> separator_candidates;
    //     std::unordered_set<Index> part1_vertices, part2_vertices;

    //     for (Index v = 0; v < n; ++v)
    //     {
    //         bool has_part0_edges = !vertex_part0_edges[v].empty();
    //         bool has_part1_edges = !vertex_part1_edges[v].empty();

    //         if (has_part0_edges && has_part1_edges)
    //         {
    //             // Vertex is incident to edges in both partitions -> separator candidate
    //             separator_candidates.insert(v);
    //         }
    //         else if (has_part0_edges)
    //         {
    //             // Vertex only connected to partition 0 edges
    //             part1_vertices.insert(v);
    //         }
    //         else if (has_part1_edges)
    //         {
    //             // Vertex only connected to partition 1 edges
    //             part2_vertices.insert(v);
    //         }
    //         // Vertices with no edges are isolated - assign arbitrarily to part1
    //         else
    //         {
    //             part1_vertices.insert(v);
    //         }
    //     }

    //     // Step 4: Build result
    //     VertexSeparatorResult result;

    //     // Remove separator vertices from parts
    //     for (Index v : separator_candidates)
    //     {
    //         part1_vertices.erase(v);
    //         part2_vertices.erase(v);
    //     }

    //     // Convert to sorted vectors
    //     result.separator.assign(separator_candidates.begin(), separator_candidates.end());
    //     result.part1.assign(part1_vertices.begin(), part1_vertices.end());
    //     result.part2.assign(part2_vertices.begin(), part2_vertices.end());

    //     std::sort(result.separator.begin(), result.separator.end());
    //     std::sort(result.part1.begin(), result.part1.end());
    //     std::sort(result.part2.begin(), result.part2.end());

    //     if (config_.verbose)
    //     {
    //         std::cout << "C2 vertex separator decoded:" << std::endl;
    //         std::cout << "  Separator size: " << result.separator.size() << std::endl;
    //         std::cout << "  Part 1 size: " << result.part1.size() << std::endl;
    //         std::cout << "  Part 2 size: " << result.part2.size() << std::endl;
    //     }

    //     return result;
    // }

    VertexSeparatorResult HypergraphOrdering::decodeC3PartitionToVertexSeparator(
        const SparseMatrix &matrix,
        const std::vector<PartitionID> &partition) const
    {

        if (config_.verbose)
        {
            std::cout << "C3 partition decoding, falling back to C2" << std::endl;
        }

        // Fallback: reconstruct as C2 hypergraph and decode
        auto upper_tri = matrix.upperTriangular();
        std::vector<Index> edge_rows, edge_cols;

        for (Index i = 0; i < matrix.rows(); ++i)
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

        // Use C2 decoding as fallback
        return decodeC2PartitionToVertexSeparator(matrix, partition, edge_rows, edge_cols);
    }

}