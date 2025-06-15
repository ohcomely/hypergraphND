#include "hypergraph_ordering.hpp"

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

        const Index n = matrix.rows();
        const Index num_edges = edge_rows.size();

        if (partition.size() != num_edges)
        {
            throw std::runtime_error("Partition size mismatch with number of edges");
        }

        // Step 1: Build edge partition mapping
        std::unordered_set<Index> part0_edges, part1_edges;

        for (Index edge_idx = 0; edge_idx < num_edges; ++edge_idx)
        {
            if (partition[edge_idx] == 0)
            {
                part0_edges.insert(edge_idx);
            }
            else if (partition[edge_idx] == 1)
            {
                part1_edges.insert(edge_idx);
            }
        }

        // Step 2: Determine vertex partitions based on edge connectivity
        std::vector<std::unordered_set<Index>> vertex_part0_edges(n), vertex_part1_edges(n);

        // Collect edges incident to each vertex by partition
        for (Index edge_idx = 0; edge_idx < num_edges; ++edge_idx)
        {
            Index u = edge_rows[edge_idx];
            Index v = edge_cols[edge_idx];

            if (partition[edge_idx] == 0)
            {
                vertex_part0_edges[u].insert(edge_idx);
                vertex_part0_edges[v].insert(edge_idx);
            }
            else if (partition[edge_idx] == 1)
            {
                vertex_part1_edges[u].insert(edge_idx);
                vertex_part1_edges[v].insert(edge_idx);
            }
        }

        // Step 3: Classify vertices
        std::unordered_set<Index> separator_candidates;
        std::unordered_set<Index> part1_vertices, part2_vertices;

        for (Index v = 0; v < n; ++v)
        {
            bool has_part0_edges = !vertex_part0_edges[v].empty();
            bool has_part1_edges = !vertex_part1_edges[v].empty();

            if (has_part0_edges && has_part1_edges)
            {
                // Vertex is incident to edges in both partitions -> separator candidate
                separator_candidates.insert(v);
            }
            else if (has_part0_edges)
            {
                // Vertex only connected to partition 0 edges
                part1_vertices.insert(v);
            }
            else if (has_part1_edges)
            {
                // Vertex only connected to partition 1 edges
                part2_vertices.insert(v);
            }
            // Vertices with no edges are isolated - assign arbitrarily to part1
            else
            {
                part1_vertices.insert(v);
            }
        }

        // Step 4: Build result
        VertexSeparatorResult result;

        // Remove separator vertices from parts
        for (Index v : separator_candidates)
        {
            part1_vertices.erase(v);
            part2_vertices.erase(v);
        }

        // Convert to sorted vectors
        result.separator.assign(separator_candidates.begin(), separator_candidates.end());
        result.part1.assign(part1_vertices.begin(), part1_vertices.end());
        result.part2.assign(part2_vertices.begin(), part2_vertices.end());

        std::sort(result.separator.begin(), result.separator.end());
        std::sort(result.part1.begin(), result.part1.end());
        std::sort(result.part2.begin(), result.part2.end());

        if (config_.verbose)
        {
            std::cout << "C2 vertex separator decoded:" << std::endl;
            std::cout << "  Separator size: " << result.separator.size() << std::endl;
            std::cout << "  Part 1 size: " << result.part1.size() << std::endl;
            std::cout << "  Part 2 size: " << result.part2.size() << std::endl;
        }

        return result;
    }

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