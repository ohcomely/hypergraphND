#include "hypergraph_ordering.hpp"

namespace hypergraph_ordering
{
    HypergraphData HypergraphOrdering::constructCliqueNodeHypergraph(const SparseMatrix &matrix) const
    {
        Timer timer;

        if (config_.clique_type == OrderingConfig::C2_CLIQUES)
        {
            auto result = constructC2CliqueNodeHypergraph(matrix);
            stats_.hypergraph_construction_time += timer.elapsed();
            return result;
        }
        else
        {
            auto result = constructC3CliqueNodeHypergraph(matrix);
            stats_.hypergraph_construction_time += timer.elapsed();
            return result;
        }
    }

    HypergraphData HypergraphOrdering::constructC2CliqueNodeHypergraph(const SparseMatrix &matrix) const
    {
        const Index n = matrix.rows();

        if (config_.verbose)
        {
            std::cout << "Constructing C2 (edge-based) clique-node hypergraph..." << std::endl;
        }

        // Step 1: Extract edges from upper triangular part
        auto upper_tri = matrix.upperTriangular();
        std::vector<std::pair<Index, Index>> edges;

        for (Index i = 0; i < n; ++i)
        {
            for (Index ptr = upper_tri.rowPtr()[i]; ptr < upper_tri.rowPtr()[i + 1]; ++ptr)
            {
                Index j = upper_tri.colInd()[ptr];
                if (j > i)
                { // Strict upper triangular
                    edges.emplace_back(i, j);
                }
            }
        }

        const Index num_edges = edges.size();
        if (num_edges == 0)
        {
            // Empty hypergraph
            HypergraphData hg;
            hg.num_nodes = 0;
            hg.num_nets = 0;
            hg.num_pins = 0;
            return hg;
        }

        if (config_.verbose)
        {
            std::cout << "Found " << num_edges << " edges (hypergraph nodes)" << std::endl;
        }

        // Step 2: Build edge-to-node mapping
        std::unordered_map<std::pair<Index, Index>, Index, PairHash> edge_to_node;
        for (Index idx = 0; idx < num_edges; ++idx)
        {
            edge_to_node[edges[idx]] = idx;
        }

        // Step 3: Build nets for each vertex (collect edges incident to each vertex)
        std::vector<std::vector<Index>> vertex_to_nets(n);
        std::unordered_set<Index> vertices_with_edges;

        for (Index edge_idx = 0; edge_idx < num_edges; ++edge_idx)
        {
            Index u = edges[edge_idx].first;
            Index v = edges[edge_idx].second;

            vertex_to_nets[u].push_back(edge_idx);
            vertex_to_nets[v].push_back(edge_idx);
            vertices_with_edges.insert(u);
            vertices_with_edges.insert(v);
        }

        // Step 4: Create nets only for vertices that have edges
        std::vector<Index> vertex_list(vertices_with_edges.begin(), vertices_with_edges.end());
        std::sort(vertex_list.begin(), vertex_list.end());

        const Index num_nets = vertex_list.size();

        if (config_.verbose)
        {
            std::cout << "Creating " << num_nets << " nets (for vertices with edges)" << std::endl;
        }

        // Step 5: Build hypergraph in CSR format
        HypergraphData hg;
        hg.num_nodes = num_edges;
        hg.num_nets = num_nets;
        hg.num_pins = 0;

        // Count pins first
        for (Index vertex : vertex_list)
        {
            hg.num_pins += vertex_to_nets[vertex].size();
        }

        // Reserve memory
        hg.hyperedge_indices.reserve(num_nets + 1);
        hg.hyperedges.reserve(hg.num_pins);
        hg.node_weights.assign(num_edges, 1); // Unit weights for edges
        hg.net_weights.assign(num_nets, 1);   // Unit weights for nets

        // Build CSR representation
        hg.hyperedge_indices.push_back(0);

        for (Index vertex : vertex_list)
        {
            // Add all edges incident to this vertex
            for (Index edge_idx : vertex_to_nets[vertex])
            {
                hg.hyperedges.push_back(static_cast<kahypar_hyperedge_id_t>(edge_idx));
            }
            hg.hyperedge_indices.push_back(static_cast<HyperedgeIndex>(hg.hyperedges.size()));
        }

        if (config_.verbose)
        {
            std::cout << "C2 hypergraph constructed:" << std::endl;
            std::cout << "  Nodes (edges): " << hg.num_nodes << std::endl;
            std::cout << "  Nets (vertices): " << hg.num_nets << std::endl;
            std::cout << "  Pins: " << hg.num_pins << std::endl;
            std::cout << "  Density: " << (double)hg.num_pins / (hg.num_nodes * hg.num_nets) << std::endl;
        }

        return hg;
    }

    HypergraphData HypergraphOrdering::constructC3CliqueNodeHypergraph(const SparseMatrix &matrix) const
    {
        const Index n = matrix.rows();

        if (config_.verbose)
        {
            std::cout << "Constructing C3 (triangle-based) clique-node hypergraph..." << std::endl;
        }

        // Step 1: Find all triangles in the graph
        std::vector<std::array<Index, 3>> triangles;

        // Build adjacency lists for efficient triangle detection
        std::vector<std::vector<Index>> adj(n);
        for (Index i = 0; i < n; ++i)
        {
            for (Index ptr = matrix.rowPtr()[i]; ptr < matrix.rowPtr()[i + 1]; ++ptr)
            {
                Index j = matrix.colInd()[ptr];
                if (j != i)
                { // No self-loops
                    adj[i].push_back(j);
                }
            }
            // Sort for binary search
            std::sort(adj[i].begin(), adj[i].end());
        }

        // Find triangles using the standard algorithm
        for (Index u = 0; u < n; ++u)
        {
            for (Index v : adj[u])
            {
                if (v > u)
                { // Avoid duplicates
                    // Find common neighbors of u and v
                    for (Index w : adj[u])
                    {
                        if (w > v)
                        { // Maintain ordering u < v < w
                            // Check if (v,w) edge exists
                            if (std::binary_search(adj[v].begin(), adj[v].end(), w))
                            {
                                triangles.push_back({u, v, w});
                            }
                        }
                    }
                }
            }
        }

        const Index num_triangles = triangles.size();
        if (num_triangles == 0)
        {
            if (config_.verbose)
            {
                std::cout << "No triangles found, falling back to C2 hypergraph" << std::endl;
            }
            // Fallback to C2 if no triangles
            return constructC2CliqueNodeHypergraph(matrix);
        }

        if (config_.verbose)
        {
            std::cout << "Found " << num_triangles << " triangles (hypergraph nodes)" << std::endl;
        }

        // Step 2: Build triangle-to-node mapping
        std::unordered_map<std::array<Index, 3>, Index, TriangleHash> triangle_to_node;
        for (Index idx = 0; idx < num_triangles; ++idx)
        {
            triangle_to_node[triangles[idx]] = idx;
        }

        // Step 3: Build nets for each vertex (collect triangles incident to each vertex)
        std::vector<std::vector<Index>> vertex_to_nets(n);
        std::unordered_set<Index> vertices_with_triangles;

        for (Index tri_idx = 0; tri_idx < num_triangles; ++tri_idx)
        {
            const auto &triangle = triangles[tri_idx];
            for (Index i = 0; i < 3; ++i)
            {
                Index vertex = triangle[i];
                vertex_to_nets[vertex].push_back(tri_idx);
                vertices_with_triangles.insert(vertex);
            }
        }

        // Step 4: Create nets only for vertices that are in triangles
        std::vector<Index> vertex_list(vertices_with_triangles.begin(), vertices_with_triangles.end());
        std::sort(vertex_list.begin(), vertex_list.end());

        const Index num_nets = vertex_list.size();

        if (config_.verbose)
        {
            std::cout << "Creating " << num_nets << " nets (for vertices in triangles)" << std::endl;
        }

        // Step 5: Build hypergraph in CSR format
        HypergraphData hg;
        hg.num_nodes = num_triangles;
        hg.num_nets = num_nets;
        hg.num_pins = 0;

        // Count pins
        for (Index vertex : vertex_list)
        {
            hg.num_pins += vertex_to_nets[vertex].size();
        }

        // Reserve memory
        hg.hyperedge_indices.reserve(num_nets + 1);
        hg.hyperedges.reserve(hg.num_pins);
        hg.node_weights.assign(num_triangles, 1); // Unit weights for triangles
        hg.net_weights.assign(num_nets, 1);       // Unit weights for nets

        // Build CSR representation
        hg.hyperedge_indices.push_back(0);

        for (Index vertex : vertex_list)
        {
            // Add all triangles incident to this vertex
            for (Index tri_idx : vertex_to_nets[vertex])
            {
                hg.hyperedges.push_back(static_cast<kahypar_hyperedge_id_t>(tri_idx));
            }
            hg.hyperedge_indices.push_back(static_cast<HyperedgeIndex>(hg.hyperedges.size()));
        }

        if (config_.verbose)
        {
            std::cout << "C3 hypergraph constructed:" << std::endl;
            std::cout << "  Nodes (triangles): " << hg.num_nodes << std::endl;
            std::cout << "  Nets (vertices): " << hg.num_nets << std::endl;
            std::cout << "  Pins: " << hg.num_pins << std::endl;
            std::cout << "  Density: " << (double)hg.num_pins / (hg.num_nodes * hg.num_nets) << std::endl;
        }

        return hg;
    }

}