#include "hypergraph_ordering.hpp"
#include <unordered_map>
#include <omp.h>
#include <numeric>
#include <iomanip>

namespace hypergraph_ordering
{
    HypergraphData HypergraphOrdering::constructCliqueNodeHypergraph(const SparseMatrix &matrix) const
    {
        Timer timer;

        if (config_.clique_type == OrderingConfig::C2_CLIQUES)
        {
            auto result = constructC2CliqueNodeHypergraph(matrix);
            // auto result = constructC2CliqueNodeHypergraphParallel(matrix);
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

    HypergraphData HypergraphOrdering::constructC2CliqueNodeHypergraphParallel(
        const SparseMatrix &matrix) const
    {

        const Index n = matrix.rows();

        // Step 1: Count edges per row in parallel
        std::vector<Index> row_edge_counts(n, 0);

#pragma omp parallel for
        for (Index i = 0; i < n; ++i)
        {
            for (Index ptr = matrix.rowPtr()[i]; ptr < matrix.rowPtr()[i + 1]; ++ptr)
            {
                Index j = matrix.colInd()[ptr];
                if (j > i)
                {
                    row_edge_counts[i]++;
                }
            }
        }

        // Step 2: Prefix sum for edge indexing
        std::vector<Index> edge_offsets(n + 1, 0);
        std::partial_sum(row_edge_counts.begin(), row_edge_counts.end(),
                         edge_offsets.begin() + 1);

        const Index total_edges = edge_offsets[n];

        // Step 3: Parallel edge collection with pre-allocated memory
        std::vector<std::pair<Index, Index>> edges(total_edges);

#pragma omp parallel for
        for (Index i = 0; i < n; ++i)
        {
            Index edge_idx = edge_offsets[i];
            for (Index ptr = matrix.rowPtr()[i]; ptr < matrix.rowPtr()[i + 1]; ++ptr)
            {
                Index j = matrix.colInd()[ptr];
                if (j > i)
                {
                    edges[edge_idx++] = {i, j};
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

    HypergraphData HypergraphOrdering::constructC4CliqueNodeHypergraph(const SparseMatrix &matrix) const
    {
        const Index n = matrix.rows();

        if (config_.verbose)
        {
            std::cout << "Constructing C4 (4-clique-based) clique-node hypergraph..." << std::endl;
        }

        // Build adjacency lists for efficient clique detection
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

        // Constants
        const Index INVALID_INDEX = static_cast<Index>(-1);

        // Data structures for C4 construction
        std::vector<std::array<Index, 4>> cliques_4;
        std::vector<std::array<Index, 3>> cliques_3;
        std::vector<std::pair<Index, Index>> cliques_2;

        // Track covered edges to avoid overlaps
        std::unordered_set<std::pair<Index, Index>, PairHash> covered_edges;

        // Parent pointers for efficient clique detection
        std::vector<Index> parent1(n, INVALID_INDEX); // For 3-cliques
        std::vector<Index> parent2(n, INVALID_INDEX); // For 4-cliques (second parent)

        // Main C4 construction algorithm
        // Process vertices in order to avoid duplicates
        for (Index vi = 0; vi < n; ++vi)
        {
            // Reset parent pointers for this iteration
            std::fill(parent1.begin(), parent1.end(), INVALID_INDEX);
            std::fill(parent2.begin(), parent2.end(), INVALID_INDEX);

            // Set first-level parents for neighbors of vi
            for (Index vj : adj[vi])
            {
                if (vj > vi) // Process each edge only once
                {
                    parent1[vj] = vi;
                }
            }

            // Look for 4-cliques starting from vi
            for (Index vj : adj[vi])
            {
                if (vj > vi)
                {
                    // Set second-level parents for neighbors of vj
                    for (Index vk : adj[vj])
                    {
                        if (vk > vj && parent1[vk] == vi) // vk is neighbor of both vi and vj
                        {
                            parent2[vk] = vj;
                        }
                    }

                    // Search for 4-cliques {vi, vj, vk, vl}
                    for (Index vk : adj[vj])
                    {
                        if (vk > vj && parent1[vk] == vi) // Triangle {vi, vj, vk} exists
                        {
                            for (Index vl : adj[vk])
                            {
                                if (vl > vk && parent2[vl] == vj && parent1[vl] == vi)
                                {
                                    // We have a 4-clique {vi, vj, vk, vl}
                                    // Check if we should include it based on coverage
                                    std::array<Index, 4> clique = {vi, vj, vk, vl};

                                    // Count how many edges are already covered
                                    Index covered_count = 0;
                                    std::vector<std::pair<Index, Index>> clique_edges = {
                                        {vi, vj}, {vi, vk}, {vi, vl}, {vj, vk}, {vj, vl}, {vk, vl}};

                                    for (const auto &edge : clique_edges)
                                    {
                                        if (covered_edges.count(edge))
                                        {
                                            covered_count++;
                                        }
                                    }

                                    // Include 4-clique if at most 3 edges are covered
                                    // This gives better pin count than using smaller cliques
                                    if (covered_count <= 3)
                                    {
                                        cliques_4.push_back(clique);

                                        // Mark all edges as covered
                                        for (const auto &edge : clique_edges)
                                        {
                                            covered_edges.insert(edge);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Find remaining 3-cliques not covered by 4-cliques
        for (Index vi = 0; vi < n; ++vi)
        {
            std::fill(parent1.begin(), parent1.end(), INVALID_INDEX);

            for (Index vj : adj[vi])
            {
                if (vj > vi)
                {
                    parent1[vj] = vi;
                }
            }

            for (Index vj : adj[vi])
            {
                if (vj > vi)
                {
                    for (Index vk : adj[vj])
                    {
                        if (vk > vj && parent1[vk] == vi)
                        {
                            // Triangle {vi, vj, vk} exists
                            std::vector<std::pair<Index, Index>> triangle_edges = {
                                {vi, vj}, {vi, vk}, {vj, vk}};

                            Index covered_count = 0;
                            for (const auto &edge : triangle_edges)
                            {
                                if (covered_edges.count(edge))
                                {
                                    covered_count++;
                                }
                            }

                            // Include 3-clique if at most 1 edge is covered
                            if (covered_count <= 1)
                            {
                                cliques_3.push_back({vi, vj, vk});

                                for (const auto &edge : triangle_edges)
                                {
                                    covered_edges.insert(edge);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Find remaining 2-cliques (edges) not covered by larger cliques
        for (Index i = 0; i < n; ++i)
        {
            for (Index ptr = matrix.rowPtr()[i]; ptr < matrix.rowPtr()[i + 1]; ++ptr)
            {
                Index j = matrix.colInd()[ptr];
                if (j > i) // Process each edge only once
                {
                    std::pair<Index, Index> edge = {i, j};
                    if (!covered_edges.count(edge))
                    {
                        cliques_2.push_back(edge);
                        covered_edges.insert(edge);
                    }
                }
            }
        }

        // Build the clique list in order: 4-cliques, 3-cliques, 2-cliques
        std::vector<std::vector<Index>> all_cliques;

        // Add 4-cliques
        for (const auto &clique : cliques_4)
        {
            all_cliques.push_back({clique[0], clique[1], clique[2], clique[3]});
        }

        // Add 3-cliques
        for (const auto &clique : cliques_3)
        {
            all_cliques.push_back({clique[0], clique[1], clique[2]});
        }

        // Add 2-cliques
        for (const auto &clique : cliques_2)
        {
            all_cliques.push_back({clique.first, clique.second});
        }

        stored_cliques_ = all_cliques;

        const Index num_cliques = all_cliques.size();
        if (num_cliques == 0)
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
            std::cout << "Found " << cliques_4.size() << " 4-cliques, "
                      << cliques_3.size() << " 3-cliques, "
                      << cliques_2.size() << " 2-cliques" << std::endl;
            std::cout << "Total cliques (hypergraph nodes): " << num_cliques << std::endl;
        }

        // Build nets for each vertex (collect cliques incident to each vertex)
        std::vector<std::vector<Index>> vertex_to_nets(n);
        std::unordered_set<Index> vertices_with_cliques;

        for (Index clique_idx = 0; clique_idx < num_cliques; ++clique_idx)
        {
            for (Index vertex : all_cliques[clique_idx])
            {
                vertex_to_nets[vertex].push_back(clique_idx);
                vertices_with_cliques.insert(vertex);
            }
        }

        // Create nets only for vertices that are in cliques
        std::vector<Index> vertex_list(vertices_with_cliques.begin(), vertices_with_cliques.end());
        std::sort(vertex_list.begin(), vertex_list.end());

        const Index num_nets = vertex_list.size();

        if (config_.verbose)
        {
            std::cout << "Creating " << num_nets << " nets (for vertices in cliques)" << std::endl;
        }

        // Build hypergraph in CSR format
        HypergraphData hg;
        hg.num_nodes = num_cliques;
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
        hg.node_weights.assign(num_cliques, 1); // Unit weights for cliques
        hg.net_weights.assign(num_nets, 1);     // Unit weights for nets

        // Build CSR representation
        hg.hyperedge_indices.push_back(0);

        for (Index vertex : vertex_list)
        {
            // Add all cliques incident to this vertex
            for (Index clique_idx : vertex_to_nets[vertex])
            {
                hg.hyperedges.push_back(static_cast<kahypar_hyperedge_id_t>(clique_idx));
            }
            hg.hyperedge_indices.push_back(static_cast<HyperedgeIndex>(hg.hyperedges.size()));
        }

        if (config_.verbose)
        {
            std::cout << "C4 hypergraph constructed:" << std::endl;
            std::cout << "  Nodes (cliques): " << hg.num_nodes << std::endl;
            std::cout << "  Nets (vertices): " << hg.num_nets << std::endl;
            std::cout << "  Pins: " << hg.num_pins << std::endl;
            std::cout << "  Density: " << (double)hg.num_pins / (hg.num_nodes * hg.num_nets) << std::endl;

            // Show clique distribution
            double avg_clique_size = 0.0;
            for (const auto &clique : all_cliques)
            {
                avg_clique_size += clique.size();
            }
            avg_clique_size /= num_cliques;
            std::cout << "  Average clique size: " << std::fixed << std::setprecision(2)
                      << avg_clique_size << std::endl;
        }

        return hg;
    }
}