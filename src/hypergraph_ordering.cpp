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
    HypergraphOrdering::HypergraphOrdering(const OrderingConfig &config)
        : config_(config), stats_(), timer_()
    {
        // Validate the configuration
        if (!config.isValid())
        {
            throw std::invalid_argument("Invalid ordering configuration");
        }

        if (config_.verbose)
        {
            std::cout << "HypergraphOrdering initialized successfully." << std::endl;
            std::cout << "Configuration:" << std::endl;
            std::cout << "  Max recursion depth: " << config_.max_recursion_depth << std::endl;
            std::cout << "  Min subproblem size: " << config_.min_subproblem_size << std::endl;
            std::cout << "  Min nodes for partitioning: " << config_.min_nodes_for_partitioning << std::endl;
            std::cout << "  Imbalance tolerance: " << config_.imbalance << std::endl;
            std::cout << "  KaHyPar config: " << config_.kahypar_config_path << std::endl;
        }
    }

    kahypar_context_t *HypergraphOrdering::createKaHyParContext() const
    {
        kahypar_context_t *context = kahypar_context_new();
        kahypar_configure_context_from_file(context, config_.kahypar_config_path.c_str());

        // TODO: Can set KaHyPar seed as random
        kahypar_set_seed(context, 42);

        kahypar_supress_output(context, true);

        return context;
    }

    void HypergraphOrdering::destroyKaHyParContext(kahypar_context_t *context) const
    {
        if (context)
        {
            kahypar_context_free(context);

            if (config_.verbose)
            {
                std::cout << "KaHyPar context destroyed." << std::endl;
            }
        }
    }

    std::vector<PartitionID> HypergraphOrdering::partitionHypergraph(const HypergraphData &hg) const
    {
        Timer timer;

        // Create KaHyPar context
        kahypar_context_t *context = createKaHyParContext();

        try
        {
            if (!kahypar_validate_input(hg.num_nodes, hg.num_nets,
                                        hg.hyperedge_indices.data(), hg.hyperedges.data(),
                                        hg.net_weights.data(), hg.node_weights.data(), false))
            {
                throw std::runtime_error("Invalid hypergraph input for KaHyPar");
            }

            std::vector<PartitionID> partition(hg.num_nodes);
            kahypar_hyperedge_weight_t objective = 0;

            kahypar_partition(
                hg.num_nodes,                // num_vertices,
                hg.num_nets,                 // num_hyperedges
                config_.imbalance,           // epsilon
                config_.num_blocks,          // num_blocks
                hg.node_weights.data(),      // vertex_weights
                hg.net_weights.data(),       // hyperedge_weights
                hg.hyperedge_indices.data(), // hyperedge_indices
                hg.hyperedges.data(),        // hyperedges
                &objective,                  // objective
                context,                     // kahypar_context
                partition.data()             // partition
            );

            stats_.partitioning_time += timer.elapsed();
            stats_.hypergraph_partitions++;

            if (config_.verbose)
            {
                std::cout << "KaHyPar partitioning completed:" << std::endl;
                std::cout << "  Objective (cut): " << objective << std::endl;
                std::cout << "  Time: " << timer.elapsed() << " seconds" << std::endl;
            }

            destroyKaHyParContext(context);

            return partition;
        }
        catch (const std::exception &e)
        {
            destroyKaHyParContext(context);
            throw std::runtime_error("Error during hypergraph partitioning: " + std::string(e.what()));
        }
    }

    // Helper hash functions for unordered_map
    struct PairHash
    {
        std::size_t operator()(const std::pair<Index, Index> &p) const
        {
            return std::hash<Index>()(p.first) ^ (std::hash<Index>()(p.second) << 1);
        }
    };

    struct TriangleHash
    {
        std::size_t operator()(const std::array<Index, 3> &t) const
        {
            return std::hash<Index>()(t[0]) ^
                   (std::hash<Index>()(t[1]) << 1) ^
                   (std::hash<Index>()(t[2]) << 2);
        }
    };

    std::vector<Index> HypergraphOrdering::minimumDegreeOrdering(
        const SparseMatrix &matrix,
        const std::vector<Index> &vertices) const
    {

        Timer timer;

        if (vertices.empty())
        {
            return {};
        }

        // For small problems, use simple minimum degree
        if (vertices.size() <= 50)
        {
            auto result = simpleMinimumDegree(matrix, vertices);
            stats_.minimum_degree_calls++;
            return result;
        }

        // For larger problems, use approximate minimum degree
        auto result = approximateMinimumDegree(matrix, vertices);
        stats_.minimum_degree_calls++;

        if (config_.verbose)
        {
            std::cout << "AMD ordering completed for " << vertices.size()
                      << " vertices in " << timer.elapsed() << "s" << std::endl;
        }

        return result;
    }

    std::vector<Index> HypergraphOrdering::simpleMinimumDegree(
        const SparseMatrix &matrix,
        const std::vector<Index> &vertices) const
    {

        const Index n = vertices.size();
        if (n <= 1)
            return vertices;

        // Create vertex mapping for subproblem
        std::unordered_map<Index, Index> global_to_local;
        for (Index i = 0; i < n; ++i)
        {
            global_to_local[vertices[i]] = i;
        }

        // Build adjacency lists for submatrix
        std::vector<std::unordered_set<Index>> adj(n);

        for (Index i = 0; i < n; ++i)
        {
            Index global_i = vertices[i];

            // Get neighbors from matrix
            for (Index ptr = matrix.rowPtr()[global_i]; ptr < matrix.rowPtr()[global_i + 1]; ++ptr)
            {
                Index global_j = matrix.colInd()[ptr];

                // Only include if both vertices are in our subset
                auto it = global_to_local.find(global_j);
                if (it != global_to_local.end() && it->second != i)
                {
                    Index local_j = it->second;
                    adj[i].insert(local_j);
                    adj[local_j].insert(i);
                }
            }
        }

        // Simple minimum degree elimination
        std::vector<Index> ordering;
        std::vector<bool> eliminated(n, false);
        ordering.reserve(n);

        for (Index step = 0; step < n; ++step)
        {
            // Find vertex with minimum degree
            Index min_degree = std::numeric_limits<Index>::max();
            Index min_vertex = 0;

            for (Index i = 0; i < n; ++i)
            {
                if (!eliminated[i])
                {
                    Index degree = 0;
                    for (Index neighbor : adj[i])
                    {
                        if (!eliminated[neighbor])
                        {
                            degree++;
                        }
                    }

                    if (degree < min_degree)
                    {
                        min_degree = degree;
                        min_vertex = i;
                    }
                }
            }

            // Eliminate the minimum degree vertex
            eliminated[min_vertex] = true;
            ordering.push_back(vertices[min_vertex]);

            // Update adjacency (make neighbors of min_vertex adjacent to each other)
            std::vector<Index> neighbors;
            for (Index neighbor : adj[min_vertex])
            {
                if (!eliminated[neighbor])
                {
                    neighbors.push_back(neighbor);
                }
            }

            // Create clique among neighbors (fill-in)
            for (Index i = 0; i < neighbors.size(); ++i)
            {
                for (Index j = i + 1; j < neighbors.size(); ++j)
                {
                    adj[neighbors[i]].insert(neighbors[j]);
                    adj[neighbors[j]].insert(neighbors[i]);
                }
            }
        }

        return ordering;
    }

    std::vector<Index> HypergraphOrdering::approximateMinimumDegree(
        const SparseMatrix &matrix,
        const std::vector<Index> &vertices) const
    {

        const Index n = vertices.size();
        if (n <= 1)
            return vertices;

        // Create vertex mapping
        std::unordered_map<Index, Index> global_to_local;
        for (Index i = 0; i < n; ++i)
        {
            global_to_local[vertices[i]] = i;
        }

        // Build initial adjacency structure with approximate degrees
        std::vector<std::unordered_set<Index>> adj(n);
        std::vector<Index> degree(n, 0);

        // Build adjacency lists
        for (Index i = 0; i < n; ++i)
        {
            Index global_i = vertices[i];

            for (Index ptr = matrix.rowPtr()[global_i]; ptr < matrix.rowPtr()[global_i + 1]; ++ptr)
            {
                Index global_j = matrix.colInd()[ptr];

                auto it = global_to_local.find(global_j);
                if (it != global_to_local.end() && it->second != i)
                {
                    Index local_j = it->second;
                    adj[i].insert(local_j);
                    adj[local_j].insert(i);
                }
            }
            degree[i] = adj[i].size();
        }

        // AMD elimination with approximations
        std::vector<Index> ordering;
        std::vector<bool> eliminated(n, false);
        std::vector<Index> weight(n, 1); // Supernode weights
        ordering.reserve(n);

        // Use a simple priority queue (could be optimized with heap)
        auto findMinDegreeVertex = [&]() -> Index
        {
            Index min_score = std::numeric_limits<Index>::max();
            Index min_vertex = 0;
            bool found = false;

            for (Index i = 0; i < n; ++i)
            {
                if (!eliminated[i])
                {
                    // Approximate degree calculation
                    Index approx_degree = 0;
                    for (Index neighbor : adj[i])
                    {
                        if (!eliminated[neighbor])
                        {
                            approx_degree += weight[neighbor];
                        }
                    }

                    // Simple scoring: degree + small random tie-breaker
                    Index score = approx_degree * 100 + (i % 13); // Simple hash for tie-breaking

                    if (!found || score < min_score)
                    {
                        min_score = score;
                        min_vertex = i;
                        found = true;
                    }
                }
            }
            return min_vertex;
        };

        // Elimination loop
        for (Index step = 0; step < n; ++step)
        {
            Index pivot = findMinDegreeVertex();

            // Eliminate pivot
            eliminated[pivot] = true;
            ordering.push_back(vertices[pivot]);

            // Collect non-eliminated neighbors
            std::vector<Index> neighbors;
            for (Index neighbor : adj[pivot])
            {
                if (!eliminated[neighbor])
                {
                    neighbors.push_back(neighbor);
                }
            }

            // Update graph structure (simplified element absorption)
            if (neighbors.size() > 1)
            {
                // Find a "principal" neighbor to absorb others
                Index principal = neighbors[0];
                Index max_weight = weight[principal];

                for (Index i = 1; i < neighbors.size(); ++i)
                {
                    if (weight[neighbors[i]] > max_weight)
                    {
                        max_weight = weight[neighbors[i]];
                        principal = neighbors[i];
                    }
                }

                // Absorb smaller neighbors into principal (approximate)
                for (Index neighbor : neighbors)
                {
                    if (neighbor != principal)
                    {
                        // Add neighbor's adjacencies to principal
                        for (Index adj_vertex : adj[neighbor])
                        {
                            if (!eliminated[adj_vertex] && adj_vertex != principal)
                            {
                                adj[principal].insert(adj_vertex);
                                adj[adj_vertex].insert(principal);
                                adj[adj_vertex].erase(neighbor);
                            }
                        }
                        // Increase weight of principal
                        weight[principal] += weight[neighbor];
                        adj[neighbor].clear();
                    }
                }

                // Create clique among remaining neighbors (limited to avoid excessive fill)
                if (neighbors.size() <= 10)
                { // Limit clique size for efficiency
                    for (Index i = 0; i < neighbors.size(); ++i)
                    {
                        if (eliminated[neighbors[i]])
                            continue;
                        for (Index j = i + 1; j < neighbors.size(); ++j)
                        {
                            if (eliminated[neighbors[j]])
                                continue;
                            adj[neighbors[i]].insert(neighbors[j]);
                            adj[neighbors[j]].insert(neighbors[i]);
                        }
                    }
                }
            }

            // Clean up pivot's adjacencies
            for (Index neighbor : adj[pivot])
            {
                adj[neighbor].erase(pivot);
            }
            adj[pivot].clear();
        }

        return ordering;
    }

    // Helper method for degree calculation
    Index HypergraphOrdering::calculateApproximateDegree(
        Index vertex,
        const std::vector<std::unordered_set<Index>> &adj,
        const std::vector<bool> &eliminated,
        const std::vector<Index> &weight) const
    {

        Index degree = 0;
        for (Index neighbor : adj[vertex])
        {
            if (!eliminated[neighbor])
            {
                degree += weight[neighbor];
            }
        }
        return degree;
    }

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

    // HypergraphData HypergraphOrdering::constructC2CliqueNodeHypergraph(const SparseMatrix &matrix) const
    // {
    //     const Index n = matrix.rows();

    //     // Single pass: count edges and build vertex sets directly
    //     std::vector<std::vector<Index>> vertex_to_edges(n);
    //     std::vector<Index> vertices_with_edges;
    //     std::vector<bool> has_edges(n, false);

    //     Index edge_count = 0;
    //     for (Index i = 0; i < n; ++i)
    //     {
    //         for (Index ptr = matrix.rowPtr()[i]; ptr < matrix.rowPtr()[i + 1]; ++ptr)
    //         {
    //             Index j = matrix.colInd()[ptr];
    //             if (j > i)
    //             { // Upper triangular only
    //                 vertex_to_edges[i].push_back(edge_count);
    //                 vertex_to_edges[j].push_back(edge_count);

    //                 if (!has_edges[i])
    //                 {
    //                     vertices_with_edges.push_back(i);
    //                     has_edges[i] = true;
    //                 }
    //                 if (!has_edges[j])
    //                 {
    //                     vertices_with_edges.push_back(j);
    //                     has_edges[j] = true;
    //                 }

    //                 ++edge_count;
    //             }
    //         }
    //     }

    //     // Build hypergraph directly (no intermediate storage)
    //     HypergraphData hg;
    //     hg.num_nodes = edge_count;
    //     hg.num_nets = vertices_with_edges.size();

    //     // Reserve exact memory
    //     Index total_pins = 0;
    //     for (Index v : vertices_with_edges)
    //     {
    //         total_pins += vertex_to_edges[v].size();
    //     }
    //     hg.num_pins = total_pins;

    //     // Build CSR directly
    //     hg.hyperedge_indices.reserve(hg.num_nets + 1);
    //     hg.hyperedges.reserve(total_pins);

    //     std::sort(vertices_with_edges.begin(), vertices_with_edges.end());
    //     hg.hyperedge_indices.push_back(0);

    //     for (Index vertex : vertices_with_edges)
    //     {
    //         hg.hyperedges.insert(hg.hyperedges.end(),
    //                              vertex_to_edges[vertex].begin(),
    //                              vertex_to_edges[vertex].end());
    //         hg.hyperedge_indices.push_back(hg.hyperedges.size());
    //     }

    //     if (config_.verbose)
    //     {
    //         std::cout << "C2 hypergraph constructed:" << std::endl;
    //         std::cout << "  Nodes (edges): " << hg.num_nodes << std::endl;
    //         std::cout << "  Nets (vertices): " << hg.num_nets << std::endl;
    //         std::cout << "  Pins: " << hg.num_pins << std::endl;
    //         std::cout << "  Density: " << (double)hg.num_pins / (hg.num_nodes * hg.num_nets) << std::endl;
    //     }

    //     return hg;
    // }

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

    // Add this to your hypergraph_ordering.cpp file

    // Add this complete block to hypergraph_ordering.cpp

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