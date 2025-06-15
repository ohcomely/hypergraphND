#include "hypergraph_ordering.hpp"

namespace hypergraph_ordering
{
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

}