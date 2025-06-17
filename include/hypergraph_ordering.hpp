#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <array>

// Forward declarations
extern "C"
{
#include <libkahypar.h>
}

namespace hypergraph_ordering
{

    // Type aliases for clarity and potential optimization
    using Index = kahypar_hyperedge_id_t;
    using Weight = kahypar_hypernode_weight_t;
    using HyperedgeIndex = size_t;
    using PartitionID = int32_t;

    /**
     * Compressed Sparse Row (CSR) sparse matrix representation
     * Optimized for symmetric matrices in Matrix Market format
     */
    class SparseMatrix
    {
    public:
        // Constructors
        SparseMatrix() = default;
        SparseMatrix(Index n, Index nnz);
        SparseMatrix(const SparseMatrix &other) = default;
        SparseMatrix(SparseMatrix &&other) noexcept = default;

        // Assignment operators
        SparseMatrix &operator=(const SparseMatrix &other) = default;
        SparseMatrix &operator=(SparseMatrix &&other) noexcept = default;

        // Basic properties
        Index rows() const { return n_; }
        Index cols() const { return n_; }
        Index nnz() const { return nnz_; }
        bool empty() const { return nnz_ == 0; }

        // Data access
        const std::vector<Index> &rowPtr() const { return row_ptr_; }
        const std::vector<Index> &colInd() const { return col_ind_; }
        const std::vector<double> &values() const { return values_; }

        // Mutable access for construction
        std::vector<Index> &rowPtr() { return row_ptr_; }
        std::vector<Index> &colInd() { return col_ind_; }
        std::vector<double> &values() { return values_; }

        // Matrix operations
        static SparseMatrix loadFromMatrixMarket(const std::string &filename);
        void makeSymmetric();
        void eliminateZeros();
        SparseMatrix extractSubmatrix(const std::vector<Index> &vertices,
                                      std::unordered_map<Index, Index> &vertex_map) const;
        SparseMatrix upperTriangular() const;
        void buildFromTriplets(const std::vector<std::tuple<Index, Index, double>> &triplets);
        bool isSymmetric(double tolerance) const;
        void getStatistics(Index &min_degree, Index &max_degree,
                           double &avg_degree, double &std_degree) const;
        std::vector<std::vector<double>> toDense() const;

        // Utility functions
        std::pair<std::vector<Index>, std::vector<Index>> getEdges() const;
        void reserve(Index nnz_estimate);
        void finalize(); // Call after manual construction

        // Debug and validation
        bool isValid() const;
        void printInfo(std::ostream &os = std::cout) const;

    private:
        Index n_ = 0;                // Matrix dimension (n x n)
        Index nnz_ = 0;              // Number of non-zeros
        std::vector<Index> row_ptr_; // Row pointers (size n+1)
        std::vector<Index> col_ind_; // Column indices (size nnz)
        std::vector<double> values_; // Matrix values (size nnz)

        void sortIndices(); // Sort column indices within each row
    };

    /**
     * Hypergraph data structure for efficient partitioning
     * Uses CSR-like representation for hyperedges
     */
    struct HypergraphData
    {
        Index num_nodes = 0;
        Index num_nets = 0;
        Index num_pins = 0;

        // CSR representation of hypergraph
        std::vector<HyperedgeIndex> hyperedge_indices;  // Net boundaries (size num_nets + 1)
        std::vector<kahypar_hyperedge_id_t> hyperedges; // Pin list (size num_pins)

        // Weights
        std::vector<Weight> node_weights; // Node weights (size num_nodes)
        std::vector<Weight> net_weights;  // Net weights (size num_nets)

        // Validation and utility
        bool isValid() const;
        void printInfo(std::ostream &os = std::cout) const;

        // Memory management
        void reserve(Index nodes, Index nets, Index pins);
        void clear();
    };

    /**
     * Result of decoding hypergraph partition to vertex separator
     */
    struct VertexSeparatorResult
    {
        std::vector<Index> separator;
        std::vector<Index> part1;
        std::vector<Index> part2;

        void clear();
        bool isValid(Index total_vertices) const;
    };

    /**
     * Configuration for the ordering algorithm
     */
    struct OrderingConfig
    {
        std::string kahypar_config_path = "../cut_kKaHyPar_sea20.ini";
        Index max_recursion_depth = 10;
        Index min_subproblem_size = 200;
        Index min_nodes_for_partitioning = 100;
        bool use_minimum_degree = true;
        bool verbose = false;

        // Hypergraph construction options
        enum CliqueType
        {
            C2_CLIQUES = 2, // 2-cliques (edges)
            C3_CLIQUES = 3, // 3-cliques (triangles)
            C4_CLIQUES = 4  // 4-cliques (tetrahedra)
        };
        CliqueType clique_type = C4_CLIQUES;

        // KaHyPar specific settings
        double imbalance = 0.001;
        PartitionID num_blocks = 2;

        bool isValid() const;
    };

    /**
     * Main hypergraph-based matrix ordering class
     * Implements the Çatalyürek-Aykanat algorithm
     */
    class HypergraphOrdering
    {
    public:
        // Constructor
        explicit HypergraphOrdering(const OrderingConfig &config = OrderingConfig{});
        ~HypergraphOrdering();
        // Disable copy constructor and assignment (due to KaHyPar context)
        HypergraphOrdering(const HypergraphOrdering &) = delete;
        HypergraphOrdering &operator=(const HypergraphOrdering &) = delete;

        // Enable move constructor and assignment
        HypergraphOrdering(HypergraphOrdering &&other) noexcept;
        HypergraphOrdering &operator=(HypergraphOrdering &&other) noexcept;

        // Main interface
        std::vector<Index> orderMatrix(const SparseMatrix &matrix);
        std::vector<Index> orderMatrix(const std::string &matrix_file);

        // Configuration
        void setConfig(const OrderingConfig &config) { config_ = config; }
        const OrderingConfig &getConfig() const { return config_; }

        // Statistics and diagnostics
        struct Statistics
        {
            double total_time = 0.0;
            double hypergraph_construction_time = 0.0;
            double partitioning_time = 0.0;
            double separator_decode_time = 0.0;
            Index max_recursion_depth_reached = 0;
            Index total_subproblems = 0;
            Index hypergraph_partitions = 0;
            Index minimum_degree_calls = 0;
        };

        const Statistics &getStatistics() const { return stats_; }
        void resetStatistics() { stats_ = Statistics{}; }

    private:
        // Core algorithm components
        HypergraphData constructCliqueNodeHypergraph(const SparseMatrix &matrix) const;
        HypergraphData constructC2CliqueNodeHypergraph(const SparseMatrix &matrix) const;
        HypergraphData constructC2CliqueNodeHypergraphParallel(
            const SparseMatrix &matrix) const;
        HypergraphData constructC3CliqueNodeHypergraph(const SparseMatrix &matrix) const;
        HypergraphData constructC4CliqueNodeHypergraph(const SparseMatrix &matrix) const;

        std::vector<PartitionID> partitionHypergraph(const HypergraphData &hg) const;
        VertexSeparatorResult decodePartitionToVertexSeparator(
            const SparseMatrix &matrix,
            const std::vector<PartitionID> &partition,
            const std::vector<Index> &edge_rows,
            const std::vector<Index> &edge_cols) const;

        VertexSeparatorResult decodeC4PartitionToVertexSeparator(
            const SparseMatrix &matrix,
            const std::vector<PartitionID> &partition,
            const std::vector<std::vector<Index>> &all_cliques) const;

        // Recursive nested dissection
        void recursiveNestedDissection(
            const SparseMatrix &matrix,
            const std::vector<Index> &vertices,
            Index depth,
            std::vector<Index> &ordering) const;

        bool shouldTerminate(const SparseMatrix &matrix, const std::vector<Index> &vertices,
                             Index depth) const;

        std::vector<std::pair<Index, Index>> extractEdgesFromSubmatrix(
            const SparseMatrix &matrix) const;

        // Fallback ordering methods
        std::vector<Index> minimumDegreeOrdering(
            const SparseMatrix &matrix,
            const std::vector<Index> &vertices) const;

        std::vector<Index> simpleMinimumDegree(
            const SparseMatrix &matrix,
            const std::vector<Index> &vertices) const;

        std::vector<Index> approximateMinimumDegree(
            const SparseMatrix &matrix,
            const std::vector<Index> &vertices) const;

        Index calculateApproximateDegree(
            Index vertex,
            const std::vector<std::unordered_set<Index>> &adj,
            const std::vector<bool> &eliminated,
            const std::vector<Index> &weight) const;

        std::vector<Index> exactMinimumDegree(
            const SparseMatrix &matrix,
            const std::vector<Index> &vertices) const;

        // Helper functions
        std::vector<Index> remapVertices(
            const std::vector<Index> &vertices,
            const std::unordered_map<Index, Index> &vertex_map) const;

        bool shouldUseHypergraphPartitioning(
            const SparseMatrix &matrix,
            const std::vector<Index> &vertices,
            Index depth) const;

        VertexSeparatorResult decodeC2PartitionToVertexSeparator(
            const SparseMatrix &matrix,
            const std::vector<PartitionID> &partition,
            const std::vector<Index> &edge_rows,
            const std::vector<Index> &edge_cols) const;

        VertexSeparatorResult decodeC3PartitionToVertexSeparator(
            const SparseMatrix &matrix,
            const std::vector<PartitionID> &partition) const;

        bool validateVertexSeparator(
            const SparseMatrix &matrix,
            const VertexSeparatorResult &result) const;

        // KaHyPar interface
        kahypar_context_t *createKaHyParContext() const;
        void destroyKaHyParContext(kahypar_context_t *context) const;

        // AMD interface
        // std::vector<Index> amdOrdering(
        //     const SparseMatrix &matrix,
        //     const std::vector<Index> &vertices) const;
        // std::vector<Index> amdOrderingFull(const SparseMatrix &matrix) const;
        // void convertToAMDFormat(
        //     const SparseMatrix &matrix,
        //     std::vector<SuiteSparse_long> &Ap,
        //     std::vector<SuiteSparse_long> &Ai) const;

        // KaHyPar interface (now reuses context)
        void initializeKaHyParContext();
        void cleanupKaHyParContext();
        bool isKaHyParContextValid() const;
        void reinitializeKaHyParContext();

        kahypar_context_t *kahypar_context_;

        // Configuration and state
        OrderingConfig config_;
        mutable Statistics stats_;
        mutable std::vector<std::vector<Index>> stored_cliques_;

        // Timing utilities
        class Timer
        {
        public:
            Timer() : start_(std::chrono::high_resolution_clock::now()) {}
            double elapsed() const
            {
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start_);
                return duration.count() / 1e6;
            }

        private:
            std::chrono::high_resolution_clock::time_point start_;
        };

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

        mutable Timer timer_;

        class EliminationGraph
        {
        private:
            struct Vertex
            {
                std::unordered_set<Index> neighbors;
                Index degree;
                bool eliminated;

                Vertex() : degree(0), eliminated(false) {}
            };

            std::vector<Vertex> vertices_;
            Index n_;

        public:
            explicit EliminationGraph(Index n) : n_(n)
            {
                vertices_.resize(n);
            }

            void addEdge(Index u, Index v)
            {
                if (u != v && !vertices_[u].eliminated && !vertices_[v].eliminated)
                {
                    vertices_[u].neighbors.insert(v);
                    vertices_[v].neighbors.insert(u);
                }
            }

            void removeEdge(Index u, Index v)
            {
                vertices_[u].neighbors.erase(v);
                vertices_[v].neighbors.erase(u);
            }

            void computeInitialDegrees()
            {
                for (Index i = 0; i < n_; ++i)
                {
                    vertices_[i].degree = vertices_[i].neighbors.size();
                }
            }

            Index getDegree(Index v) const
            {
                return vertices_[v].degree;
            }

            bool isEliminated(Index v) const
            {
                return vertices_[v].eliminated;
            }

            Index findMinimumDegreeVertex() const
            {
                Index min_degree = std::numeric_limits<Index>::max();
                Index min_vertex = 0;
                bool found = false;

                for (Index i = 0; i < n_; ++i)
                {
                    if (!vertices_[i].eliminated && vertices_[i].degree < min_degree)
                    {
                        min_degree = vertices_[i].degree;
                        min_vertex = i;
                        found = true;
                    }
                }

                if (!found)
                {
                    throw std::runtime_error("No non-eliminated vertex found");
                }

                return min_vertex;
            }

            void eliminateVertex(Index pivot)
            {
                if (vertices_[pivot].eliminated)
                {
                    throw std::runtime_error("Vertex already eliminated");
                }

                // Get neighbors before elimination
                std::vector<Index> neighbors(vertices_[pivot].neighbors.begin(),
                                             vertices_[pivot].neighbors.end());

                // Remove pivot from all neighbor lists
                for (Index neighbor : neighbors)
                {
                    vertices_[neighbor].neighbors.erase(pivot);
                }

                // Create fill-in edges (clique among neighbors)
                for (size_t i = 0; i < neighbors.size(); ++i)
                {
                    for (size_t j = i + 1; j < neighbors.size(); ++j)
                    {
                        Index u = neighbors[i];
                        Index v = neighbors[j];

                        if (!vertices_[u].eliminated && !vertices_[v].eliminated)
                        {
                            // Add fill-in edge if it doesn't exist
                            if (vertices_[u].neighbors.find(v) == vertices_[u].neighbors.end())
                            {
                                addEdge(u, v);
                            }
                        }
                    }
                }

                // Mark pivot as eliminated
                vertices_[pivot].eliminated = true;
                vertices_[pivot].neighbors.clear();
                vertices_[pivot].degree = 0;

                // Update degrees of affected vertices
                std::unordered_set<Index> affected;
                for (Index neighbor : neighbors)
                {
                    if (!vertices_[neighbor].eliminated)
                    {
                        affected.insert(neighbor);
                        // Add all neighbors of this neighbor to affected set
                        for (Index nn : vertices_[neighbor].neighbors)
                        {
                            if (!vertices_[nn].eliminated)
                            {
                                affected.insert(nn);
                            }
                        }
                    }
                }

                // Recompute degrees for affected vertices
                for (Index v : affected)
                {
                    if (!vertices_[v].eliminated)
                    {
                        Index new_degree = 0;
                        for (Index neighbor : vertices_[v].neighbors)
                        {
                            if (!vertices_[neighbor].eliminated)
                            {
                                new_degree++;
                            }
                        }
                        vertices_[v].degree = new_degree;
                    }
                }
            }

            // Debug method
            void printGraph() const
            {
                for (Index i = 0; i < n_; ++i)
                {
                    if (!vertices_[i].eliminated)
                    {
                        std::cout << "Vertex " << i << " (degree " << vertices_[i].degree << "): ";
                        for (Index neighbor : vertices_[i].neighbors)
                        {
                            if (!vertices_[neighbor].eliminated)
                            {
                                std::cout << neighbor << " ";
                            }
                        }
                        std::cout << std::endl;
                    }
                }
            }

            // Get fill-in statistics
            Index getTotalFillIn() const
            {
                Index original_edges = 0;
                Index current_edges = 0;

                for (Index i = 0; i < n_; ++i)
                {
                    if (!vertices_[i].eliminated)
                    {
                        current_edges += vertices_[i].neighbors.size();
                    }
                }

                return (current_edges - original_edges) / 2; // Each edge counted twice
            }
        };
    };

    // Utility functions
    namespace utils
    {
        // File I/O helpers
        bool fileExists(const std::string &filename);
        std::string getFileExtension(const std::string &filename);

        // Vector utilities
        template <typename T>
        std::vector<T> range(T start, T end)
        {
            std::vector<T> result;
            result.reserve(end - start);
            for (T i = start; i < end; ++i)
            {
                result.push_back(i);
            }
            return result;
        }

        template <typename T>
        void removeDuplicates(std::vector<T> &vec)
        {
            std::sort(vec.begin(), vec.end());
            vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
        }

        // Set operations
        std::vector<Index> setDifference(
            const std::vector<Index> &set1,
            const std::unordered_set<Index> &set2);

        // Validation helpers
        inline bool isValidOrdering(const std::vector<Index> &ordering, Index matrix_size)
        {
            if (ordering.size() != matrix_size)
                return false;

            std::vector<bool> seen(matrix_size, false);
            for (Index v : ordering)
            {
                if (v >= matrix_size || seen[v])
                    return false;
                seen[v] = true;
            }
            return true;
        }
        bool isValidPartition(const std::vector<PartitionID> &partition, PartitionID num_parts);
    }

} // namespace hypergraph_ordering