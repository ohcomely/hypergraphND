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
        std::string kahypar_config_path = "km1_kKaHyPar_sea20.ini";
        Index max_recursion_depth = 10;
        Index min_subproblem_size = 200;
        Index min_nodes_for_partitioning = 100;
        bool use_minimum_degree = true;
        bool verbose = false;

        // Hypergraph construction options
        enum CliqueType
        {
            C2_CLIQUES = 2, // 2-cliques (edges)
            C3_CLIQUES = 3  // 3-cliques (triangles)
        };
        CliqueType clique_type = C2_CLIQUES;

        // KaHyPar specific settings
        double imbalance = 0.03;
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
        HypergraphData constructC3CliqueNodeHypergraph(const SparseMatrix &matrix) const;

        std::vector<PartitionID> partitionHypergraph(const HypergraphData &hg) const;
        VertexSeparatorResult decodePartitionToVertexSeparator(
            const SparseMatrix &matrix,
            const std::vector<PartitionID> &partition,
            const std::vector<Index> &edge_rows,
            const std::vector<Index> &edge_cols) const;

        // Recursive nested dissection
        void recursiveNestedDissection(
            const SparseMatrix &matrix,
            const std::vector<Index> &vertices,
            Index depth,
            std::vector<Index> &ordering) const;

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

        // Helper functions
        std::vector<Index> remapVertices(
            const std::vector<Index> &vertices,
            const std::unordered_map<Index, Index> &vertex_map) const;

        bool shouldUseHypergraphPartitioning(
            const SparseMatrix &matrix,
            const std::vector<Index> &vertices,
            Index depth) const;

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

        // Configuration and state
        OrderingConfig config_;
        mutable Statistics stats_;

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

        mutable Timer timer_;
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
        bool isValidOrdering(const std::vector<Index> &ordering, Index matrix_size);
        bool isValidPartition(const std::vector<PartitionID> &partition, PartitionID num_parts);
    }

} // namespace hypergraph_ordering