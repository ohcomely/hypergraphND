#include "hypergraph_ordering.hpp"
#include <numeric>

namespace hypergraph_ordering
{
    std::vector<Index> HypergraphOrdering::minimumDegreeOrdering(
        const SparseMatrix &matrix,
        const std::vector<Index> &vertices) const
    {
        stats_.minimum_degree_calls++;

        // if (config_.verbose)
        if (true)
        {
            std::cout << "Using AMD ordering for " << vertices.size() << " vertices" << std::endl;
        }

        // If vertices represent the full matrix, use optimized full ordering
        if (vertices.size() == matrix.rows())
        {
            bool is_full_matrix = true;
            for (size_t i = 0; i < vertices.size(); ++i)
            {
                if (vertices[i] != static_cast<Index>(i))
                {
                    is_full_matrix = false;
                    break;
                }
            }
            if (is_full_matrix)
            {
                return amdOrderingFull(matrix);
            }
        }

        // Use submatrix AMD ordering
        return amdOrdering(matrix, vertices);
    }

    std::vector<Index> HypergraphOrdering::amdOrderingFull(const SparseMatrix &matrix) const
    {
        const Index n = matrix.rows();
        if (n == 0)
            return {};

        // Convert to AMD format
        std::vector<SuiteSparse_long> Ap, Ai;
        convertToAMDFormat(matrix, Ap, Ai);

        // Prepare AMD arrays
        std::vector<SuiteSparse_long> P(n);
        std::vector<double> Control(AMD_CONTROL), Info(AMD_INFO);

        // Set AMD control parameters for better performance
        amd_defaults(Control.data());
        Control[AMD_DENSE] = 10.0;     // Treat rows/cols as dense if > 10*sqrt(n) entries
        Control[AMD_AGGRESSIVE] = 1.0; // Perform aggressive absorption

        // Call AMD
        SuiteSparse_long result = amd_l_order(
            static_cast<SuiteSparse_long>(n),
            Ap.data(),
            Ai.data(),
            P.data(),
            Control.data(),
            Info.data());

        if (result != AMD_OK && result != AMD_OK_BUT_JUMBLED)
        {
            if (config_.verbose)
            {
                std::cout << "AMD failed with status " << result << ", using simple ordering" << std::endl;
            }
            // Fallback to identity ordering
            std::vector<Index> ordering(n);
            std::iota(ordering.begin(), ordering.end(), 0);
            return ordering;
        }

        // Convert result back to Index type
        std::vector<Index> ordering(n);
        for (Index i = 0; i < n; ++i)
        {
            ordering[i] = static_cast<Index>(P[i]);
        }

        if (config_.verbose)
        {
            std::cout << "AMD completed successfully. Fill-in: "
                      << Info[AMD_LNZ] << ", Flops: " << Info[AMD_NDIV] << std::endl;
        }

        return ordering;
    }

    std::vector<Index> HypergraphOrdering::amdOrdering(
        const SparseMatrix &matrix,
        const std::vector<Index> &vertices) const
    {
        const Index n = vertices.size();
        if (n == 0)
            return {};
        if (n == 1)
            return {vertices[0]};

        // Convert submatrix to AMD format
        std::vector<SuiteSparse_long> Ap, Ai;
        std::unordered_map<Index, SuiteSparse_long> vertex_map;
        convertSubmatrixToAMDFormat(matrix, vertices, Ap, Ai, vertex_map);

        // Prepare AMD arrays
        std::vector<SuiteSparse_long> P(n);
        std::vector<double> Control(AMD_CONTROL), Info(AMD_INFO);

        // Set AMD control parameters
        amd_defaults(Control.data());
        Control[AMD_DENSE] = 10.0;
        Control[AMD_AGGRESSIVE] = 1.0;

        // Call AMD
        SuiteSparse_long result = amd_l_order(
            static_cast<SuiteSparse_long>(n),
            Ap.data(),
            Ai.data(),
            P.data(),
            Control.data(),
            Info.data());

        if (result != AMD_OK && result != AMD_OK_BUT_JUMBLED)
        {
            if (config_.verbose)
            {
                std::cout << "AMD failed for submatrix, using simple ordering" << std::endl;
            }
            return vertices; // Return original ordering as fallback
        }

        // Map back to original vertex indices
        std::vector<Index> ordering(n);
        for (Index i = 0; i < n; ++i)
        {
            ordering[i] = vertices[P[i]];
        }

        return ordering;
    }

    void HypergraphOrdering::convertToAMDFormat(
        const SparseMatrix &matrix,
        std::vector<SuiteSparse_long> &Ap,
        std::vector<SuiteSparse_long> &Ai) const
    {
        const Index n = matrix.rows();
        const auto &row_ptr = matrix.rowPtr();
        const auto &col_ind = matrix.colInd();

        // AMD expects symmetric matrix in compressed column format
        // Since our matrix is symmetric, we can work with upper triangular part

        // First pass: count entries in each column
        std::vector<SuiteSparse_long> col_counts(n, 0);

        for (Index row = 0; row < n; ++row)
        {
            for (Index idx = row_ptr[row]; idx < row_ptr[row + 1]; ++idx)
            {
                Index col = col_ind[idx];
                if (col >= row)
                { // Upper triangular part
                    col_counts[col]++;
                    if (col != row)
                    { // Off-diagonal entries appear in both row and column
                        col_counts[row]++;
                    }
                }
            }
        }

        // Build column pointers
        Ap.resize(n + 1);
        Ap[0] = 0;
        for (Index col = 0; col < n; ++col)
        {
            Ap[col + 1] = Ap[col] + col_counts[col];
        }

        // Allocate row indices array
        Ai.resize(Ap[n]);

        // Second pass: fill in the row indices
        std::vector<SuiteSparse_long> col_pos = Ap; // Copy starting positions

        for (Index row = 0; row < n; ++row)
        {
            for (Index idx = row_ptr[row]; idx < row_ptr[row + 1]; ++idx)
            {
                Index col = col_ind[idx];
                if (col >= row)
                {
                    // Add entry (row, col)
                    Ai[col_pos[col]++] = static_cast<SuiteSparse_long>(row);
                    if (col != row)
                    {
                        // Add symmetric entry (col, row)
                        Ai[col_pos[row]++] = static_cast<SuiteSparse_long>(col);
                    }
                }
            }
        }

        // Sort row indices within each column (AMD prefers this)
        for (Index col = 0; col < n; ++col)
        {
            std::sort(Ai.begin() + Ap[col], Ai.begin() + Ap[col + 1]);
        }
    }

    void HypergraphOrdering::convertSubmatrixToAMDFormat(
        const SparseMatrix &matrix,
        const std::vector<Index> &vertices,
        std::vector<SuiteSparse_long> &Ap,
        std::vector<SuiteSparse_long> &Ai,
        std::unordered_map<Index, SuiteSparse_long> &vertex_map) const
    {
        const Index n = vertices.size();

        // Create mapping from original to local indices
        vertex_map.clear();
        for (Index i = 0; i < n; ++i)
        {
            vertex_map[vertices[i]] = static_cast<SuiteSparse_long>(i);
        }

        // Create vertex set for fast lookup
        std::unordered_set<Index> vertex_set(vertices.begin(), vertices.end());

        // First pass: count entries in each column of submatrix
        std::vector<SuiteSparse_long> col_counts(n, 0);

        const auto &row_ptr = matrix.rowPtr();
        const auto &col_ind = matrix.colInd();

        for (Index orig_row : vertices)
        {
            SuiteSparse_long local_row = vertex_map[orig_row];

            for (Index idx = row_ptr[orig_row]; idx < row_ptr[orig_row + 1]; ++idx)
            {
                Index orig_col = col_ind[idx];

                // Only include edges within the vertex subset
                if (vertex_set.find(orig_col) != vertex_set.end())
                {
                    SuiteSparse_long local_col = vertex_map[orig_col];

                    if (local_col >= local_row)
                    { // Upper triangular
                        col_counts[local_col]++;
                        if (local_col != local_row)
                        {
                            col_counts[local_row]++;
                        }
                    }
                }
            }
        }

        // Build column pointers
        Ap.resize(n + 1);
        Ap[0] = 0;
        for (Index col = 0; col < n; ++col)
        {
            Ap[col + 1] = Ap[col] + col_counts[col];
        }

        // Allocate and fill row indices
        Ai.resize(Ap[n]);
        std::vector<SuiteSparse_long> col_pos = Ap;

        for (Index orig_row : vertices)
        {
            SuiteSparse_long local_row = vertex_map[orig_row];

            for (Index idx = row_ptr[orig_row]; idx < row_ptr[orig_row + 1]; ++idx)
            {
                Index orig_col = col_ind[idx];

                if (vertex_set.find(orig_col) != vertex_set.end())
                {
                    SuiteSparse_long local_col = vertex_map[orig_col];

                    if (local_col >= local_row)
                    {
                        Ai[col_pos[local_col]++] = local_row;
                        if (local_col != local_row)
                        {
                            Ai[col_pos[local_row]++] = local_col;
                        }
                    }
                }
            }
        }

        // Sort row indices within each column
        for (Index col = 0; col < n; ++col)
        {
            std::sort(Ai.begin() + Ap[col], Ai.begin() + Ap[col + 1]);
        }
    }
}