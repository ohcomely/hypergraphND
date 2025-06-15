#include "hypergraph_ordering.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <iostream>
#include <iomanip>
#include <map>
#include <unordered_set>
#include <math.h>

namespace hypergraph_ordering
{

    // Constructors
    SparseMatrix::SparseMatrix(Index n, Index nnz)
        : n_(n), nnz_(nnz)
    {
        row_ptr_.resize(n + 1, 0);
        col_ind_.reserve(nnz);
        values_.reserve(nnz);
    }

    // Static method to load from Matrix Market format
    SparseMatrix SparseMatrix::loadFromMatrixMarket(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        std::string line;

        // Skip header comments
        while (std::getline(file, line))
        {
            if (line.empty() || line[0] != '%')
            {
                break;
            }
        }

        // Parse matrix dimensions
        std::istringstream iss(line);
        Index rows, cols, nnz;
        if (!(iss >> rows >> cols >> nnz))
        {
            throw std::runtime_error("Invalid Matrix Market header");
        }

        if (rows != cols)
        {
            throw std::runtime_error("Matrix must be square");
        }

        // Create matrix
        SparseMatrix matrix(rows, nnz);

        // Temporary storage for COO format
        std::vector<std::tuple<Index, Index, double>> triplets;
        triplets.reserve(nnz);

        // Read matrix entries
        Index entries_read = 0;
        while (std::getline(file, line) && entries_read < nnz)
        {
            if (line.empty() || line[0] == '%')
                continue;

            std::istringstream line_stream(line);
            Index i, j;
            double val = 1.0; // Default value for pattern matrices

            if (!(line_stream >> i >> j))
            {
                throw std::runtime_error("Invalid matrix entry format");
            }

            // Try to read value (optional for pattern matrices)
            line_stream >> val;

            // Convert to 0-based indexing
            i--;
            j--;

            if (i >= rows || j >= cols || i < 0 || j < 0)
            {
                throw std::runtime_error("Matrix indices out of bounds");
            }

            triplets.emplace_back(i, j, val);
            entries_read++;
        }

        // Convert COO to CSR
        matrix.buildFromTriplets(triplets);

        // Make symmetric if not already
        matrix.makeSymmetric();
        matrix.eliminateZeros();

        return matrix;
    }

    // Build CSR matrix from COO triplets
    void SparseMatrix::buildFromTriplets(const std::vector<std::tuple<Index, Index, double>> &triplets)
    {
        // Count entries per row
        std::vector<Index> row_count(n_, 0);
        for (const auto &triplet : triplets)
        {
            Index i = std::get<0>(triplet);
            row_count[i]++;
        }

        // Build row pointers
        row_ptr_[0] = 0;
        for (Index i = 0; i < n_; ++i)
        {
            row_ptr_[i + 1] = row_ptr_[i] + row_count[i];
        }

        // Resize storage
        col_ind_.resize(row_ptr_[n_]);
        values_.resize(row_ptr_[n_]);

        // Fill in entries
        std::vector<Index> current_pos = row_ptr_;
        for (const auto &triplet : triplets)
        {
            Index i = std::get<0>(triplet);
            Index j = std::get<1>(triplet);
            double val = std::get<2>(triplet);

            Index pos = current_pos[i]++;
            col_ind_[pos] = j;
            values_[pos] = val;
        }

        // Sort indices within each row
        sortIndices();

        // Update nnz
        nnz_ = col_ind_.size();
    }

    // Make matrix symmetric
    void SparseMatrix::makeSymmetric()
    {
        std::vector<std::tuple<Index, Index, double>> all_triplets;

        // Collect existing entries
        for (Index i = 0; i < n_; ++i)
        {
            for (Index ptr = row_ptr_[i]; ptr < row_ptr_[i + 1]; ++ptr)
            {
                Index j = col_ind_[ptr];
                double val = values_[ptr];
                all_triplets.emplace_back(i, j, val);

                // Add symmetric entry if not diagonal
                if (i != j)
                {
                    all_triplets.emplace_back(j, i, val);
                }
            }
        }

        // Remove duplicates and sum values
        std::map<std::pair<Index, Index>, double> entry_map;
        for (const auto &triplet : all_triplets)
        {
            Index i = std::get<0>(triplet);
            Index j = std::get<1>(triplet);
            double val = std::get<2>(triplet);

            entry_map[{i, j}] += val;
        }

        // Convert back to triplets
        std::vector<std::tuple<Index, Index, double>> sym_triplets;
        for (const auto &entry : entry_map)
        {
            if (std::abs(entry.second) > 1e-15)
            { // Skip near-zero entries
                sym_triplets.emplace_back(entry.first.first, entry.first.second, entry.second);
            }
        }

        // Rebuild matrix
        buildFromTriplets(sym_triplets);
    }

    // Eliminate zero entries
    void SparseMatrix::eliminateZeros()
    {
        std::vector<Index> new_col_ind;
        std::vector<double> new_values;
        std::vector<Index> new_row_ptr(n_ + 1, 0);

        for (Index i = 0; i < n_; ++i)
        {
            new_row_ptr[i] = new_col_ind.size();

            for (Index ptr = row_ptr_[i]; ptr < row_ptr_[i + 1]; ++ptr)
            {
                if (std::abs(values_[ptr]) > 1e-15)
                {
                    new_col_ind.push_back(col_ind_[ptr]);
                    new_values.push_back(values_[ptr]);
                }
            }
        }
        new_row_ptr[n_] = new_col_ind.size();

        // Update matrix
        col_ind_ = std::move(new_col_ind);
        values_ = std::move(new_values);
        row_ptr_ = std::move(new_row_ptr);
        nnz_ = col_ind_.size();
    }

    // Extract submatrix for given vertices
    SparseMatrix SparseMatrix::extractSubmatrix(const std::vector<Index> &vertices,
                                                std::unordered_map<Index, Index> &vertex_map) const
    {
        const Index sub_n = vertices.size();

        // Create vertex mapping
        vertex_map.clear();
        for (Index i = 0; i < sub_n; ++i)
        {
            vertex_map[vertices[i]] = i;
        }

        // Create vertex set for fast lookup
        std::unordered_set<Index> vertex_set(vertices.begin(), vertices.end());

        // Collect submatrix entries
        std::vector<std::tuple<Index, Index, double>> sub_triplets;

        for (Index i = 0; i < sub_n; ++i)
        {
            Index global_i = vertices[i];

            for (Index ptr = row_ptr_[global_i]; ptr < row_ptr_[global_i + 1]; ++ptr)
            {
                Index global_j = col_ind_[ptr];
                double val = values_[ptr];

                // Only include if both vertices are in the submatrix
                if (vertex_set.count(global_j))
                {
                    Index local_j = vertex_map[global_j];
                    sub_triplets.emplace_back(i, local_j, val);
                }
            }
        }

        // Build submatrix
        SparseMatrix submatrix(sub_n, sub_triplets.size());
        submatrix.buildFromTriplets(sub_triplets);

        return submatrix;
    }

    // Get upper triangular part
    SparseMatrix SparseMatrix::upperTriangular() const
    {
        std::vector<std::tuple<Index, Index, double>> upper_triplets;

        for (Index i = 0; i < n_; ++i)
        {
            for (Index ptr = row_ptr_[i]; ptr < row_ptr_[i + 1]; ++ptr)
            {
                Index j = col_ind_[ptr];
                double val = values_[ptr];

                if (j >= i)
                { // Upper triangular (including diagonal)
                    upper_triplets.emplace_back(i, j, val);
                }
            }
        }

        SparseMatrix upper(n_, upper_triplets.size());
        upper.buildFromTriplets(upper_triplets);

        return upper;
    }

    // Get edges as (row, col) pairs
    std::pair<std::vector<Index>, std::vector<Index>> SparseMatrix::getEdges() const
    {
        std::vector<Index> edge_rows, edge_cols;

        // Get upper triangular edges only
        for (Index i = 0; i < n_; ++i)
        {
            for (Index ptr = row_ptr_[i]; ptr < row_ptr_[i + 1]; ++ptr)
            {
                Index j = col_ind_[ptr];
                if (j > i)
                { // Strict upper triangular
                    edge_rows.push_back(i);
                    edge_cols.push_back(j);
                }
            }
        }

        return {edge_rows, edge_cols};
    }

    // Reserve memory
    void SparseMatrix::reserve(Index nnz_estimate)
    {
        col_ind_.reserve(nnz_estimate);
        values_.reserve(nnz_estimate);
    }

    // Finalize after manual construction
    void SparseMatrix::finalize()
    {
        sortIndices();
        eliminateZeros();
    }

    // Validation
    bool SparseMatrix::isValid() const
    {
        // Check dimensions
        if (n_ <= 0 || nnz_ < 0)
            return false;
        if (row_ptr_.size() != n_ + 1)
            return false;
        if (col_ind_.size() != nnz_ || values_.size() != nnz_)
            return false;

        // Check row pointers
        if (row_ptr_[0] != 0 || row_ptr_[n_] != nnz_)
            return false;

        for (Index i = 0; i < n_; ++i)
        {
            if (row_ptr_[i] > row_ptr_[i + 1])
                return false;

            // Check column indices
            for (Index ptr = row_ptr_[i]; ptr < row_ptr_[i + 1]; ++ptr)
            {
                if (col_ind_[ptr] >= n_)
                    return false;

                // Check if sorted (optional)
                if (ptr > row_ptr_[i] && col_ind_[ptr] <= col_ind_[ptr - 1])
                {
                    // Not sorted within row
                }
            }
        }

        return true;
    }

    // Print matrix information
    void SparseMatrix::printInfo(std::ostream &os) const
    {
        os << "SparseMatrix Info:" << std::endl;
        os << "  Dimensions: " << n_ << " x " << n_ << std::endl;
        os << "  Non-zeros: " << nnz_ << std::endl;
        os << "  Density: " << std::fixed << std::setprecision(6)
           << (double)nnz_ / (n_ * n_) << std::endl;

        if (n_ > 0)
        {
            os << "  Avg nnz per row: " << (double)nnz_ / n_ << std::endl;

            // Find min/max row nnz
            Index min_nnz = n_, max_nnz = 0;
            for (Index i = 0; i < n_; ++i)
            {
                Index row_nnz = row_ptr_[i + 1] - row_ptr_[i];
                min_nnz = std::min(min_nnz, row_nnz);
                max_nnz = std::max(max_nnz, row_nnz);
            }
            os << "  Min/Max row nnz: " << min_nnz << " / " << max_nnz << std::endl;
        }

        os << "  Valid: " << (isValid() ? "Yes" : "No") << std::endl;
    }

    // Sort column indices within each row
    void SparseMatrix::sortIndices()
    {
        for (Index i = 0; i < n_; ++i)
        {
            Index start = row_ptr_[i];
            Index end = row_ptr_[i + 1];

            if (end > start)
            {
                // Create index pairs for sorting
                std::vector<std::pair<Index, double>> pairs;
                pairs.reserve(end - start);

                for (Index ptr = start; ptr < end; ++ptr)
                {
                    pairs.emplace_back(col_ind_[ptr], values_[ptr]);
                }

                // Sort by column index
                std::sort(pairs.begin(), pairs.end());

                // Write back
                for (Index k = 0; k < pairs.size(); ++k)
                {
                    col_ind_[start + k] = pairs[k].first;
                    values_[start + k] = pairs[k].second;
                }
            }
        }
    }

    // Additional utility methods

    // Check if matrix is symmetric (approximately)
    bool SparseMatrix::isSymmetric(double tolerance) const
    {
        for (Index i = 0; i < n_; ++i)
        {
            for (Index ptr = row_ptr_[i]; ptr < row_ptr_[i + 1]; ++ptr)
            {
                Index j = col_ind_[ptr];
                double val = values_[ptr];

                if (i != j)
                {
                    // Find (j,i) entry
                    bool found = false;
                    for (Index ptr2 = row_ptr_[j]; ptr2 < row_ptr_[j + 1]; ++ptr2)
                    {
                        if (col_ind_[ptr2] == i)
                        {
                            if (std::abs(values_[ptr2] - val) > tolerance)
                            {
                                return false;
                            }
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                        return false;
                }
            }
        }
        return true;
    }

    // Get matrix statistics
    void SparseMatrix::getStatistics(Index &min_degree, Index &max_degree,
                                     double &avg_degree, double &std_degree) const
    {
        if (n_ == 0)
        {
            min_degree = max_degree = 0;
            avg_degree = std_degree = 0.0;
            return;
        }

        std::vector<Index> degrees(n_);
        min_degree = n_;
        max_degree = 0;
        Index total_degree = 0;

        for (Index i = 0; i < n_; ++i)
        {
            degrees[i] = row_ptr_[i + 1] - row_ptr_[i];
            min_degree = std::min(min_degree, degrees[i]);
            max_degree = std::max(max_degree, degrees[i]);
            total_degree += degrees[i];
        }

        avg_degree = (double)total_degree / n_;

        // Calculate standard deviation
        double variance = 0.0;
        for (Index i = 0; i < n_; ++i)
        {
            double diff = degrees[i] - avg_degree;
            variance += diff * diff;
        }
        std_degree = std::sqrt(variance / n_);
    }

    // Convert to dense matrix (for debugging small matrices)
    std::vector<std::vector<double>> SparseMatrix::toDense() const
    {
        std::vector<std::vector<double>> dense(n_, std::vector<double>(n_, 0.0));

        for (Index i = 0; i < n_; ++i)
        {
            for (Index ptr = row_ptr_[i]; ptr < row_ptr_[i + 1]; ++ptr)
            {
                dense[i][col_ind_[ptr]] = values_[ptr];
            }
        }

        return dense;
    }

} // namespace hypergraph_ordering