#ifndef EXACT_HUBBARD_MATRIX_HPP
#define EXACT_HUBBARD_MATRIX_HPP


#include <cassert>
#include <vector>


template <typename T>
class Matrix
{
    using ValueType = T;


    std::vector<T> storage_;
    std::size_t nrow_;
    std::size_t ncol_;


public:
    Matrix(std::size_t const nrow, std::size_t const ncol)
            : storage_(nrow*ncol), nrow_{nrow}, ncol_{ncol}
    {
        assert(nrow > 0);
        assert(ncol > 0);
    }


    T *data() noexcept
    {
        return storage_.data();
    }


    T const *data() const noexcept
    {
        return storage_.data();
    }


    [[nodiscard]] std::size_t nrow() const noexcept
    {
        return nrow_;
    }


    [[nodiscard]] std::size_t ncol() const noexcept
    {
        return ncol_;
    }


    std::size_t totalIndex(std::size_t const row, std::size_t const col) noexcept
    {
        return row * ncol_ + col;
    }


    T &operator()(std::size_t const row, std::size_t const col) noexcept
    {
        assert(row < nrow_);
        assert(col < ncol_);
        return storage_[totalIndex(row, col)];
    }


    T const &operator()(std::size_t const row, std::size_t const col) const noexcept
    {
        assert(row < nrow_);
        assert(col < ncol_);
        return storage_[totalIndex(row, col)];
    }
};


using DMatrix = Matrix<double>;

#endif //EXACT_HUBBARD_MATRIX_HPP
