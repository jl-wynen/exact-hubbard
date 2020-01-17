#ifndef EXACT_HUBBARD_LINALG_HPP
#define EXACT_HUBBARD_LINALG_HPP


/** \file
 * \brief Include blaze headers for required vector/matrix types
 * and define convenience aliases.
 */

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>


using DMatrix = blaze::DynamicMatrix<double>;
using DSparseMatrix = blaze::CompressedMatrix<double>;

using DVector = blaze::DynamicVector<double>;
using IVector = blaze::DynamicVector<int>;

#endif //EXACT_HUBBARD_LINALG_HPP
