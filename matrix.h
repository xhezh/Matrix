#pragma once
#include <iostream>
#include <stdexcept>
#include <initializer_list>
class MatrixIsDegenerateError : public std::runtime_error {
public:
  MatrixIsDegenerateError() : std::runtime_error("MatrixIsDegenerateError") {
  }
};
class MatrixOutOfRange : public std::out_of_range {
public:
  MatrixOutOfRange() : std::out_of_range("MatrixOutOfRange") {
  }
};
template<typename T, size_t Rows, size_t Cols>
class Matrix {
public:
  T data[Rows][Cols];
  size_t RowsNumber() const {
    return Rows;
  }
  size_t ColumnsNumber() const {
    return Cols;
  }
  T& operator()(size_t row, size_t col) {
    return data[row][col];
  }
  const T& operator()(size_t row, size_t col) const {
    return data[row][col];
  }
  T& At(size_t row, size_t col) {
    if (row >= Rows || col >= Cols) {
      throw MatrixOutOfRange{};
    }
    return data[row][col];
  }
  const T& At(size_t row, size_t col) const {
    if (row >= Rows || col >= Cols) {
      throw MatrixOutOfRange{};
    }
    return data[row][col];
  }
  Matrix<T, Rows, Cols>& operator+=(const Matrix<T, Rows, Cols>& other) {
    for (size_t i = 0; i < Rows; ++i) {
      for (size_t j = 0; j < Cols; ++j) {
        data[i][j] += other.data[i][j];
      }
    }
    return *this;
  }
  Matrix<T, Rows, Cols>& operator-=(const Matrix<T, Rows, Cols>& other) {
    for (size_t i = 0; i < Rows; ++i) {
      for (size_t j = 0; j < Cols; ++j) {
        data[i][j] -= other.data[i][j];
      }
    }
    return *this;
  }
  Matrix<T, Rows, Cols>& operator*=(const T& scalar) {
    for (size_t i = 0; i < Rows; ++i) {
      for (size_t j = 0; j < Cols; ++j) {
        data[i][j] *= scalar;
      }
    }
    return *this;
  }
  Matrix<T, Rows, Cols>& operator/=(const T& scalar) {
    for (size_t i = 0; i < Rows; ++i) {
      for (size_t j = 0; j < Cols; ++j) {
        data[i][j] /= scalar;
      }
    }
    return *this;
  }
  template<size_t OtherCols>
  Matrix<T, Rows, OtherCols>& operator*=(const Matrix<T, Cols, OtherCols>& other) {
    *this = *this * other;
    return *this;
  }
  template<size_t OtherCols>
  Matrix<T, Rows, OtherCols> operator*(const Matrix<T, Cols, OtherCols>& other) const {
    Matrix<T, Rows, OtherCols> result{};
    for (size_t i = 0; i < Rows; ++i) {
      for (size_t j = 0; j < OtherCols; ++j) {
        for (size_t k = 0; k < Cols; ++k) {
          result.data[i][j] += data[i][k] * other.data[k][j];
        }
      }
    }
    return result;
  }
  friend Matrix<T, Rows, Cols> operator+(Matrix<T, Rows, Cols> lhs, const Matrix<T, Rows, Cols>& rhs) {
    lhs += rhs;
    return lhs;
  }
  friend Matrix<T, Rows, Cols> operator-(Matrix<T, Rows, Cols> lhs, const Matrix<T, Rows, Cols>& rhs) {
    lhs -= rhs;
    return lhs;
  }
  friend Matrix<T, Rows, Cols> operator*(Matrix<T, Rows, Cols> lhs, const T& scalar) {
    lhs *= scalar;
    return lhs;
  }
  friend Matrix<T, Rows, Cols> operator*(const T& scalar, Matrix<T, Rows, Cols> rhs) {
    rhs *= scalar;
    return rhs;
  }
  friend Matrix<T, Rows, Cols> operator/(Matrix<T, Rows, Cols> lhs, const T& scalar) {
    lhs /= scalar;
    return lhs;
  }
  friend bool operator==(const Matrix<T, Rows, Cols>& lhs, const Matrix<T, Rows, Cols>& rhs) {
    for (size_t i = 0; i < Rows; ++i) {
      for (size_t j = 0; j < Cols; ++j) {
        if (lhs.data[i][j] != rhs.data[i][j]) {
          return false;
        }
      }
    }
    return true;
  }
  friend bool operator!=(const Matrix<T, Rows, Cols>& lhs, const Matrix<T, Rows, Cols>& rhs) {
    return !(lhs == rhs);
  }
};
template<typename T, size_t Rows, size_t Cols>
std::istream& operator>>(std::istream& is, Matrix<T, Rows, Cols>& matrix) {
  for (size_t row = 0; row < Rows; row++) {
    for (size_t column = 0; column < Cols; column++) {
      is >> matrix(row, column);
    }
  }
  return is;
}
template<typename T, size_t Rows, size_t Cols>
std::ostream& operator<<(std::ostream& os, const Matrix<T, Rows, Cols>& matrix) {
  for (size_t row = 0; row < Rows; row++) {
    for (size_t column = 0; column < Cols; column++) {
      if (column > 0) {
        os << ' ';
      }
      os << matrix(row, column);
    }
    os << '\n';
  }
  return os;
}
template<typename T, size_t Rows, size_t Cols>
Matrix<T, Cols, Rows> GetTransposed(const Matrix<T, Rows, Cols>& matrix) {
  Matrix<T, Cols, Rows> result;
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Cols; ++j) {
      result.data[j][i] = matrix.data[i][j];
    }
  }
  return result;
}
#define MATRIX_SQUARE_MATRIX_IMPLEMENTED
template<typename T, size_t Size>
T Trace(const Matrix<T, Size, Size>& matrix) {
  T trace = 0;
  for (size_t i = 0; i < Size; ++i) {
    trace += matrix.data[i][i];
  }
  return trace;
}
template<typename T, size_t Size>
T Determinant(const Matrix<T, Size, Size>& matrix) {
  if constexpr (Size == 1) {
    return matrix.data[0][0];
  }
  else if constexpr (Size == 2) {
    return matrix.data[0][0] * matrix.data[1][1] - matrix.data[0][1] * matrix.data[1][0];
  }
  else {
    T det = 0;
    for (size_t p = 0; p < Size; ++p) {
      Matrix<T, Size - 1, Size - 1> sub_matrix;
      for (size_t i = 1; i < Size; ++i) {
        size_t sub_col = 0;
        for (size_t j = 0; j < Size; ++j) {
          if (j == p) {
            continue;
          }
          sub_matrix.data[i - 1][sub_col] = matrix.data[i][j];
          ++sub_col;
        }
      }
      det += (p % 2 == 0 ? 1 : -1) * matrix.data[0][p] * Determinant(sub_matrix);
    }
    return det;
  }
}
template<typename T, size_t Size>
Matrix<T, Size, Size> GetInversed(const Matrix<T, Size, Size>& matrix) {
  T det = Determinant(matrix);
  if (det == 0) {
    throw MatrixIsDegenerateError{};
  }
  Matrix<T, Size, Size> adj;
  if constexpr (Size == 1) {
    adj.data[0][0] = 1;
  }
  else {
    for (size_t i = 0; i < Size; ++i) {
      for (size_t j = 0; j < Size; ++j) {
        Matrix<T, Size - 1, Size - 1> sub_matrix;
        size_t sub_row = 0;
        for (size_t m = 0; m < Size; ++m) {
          if (m == i) {
            continue;
          }
          size_t sub_col = 0;
          for (size_t n = 0; n < Size; ++n) {
            if (n == j) {
              continue;
            }
            sub_matrix.data[sub_row][sub_col] = matrix.data[m][n];
            ++sub_col;
          }
          ++sub_row;
        }
        adj.data[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * Determinant(sub_matrix);
      }
    }
  }
  return adj / det;
}
template<typename T, size_t Size>
void Inverse(Matrix<T, Size, Size>& matrix) {
  matrix = GetInversed(matrix);
}
template<typename T, size_t Size>
void Transpose(Matrix<T, Size, Size>& matrix) {
  for (size_t i = 0; i < Size; ++i) {
    for (size_t j = i + 1; j < Size; ++j) {
      std::swap(matrix.data[i][j], matrix.data[j][i]);
    }
  }
}
