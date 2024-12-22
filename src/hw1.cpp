#include "hw1.h"
#include <iostream>
#include <random>
namespace algebra {
Matrix zeros(size_t n, size_t m) {
  Matrix out(n, std::vector<double>(m));
  return out;
}
Matrix ones(size_t n, size_t m) {
  Matrix out(n, std::vector<double>(m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      out[i][j] = 1;
    }
  }
  return out;
}
Matrix random(size_t n, size_t m, double min, double max) {
  if (min >= max) {
    throw std::logic_error("min is bigger than max");
  }
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(min, max);
  Matrix out(n, std::vector<double>(m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      out[i][j] = dist(mt);
    }
  }
  return out;
}
void show(const Matrix &matrix) {
  for (auto const &row : matrix) {
    for (auto const &ele : row) {
      std::cout << " " << ele;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}
Matrix multiply(const Matrix &matrix, double c) {
  Matrix out(matrix.size(), std::vector<double>(matrix[0].size()));
  for (int i = 0; i < matrix.size(); i++) {
    for (int j = 0; j < matrix[0].size(); j++) {
      out[i][j] = matrix[i][j] * c;
    }
  }
  return out;
}
Matrix multiply(const Matrix &matrix1, const Matrix &matrix2) {
  if (matrix1.size() == 0 || matrix2.size() == 0) {
    Matrix out{};
    return out;
  }
  if (matrix1[0].size() != matrix2.size()) {
    throw std::logic_error("matrix1 and matrix2 doesn't match");
  }
  Matrix out(matrix1.size(), std::vector<double>(matrix2[0].size()));
  for (int i = 0; i < matrix1.size(); i++) {
    for (int j = 0; j < matrix2[0].size(); j++) {
      double ele = 0;
      for (int k = 0; k < matrix1[0].size(); k++) {
        ele += matrix1[i][k] * matrix2[k][j];
      }
      out[i][j] = ele;
    }
  }
  return out;
}
Matrix sum(const Matrix &matrix, double c) {
  if (matrix.size() == 0) {
    Matrix out{};
    return out;
  }
  Matrix out(matrix.size(), std::vector<double>(matrix[0].size()));
  for (int i = 0; i < matrix.size(); i++) {
    for (int j = 0; j < matrix[0].size(); j++) {
      out[i][j] = matrix[i][j] + c;
    }
  }
  return out;
}
Matrix sum(const Matrix &matrix1, const Matrix &matrix2) {
  if (matrix1.size() == 0 && matrix2.size() == 0) {
    Matrix out{};
    return out;
  }
  if (matrix1.size() != matrix2.size() ||
      matrix1[0].size() != matrix2[0].size()) {
    throw std::logic_error("matrix1 and matrix2 doesn't match");
  }
  Matrix out(matrix1.size(), std::vector<double>(matrix1[0].size()));
  for (int i = 0; i < matrix1.size(); i++) {
    for (int j = 0; j < matrix1[0].size(); j++) {
      out[i][j] = matrix1[i][j] + matrix2[i][j];
    }
  }
  return out;
}
Matrix transpose(const Matrix &matrix) {
  if (matrix.size() == 0) {
    Matrix out{};
    return out;
  }
  Matrix out(matrix[0].size(), std::vector<double>(matrix.size()));
  for (int i = 0; i < matrix.size(); i++) {
    for (int j = 0; j < matrix[0].size(); j++) {
      out[j][i] = matrix[i][j];
    }
  }
  return out;
}
Matrix minor(const Matrix &matrix, size_t n, size_t m) {
  Matrix out(matrix.size() - 1, std::vector<double>(matrix[0].size() - 1));
  int k = 0, l = 0;
  for (int i = 0; i < matrix.size(); i++) {
    if (i == n)
      continue;
    l = 0;
    for (int j = 0; j < matrix[0].size(); j++) {
      if (j == m)
        continue;
      out[k][l] = matrix[i][j];
      l += 1;
    }
    k += 1;
  }
  return out;
}
double determinant(const Matrix &matrix) { return 0; }
Matrix inverse(const Matrix &matrix) {
  Matrix out(matrix[0].size(), std::vector<double>(matrix.size()));
  return out;
}
Matrix concatenate(const Matrix &matrix1, const Matrix &matrix2, int axis) {
  Matrix out(matrix1[0].size(), std::vector<double>(matrix1.size()));
  return out;
}
Matrix ero_swap(const Matrix &matrix, size_t r1, size_t r2) {
  Matrix out(matrix.size(), std::vector<double>(matrix[0].size()));
  return out;
}
Matrix ero_multiply(const Matrix &matrix, size_t r, double c) {
  Matrix out(matrix.size(), std::vector<double>(matrix[0].size()));
  return out;
}
Matrix ero_sum(const Matrix &matrix, size_t r1, double c, size_t r2) {
  Matrix out(matrix.size(), std::vector<double>(matrix[0].size()));
  return out;
}
Matrix upper_triangular(const Matrix &matrix) {
  Matrix out(matrix.size(), std::vector<double>(matrix[0].size()));
  return out;
}
} // namespace algebra
