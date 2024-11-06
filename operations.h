#pragma once
#include <vector>

using std::vector;

template<typename T>
vector<vector<T>> operator *(const vector<vector<T>> matrix, T variable) // vector<vector<T>> * T умножение матрицы на число
{
    vector<vector<T>> result(matrix.size());
    for (size_t i = 0; i < matrix.size(); i++)
    {
        result[i].resize(matrix[i].size());
        for (size_t j = 0; j < matrix[i].size(); j++)
        {
            result[i][j] = matrix[i][j] * variable;
        }
    }
    return result;
}

template<typename T>
vector<T> operator *(const vector<vector<T>> matrix, const vector<T> vec) // vector<vector<T>> * vector<T> умножение матрицы на вектор
{
    vector<T> result(matrix.size());
    for (size_t i = 0; i < matrix.size(); i++)
    {
        for (size_t j = 0; j < matrix[i].size(); j++)
        {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}

template<typename T>
vector<T> operator *(const vector<T> vec, const T variable) // vector<vector<T>> * vector<T> умножение матрицы на число
{
    vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); i++)
    {
        result[i] = vec[i] * variable;
    }
    return result;
}

template<typename T>
vector<vector<T>> operator +(const vector<vector<T>> matrix1, const vector<vector<T>> matrix2)
{
    vector<vector<T>> result(matrix1.size());
    for (size_t i = 0; i < matrix1.size(); i++)
    {
        result[i].resize(matrix1[i].size());
        for (size_t j = 0; j < matrix1[i].size(); j++)
        {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return result;
}