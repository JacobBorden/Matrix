#pragma once
#include <cstddef>
#include <memory>
#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <time.h>
#include <iterator>
namespace Matrix
{
template <typename MatrixRow>
class MatrixRowIterator {
public:
    // Member type definitions to conform to iterator requirements
    using value_type = typename MatrixRow::value_type; // Type of elements the iterator refers to
    using pointer = value_type*; // Pointer to the element type
    using reference = value_type&; // Reference to the element type
    using iterator_category = std::random_access_iterator_tag; // Iterator category to support random access
    using difference_type = std::ptrdiff_t; // Type to express the difference between two iterators

    // Constructor initializes the iterator with a pointer to a matrix row element
    MatrixRowIterator(pointer ptr) : m_ptr(ptr) {}

    // Pre-increment operator advances the iterator to the next element and returns a reference to the updated iterator
    MatrixRowIterator& operator++() {
        m_ptr++;
        return *this;
    }

    // Post-increment operator advances the iterator to the next element and returns the iterator before advancement
    MatrixRowIterator operator++(int) {
        MatrixRowIterator it = *this;
        ++*this;
        return it;
    }

    // Addition operator returns a new iterator advanced by 'n' positions
    MatrixRowIterator operator+(difference_type n) const { return MatrixRowIterator(m_ptr + n); }

    // Compound addition operator advances the iterator by 'n' positions and returns a reference to the updated iterator
    MatrixRowIterator& operator+=(difference_type n) {
        m_ptr += n;
        return *this;
    }

    // Pre-decrement operator moves the iterator to the previous element and returns a reference to the updated iterator
    MatrixRowIterator& operator--() {
        m_ptr--;
        return *this;
    }

    // Post-decrement operator moves the iterator to the previous element and returns the iterator before movement
    MatrixRowIterator operator--(int) {
        MatrixRowIterator it = *this;
        --*this;
        return it;
    }

    // Subtraction operator returns a new iterator moved back by 'n' positions
    MatrixRowIterator operator-(difference_type n) const { return MatrixRowIterator(m_ptr - n); }

    // Compound subtraction operator moves the iterator back by 'n' positions and returns a reference to the updated iterator
    MatrixRowIterator& operator-=(difference_type n) {
        m_ptr -= n;
        return *this;
    }

    // Subtraction operator calculates the difference between two iterators
    difference_type operator-(const MatrixRowIterator& other) const { return m_ptr - other.m_ptr; }

    // Arrow operator provides access to the element's members the iterator points to
    pointer operator->() const { return m_ptr; }

    // Dereference operators return a (const) reference to the element the iterator points to
    reference operator*() { return *m_ptr; }
    const reference operator*() const { return *m_ptr; }

    // Comparison operators for equality and inequality checks between iterators
    bool operator==(const MatrixRowIterator& other) const { return m_ptr == other.m_ptr; }
    bool operator!=(const MatrixRowIterator& other) const { return m_ptr != other.m_ptr; }

    // Relational operators compare the positions of two iterators
    bool operator<(const MatrixRowIterator& other) const { return m_ptr < other.m_ptr; }
    bool operator<=(const MatrixRowIterator& other) const { return m_ptr <= other.m_ptr; }
    bool operator>(const MatrixRowIterator& other) const { return m_ptr > other.m_ptr; }
    bool operator>=(const MatrixRowIterator& other) const { return m_ptr >= other.m_ptr; }

    // Subscript operator provides random access to elements relative to the current iterator position
    reference operator[](difference_type n) const { return *(*this + n); }

private:
    pointer m_ptr; // Internal pointer to the current element
};
//--------------------------------------------------------------------------

template <typename T>
class MatrixColumnIterator {
public:
    // Type aliases for iterator traits
    using value_type = T; // Type of elements the iterator can dereference
    using pointer = T*; // Pointer to the element type
    using reference = T&; // Reference to the element type
    using iterator_category = std::random_access_iterator_tag; // Iterator category defining the capabilities of the iterator
    using difference_type = std::ptrdiff_t; // Type to express the difference between two iterators

    // Constructor initializes the iterator with a pointer to a matrix element and the total number of columns in the matrix
    MatrixColumnIterator(pointer ptr, size_t totalColumns) : m_ptr(ptr), m_totalColumns(totalColumns) {}

    // Pre-increment operator advances the iterator to the next element in the column and returns a reference to the updated iterator
    MatrixColumnIterator& operator++() {
        m_ptr += m_totalColumns; // Move pointer down one row in the current column
        return *this;
    }

    // Post-increment operator advances the iterator to the next element in the column and returns the iterator before the increment
    MatrixColumnIterator operator++(int) {
        MatrixColumnIterator it = *this; // Make a copy of the current iterator
        m_ptr += m_totalColumns; // Move pointer down one row in the current column
        return it; // Return the copy representing the iterator before increment
    }

    // Addition operator returns a new iterator advanced by 'n' positions in the column
    MatrixColumnIterator operator+(difference_type n) const {
        return MatrixColumnIterator(m_ptr + (n * m_totalColumns), m_totalColumns); // Calculate new position and create a new iterator
    }

    // Compound addition operator advances the iterator by 'n' positions in the column and returns a reference to the updated iterator
    MatrixColumnIterator& operator+=(difference_type n) {
        m_ptr += (n * m_totalColumns); // Adjust pointer by 'n' rows down in the current column
        return *this;
    }

    // Pre-decrement operator moves the iterator to the previous element in the column and returns a reference to the updated iterator
    MatrixColumnIterator& operator--() {
        m_ptr -= m_totalColumns; // Move pointer up one row in the current column
        return *this;
    }

    // Post-decrement operator moves the iterator to the previous element in the column and returns the iterator before the decrement
    MatrixColumnIterator operator--(int) {
        MatrixColumnIterator it = *this; // Make a copy of the current iterator
        m_ptr -= m_totalColumns; // Move pointer up one row in the current column
        return it; // Return the copy representing the iterator before decrement
    }

    // Subtraction operator returns a new iterator moved back by 'n' positions in the column
    MatrixColumnIterator operator-(difference_type n) const {
        return MatrixColumnIterator(m_ptr - (n * m_totalColumns), m_totalColumns); // Calculate new position and create a new iterator
    }

    // Compound subtraction operator moves the iterator back by 'n' positions in the column and returns a reference to the updated iterator
    MatrixColumnIterator& operator-=(difference_type n) {
        m_ptr -= (n * m_totalColumns); // Adjust pointer by 'n' rows up in the current column
        return *this;
    }

    // Subtraction operator calculates the difference between two iterators in terms of column positions
    difference_type operator-(const MatrixColumnIterator& other) const {
        return (m_ptr - other.m_ptr) / m_totalColumns; // Calculate element-wise distance between iterators
    }

    // Comparison operators for checking equality and inequality between iterators
    bool operator==(const MatrixColumnIterator& other) const { return m_ptr == other.m_ptr; }
    bool operator!=(const MatrixColumnIterator& other) const { return m_ptr != other.m_ptr; }

    // Relational operators for ordering iterators
    bool operator<(const MatrixColumnIterator& other) const { return m_ptr < other.m_ptr; }
    bool operator<=(const MatrixColumnIterator& other) const { return m_ptr <= other.m_ptr; }
    bool operator>(const MatrixColumnIterator& other) const { return m_ptr > other.m_ptr; }
    bool operator>=(const MatrixColumnIterator& other) const { return m_ptr >= other.m_ptr; }

    // Dereference operator provides access to the current element the iterator points to
    reference operator*() const { return *m_ptr; }

    // Member access operator allows access to the element's members
    pointer operator->() const { return m_ptr; }

    // Subscript operator provides random access to elements relative to the current iterator position
    reference operator[](difference_type n) const { return *(*this + n); }

private:
    pointer m_ptr; // Pointer to the current element in the matrix
    size_t m_totalColumns; // Total number of columns in the matrix, used for column-wise navigation
};

//--------------------------------------------------------------------------
    template <typename Matrix>
    class MatrixIterator
    {
    public:
        using value_type = typename Matrix::value_type;
        using pointer = value_type *;
        using reference = value_type &;

        MatrixIterator(pointer ptr) : m_ptr(ptr) {}

        MatrixIterator &operator++()
        {
            m_ptr++;
            return *this;
        }

        MatrixIterator operator++(int)
        {
            MatrixIterator it = *this;
            ++(this);
            return it;
        }

        MatrixIterator &operator--()
        {
            m_ptr--; 
	    return *this;
        }

        MatrixIterator operator--(int)
        {
            MatrixIterator it = *this;
            --(this);
            return it;
        }

        pointer operator->()
        {
            return m_ptr;
        }

        reference operator*()
        {
            return *m_ptr;
        }

        bool operator==(MatrixIterator other)
        {
            return this->m_ptr == other.m_ptr;
        }

        bool operator!=(MatrixIterator other)
        {
            return this->m_ptr != other.m_ptr;
        }

    private:
        pointer m_ptr;
    };
//-------------------------------------------------------------------------
    template <typename T>
    class MatrixRow
    {
    public:
        using value_type = T;
        using Iterator = MatrixRowIterator<MatrixRow<T>>;

    public:
        MatrixRow();
        MatrixRow(size_t);
        ~MatrixRow();
        void resize(size_t size);
        void assign(size_t size, T val);
        size_t size();
        size_t capacity();
        T &operator[](int i);
        Iterator begin() { return Iterator(m_Data); }
        Iterator end() { return Iterator(m_Data + m_Size); }
	Iterator begin() const {return Iterator(m_Data);}
	Iterator end() const {return Iterator(m_Data + m_Size);}
    private:
        size_t m_Size;
        size_t m_Capacity;
        //T *m_Data;
	std::unique_ptr<T[]> m_Data;
    };

    template <typename T>
    inline MatrixRow<T>::MatrixRow()
    {
        m_Data = new T[0];
        m_Size = 0;
        m_Capacity = sizeof(T) * m_Size;
    }

    template <typename T>
    inline MatrixRow<T>::MatrixRow(size_t size)
    {
        m_Data = new T[size];
        m_Size = size;
        m_Capacity = sizeof(T) * m_Size;
    }
	
    template <typename T>
    inline MatrixRow<T>::~MatrixRow()
    {
	    delete[] m_Data;
    }

    template <typename T>
    inline void MatrixRow<T>::resize(size_t newSize)
    {
	    T* newBlock = new T[newSize];
	    size_t minSize = std::min(newSize, m_Size);
	    for(size_t i = 0; i < minSize; ++i){
		    newBlock[i] = m_Data[i];
	    }
	    for (size_t i = minSize; i < newSize; ++i){
		    newBlock[i] = T();
	    }
	    delete[] m_Data;
	    m_Data = newBlock;
	    m_Size = newSize;
	    m_Capacity = newSize;
    }

    template <typename T>
    inline void MatrixRow<T>::assign(size_t size, T val)
    {
        resize(size);
        std::fill(begin(), end(), val);
    }

    template <typename T>
    inline size_t MatrixRow<T>::size()
    {
        return m_Size;
    }

    template <typename T>
    inline size_t MatrixRow<T>::capacity()
    {
        return m_Capacity;
    }

    template <typename T>
    inline T &MatrixRow<T>::operator[](int i)
    {
        return m_Data[i];
    }
//-------------------------------------------------------------------------
    template <typename T>
    class Matrix
    {
    public:
        using value_type = MatrixRow<T>;
        using Iterator = MatrixIterator<Matrix<T>>;

    public:
        Matrix<T>();
        Matrix<T>(int row_count, int column_count);
        ~Matrix<T>();
        int size();
        int rows();
        int cols();
        size_t capacity();
        void resize(int row_count, int col_count);
        void assign(int row_count, int col_count, T val);
        void assign(T val);
        Matrix<T> MergeVertical(Matrix<T> &b);
        Matrix<T> MergeHorizontal(Matrix<T> &b);
        Matrix<T> *SplitVertical();
        Matrix<T> *SplitVertical(int num);
        Matrix<T> *SplitHorizontal();
        Matrix<T> *SplitHorizontal(int num);
        Matrix<T> SigmoidMatrix();
        Matrix<T> Randomize();
        Matrix<T> CreateIdentityMatrix();
        Matrix<T> ZeroMatrix();
        MatrixRow<T> &operator[](int i);
        Matrix<T> operator*(Matrix<T> &b);
        Matrix<T> operator*(T b);
        Matrix<T> operator*=(T b);
        Matrix<T> operator+(Matrix<T> b);
        Matrix<T> operator+(T b);
        Matrix<T> operator+=(Matrix<T> b);
        Matrix<T> operator+=(T b);
        Matrix<T> operator-(Matrix<T> b);
        Matrix<T> operator-(T b);
        Matrix<T> operator-=(Matrix<T> b);
        Matrix<T> operator-=(T b);
        Matrix<T> operator/(T b);
        Matrix<T> operator/=(T b);
        Iterator begin() { return Iterator(m_Data); }
        Iterator end() { return Iterator(m_Data + m_Rows); }

    private:
        int m_Rows;
        int m_Cols;
        size_t m_Size;
        size_t m_Capacity;
        MatrixRow<T> *m_Data;
    };

    template <typename T>
    inline Matrix<T>::Matrix()
    {
        m_Rows = 0;
        m_Cols = 0;
        m_Size = m_Rows * m_Cols;
        m_Capacity = m_Size * sizeof(T);
        m_Data = new MatrixRow<T>[m_Rows];
        for (int i = 0; i < m_Rows; i++)
            m_Data[i] = MatrixRow<T>(m_Cols);
    }

    template <typename T>
    inline Matrix<T>::Matrix(int row_count, int col_count)
    {
        m_Rows = row_count;
        m_Cols = col_count;
        m_Size = m_Rows * m_Cols;
        m_Capacity = m_Size * sizeof(T);
        m_Data = new MatrixRow<T>[m_Rows];
        for (int i = 0; i < m_Rows; i++)
            m_Data[i] = MatrixRow<T>(m_Cols);
    }

    template <typename T>
    inline Matrix<T>::~Matrix()
    {
    }

    template <typename T>
    inline int Matrix<T>::size()
    {
        return m_Size;
    }

    template <typename T>
    inline int Matrix<T>::rows()
    {
        return m_Rows;
    }

    template <typename T>
    inline int Matrix<T>::cols()
    {
        return m_Cols;
    }

    template <typename T>
    inline size_t Matrix<T>::capacity()
    {
        return m_Capacity;
    }

    template <typename T>
    inline void Matrix<T>::resize(int row_count, int col_count)
    {
        MatrixRow<T> *newBlock = new MatrixRow<T>[row_count];
        for (int i = 0; i < row_count; i++)
            newBlock[i] = NULL;
        if (row_count > m_Rows)
            for (int i = 0; i < m_Rows; i++)
                newBlock[i] = std::move(m_Data[i]);
        else
            newBlock = std::move(m_Data);
        m_Data = newBlock;
        m_Rows = row_count;
        m_Cols = col_count;
        m_Size = m_Rows * m_Cols;
        m_Capacity = m_Size * sizeof(T);
        for (auto &i : *this)
            i.resize(m_Cols);
    }

    template <typename T>
    inline void Matrix<T>::assign(int row_count, int col_count, T val)
    {
        this->resize(row_count, col_count);
        for (auto i : *this)
            for (auto &j : i)
                j = val;
    }

    template <typename T>
    inline void Matrix<T>::assign(T val)
    {
        for (auto i : *this)
            for (auto &j : i)
                j = val;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::MergeVertical(Matrix<T> &b)
    {
        if (m_Cols == b.m_Cols)
        {
            this->resize(m_Rows * 2, m_Cols);
            for (int i = (m_Rows / 2); i < m_Rows; i++)
                m_Data[i] = std::move(b.m_Data[i - (m_Rows / 2)]);
        }
        return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::MergeHorizontal(Matrix<T> &b)
    {
        if (m_Rows == b.m_Rows)
        {
            this->resize(m_Rows, m_Cols * 2);
            for (int i = 0; i < m_Rows; i++)
                for (int j = (m_Cols / 2); j < m_Cols; j++)
                    m_Data[i][j] = std::move(b.m_Data[i][j - (m_Cols / 2)]);
        }
        return *this;
    }

    template <typename T>
    inline Matrix<T> *Matrix<T>::SplitVertical()
    {
        Matrix<T> *split_matrix = new Matrix<T>[2];
        split_matrix[0].resize(m_Rows / 2, m_Cols);
        split_matrix[1].resize(m_Rows / 2, m_Cols);
        for (int i = 0; i < m_Rows / 2; i++)
            split_matrix[0][i] = std::move(m_Data[i]);
        for (int i = m_Rows / 2; i < m_Rows; i++)
            split_matrix[1][i - (m_Rows / 2)] = std::move(m_Data[i]);
        return split_matrix;
    }

    template <typename T>
    inline Matrix<T> *Matrix<T>::SplitVertical(int num)
    {
        Matrix<T> *split_matrix = new Matrix<T>[num];

        for (int i = 0; i < num; i++)
            split_matrix[i].resize(m_Rows / num, m_Cols);

        for (int j = 0; j < num; j++)
            for (int i = 0; i < m_Rows / num; i++)
                split_matrix[j][i] = std::move(m_Data[(j * (m_Rows / num)) + i]);

        return split_matrix;
    }

    template <typename T>
    inline Matrix<T> *Matrix<T>::SplitHorizontal()
    {
        Matrix<T> *split_matrix = new Matrix<T>[2];
        split_matrix[0].resize(m_Rows, m_Cols / 2);
        split_matrix[1].resize(m_Rows, m_Cols / 2);
        for (int i = 0; i < m_Rows; i++)
        {
            for (int j = 0; j < m_Cols / 2; j++)
                split_matrix[0][i][j] = std::move(m_Data[i][j]);
            for (int j = m_Cols / 2; j < m_Cols; j++)
                split_matrix[0][i][j - (m_Cols / 2)] = std::move(m_Data[i][j]);
        }
        return split_matrix;
    }

    template <typename T>
    inline Matrix<T> *Matrix<T>::SplitHorizontal(int num)
    {
        Matrix<T> *split_matrix = new Matrix<T>[num];
        for (int i = 0; i < num; i++)
            split_matrix[i].resize(m_Rows, m_Cols / num);

        for (int i = 0; i < num; i++)
            for (int j = 0; j < m_Rows; j++)
                for (int k = 0; k < m_Cols / num; k++)
                    split_matrix[i][j][k] = std::move(m_Data[j][(i * (m_Cols / num)) + k]);

        return split_matrix;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::SigmoidMatrix()
    {
        for (auto i : *this)
            for (auto &j : i)
                j = 1 / (1 + exp(j));

        return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::Randomize()
    {
        srand(time(0));
        for (auto i : *this)
            for (auto &j : i)
            {
                int num = rand() % 200 + 1 - 100;
                j = num / 100.0;
            }

        return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::CreateIdentityMatrix()
    {
        if (m_Rows == m_Cols)
        {
            this->ZeroMatrix();
            for (int i = 0; i < m_Rows; i++)
                m_Data[i][i] = 1;
        }
        return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::ZeroMatrix()
    {
        this->assign(0);
        return *this;
    }

    template <typename T>
    inline MatrixRow<T> &Matrix<T>::operator[](int i)
    {
        return m_Data[i];
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator+(Matrix<T> b)
    {
        if ((m_Rows == b.m_Rows) && (m_Cols == b.m_Cols))
        {
            Matrix<T> c(m_Rows, m_Cols);
            for (int i = 0; i < m_Rows; i++)
                for (int j = 0; j < m_Cols; j++)
                    c[i][j] = m_Data[i][j] + b[i][j];
            return c;
        }

        else
            return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator+(T b)
    {
        Matrix<T> c(m_Rows, m_Cols);
        for (int i = 0; i < m_Rows; i++)
            for (int j = 0; j < m_Cols; j++)
                c[i][j] = m_Data[i][j] + b;
        return c;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator+=(Matrix<T> b)
    {
        if ((m_Rows == b.m_Rows) && (m_Cols == b.m_Cols))
        {
            Matrix<T> c(m_Rows, m_Cols);
            for (int i = 0; i < m_Rows; i++)
                for (int j = 0; j < m_Cols; j++)
                    c[i][j] = m_Data[i][j] + b[i][j];
            *this = c;
        }

        return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator+=(T b)
    {
        Matrix<T> c(m_Rows, m_Cols);
        for (int i = 0; i < m_Rows; i++)
            for (int j = 0; j < m_Cols; j++)
                c[i][j] = m_Data[i][j] + b;
        *this = c;
        return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator*(Matrix<T> &b)
    {
        if (m_Cols == b.m_Rows)
        {
            Matrix<T> c(m_Rows, b.m_Cols);
            for (int i = 0; i < m_Rows; i++)
                for (int j = 0; j < m_Cols; j++)
                    for (int k = 0; k < b.m_Cols; k++)
                        c[i][k] += m_Data[i][j] * b[j][k];
            return c;
        }

        else
            return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator*(T b)
    {
        Matrix<T> c(m_Rows, m_Cols);
        for (int i = 0; i < m_Rows; i++)
            for (int j = 0; j < m_Cols; j++)
                c[i][j] = m_Data[i][j] * b;
        return c;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator*=(T b)
    {
        Matrix<T> c(m_Rows, m_Cols);
        for (int i = 0; i < m_Rows; i++)
            for (int j = 0; j < m_Cols; j++)
                c[i][j] = m_Data[i][j] * b;
        *this = c;
        return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator-(Matrix<T> b)
    {
        if ((m_Rows == b.m_Rows) && (m_Cols == b.m_Cols))
        {
            Matrix<T> c(m_Rows, m_Cols);
            for (int i = 0; i < m_Rows; i++)
                for (int j = 0; j < m_Cols; j++)
                    c[i][j] = m_Data[i][j] - b[i][j];
            return c;
        }

        else
            return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator-(T b)
    {
        Matrix<T> c(m_Rows, m_Cols);
        for (int i = 0; i < m_Rows; i++)
            for (int j = 0; j < m_Cols; j++)
                c[i][j] = m_Data[i][j] - b;
        return c;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator-=(Matrix<T> b)
    {
        if ((m_Rows == b.m_Rows) && (m_Cols == b.m_Cols))
        {
            Matrix<T> c(m_Rows, m_Cols);
            for (int i = 0; i < m_Rows; i++)
                for (int j = 0; j < m_Cols; j++)
                    c[i][j] = m_Data[i][j] - b[i][j];
            *this = c;
        }

        return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator-=(T b)
    {
        Matrix<T> c(m_Rows, m_Cols);
        for (int i = 0; i < m_Rows; i++)
            for (int j = 0; j < m_Cols; j++)
                c[i][j] = m_Data[i][j] - b;
        *this = c;
        return *this;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator/(T b)
    {
        Matrix<T> c(m_Rows, m_Cols);
        for (int i = 0; i < m_Rows; i++)
            for (int j = 0; j < m_Cols; j++)
                c[i][j] = m_Data[i][j] / b;
        return c;
    }

    template <typename T>
    inline Matrix<T> Matrix<T>::operator/=(T b)
    {
        Matrix<T> c(m_Rows, m_Cols);
        for (int i = 0; i < m_Rows; i++)
            for (int j = 0; j < m_Cols; j++)
                c[i][j] = m_Data[i][j] / b;
        *this = c;
        return *this;
    }
}
#endif
