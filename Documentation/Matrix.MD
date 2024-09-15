
# Matrix Library Documentation

## Matrix Class

### Overview

The `Matrix` class is a templated data structure that supports various matrix operations such as addition, subtraction, multiplication, and division. It is designed for flexibility and supports any data type that allows arithmetic operations.

### Class Template

```cpp
template <typename T>
class Matrix;
```

### Constructors

#### `Matrix(int rows, int columns)`

Creates a new matrix with the specified number of rows and columns.

- **Parameters:**
  - `rows` - The number of rows in the matrix.
  - `columns` - The number of columns in the matrix.

### Methods

#### `int getRows() const`

Returns the number of rows in the matrix.

- **Return Value:** The number of rows.

#### `int getCols() const`

Returns the number of columns in the matrix.

- **Return Value:** The number of columns.

#### `T& operator()(int row, int column)`

Accesses the element at the specified position in the matrix.

- **Parameters:**
  - `row` - The row index of the element.
  - `column` - The column index of the element.
- **Return Value:** A reference to the matrix element at the specified position.

#### `Matrix<T> operator+(const Matrix<T>& other)`

Performs element-wise addition of two matrices if their dimensions match.

- **Parameters:**
  - `other` - The matrix to be added.
- **Return Value:** A new matrix containing the sum.

#### `Matrix<T> operator-(const Matrix<T>& other)`

Performs element-wise subtraction of two matrices if their dimensions match.

- **Parameters:**
  - `other` - The matrix to be subtracted.
- **Return Value:** A new matrix containing the difference.

#### `Matrix<T> operator*(const T scalar)`

Multiplies each element of the matrix by a scalar value.

- **Parameters:**
  - `scalar` - The scalar value.
- **Return Value:** A new matrix containing the scaled values.

#### `Matrix<T> operator/(const T scalar)`

Divides each element of the matrix by a scalar value.

- **Parameters:**
  - `scalar` - The scalar value.
- **Return Value:** A new matrix containing the divided values.

#### `Matrix<T> transpose() const`

Returns the transpose of the matrix.

- **Return Value:** A new matrix that is the transpose of the original.

### Iterators

The `Matrix` class provides custom iterators for traversing rows and columns.

#### `MatrixRowIterator<T> row_begin(int row)`

Returns an iterator to the beginning of the specified row.

- **Parameters:**
  - `row` - The row index.
- **Return Value:** An iterator to the beginning of the row.

#### `MatrixRowIterator<T> row_end(int row)`

Returns an iterator to the end of the specified row.

- **Parameters:**
  - `row` - The row index.
- **Return Value:** An iterator to the end of the row.

#### `MatrixColumnIterator<T> column_begin(int col)`

Returns an iterator to the beginning of the specified column.

- **Parameters:**
  - `col` - The column index.
- **Return Value:** An iterator to the beginning of the column.

#### `MatrixColumnIterator<T> column_end(int col)`

Returns an iterator to the end of the specified column.

- **Parameters:**
  - `col` - The column index.
- **Return Value:** An iterator to the end of the column.

---

## MatrixRowIterator Class

### Overview

The `MatrixRowIterator` class is a custom iterator for row-wise traversal of matrix elements. It supports random access operations.

### Class Template

```cpp
template <typename T>
class MatrixRowIterator;
```

### Methods

#### `MatrixRowIterator(T* ptr)`

Constructor that initializes the iterator with a pointer to a matrix element.

- **Parameters:**
  - `ptr` - Pointer to the matrix element.

#### `MatrixRowIterator<T>& operator++()`

Moves the iterator to the next element in the row (prefix increment).

- **Return Value:** The iterator after incrementing.

#### `T& operator*()`

Dereferences the iterator to access the element it points to.

- **Return Value:** A reference to the matrix element.

#### `bool operator==(const MatrixRowIterator<T>& other) const`

Compares two iterators for equality.

- **Parameters:**
  - `other` - The iterator to compare.
- **Return Value:** `true` if both iterators point to the same element, otherwise `false`.

#### `bool operator!=(const MatrixRowIterator<T>& other) const`

Compares two iterators for inequality.

- **Parameters:**
  - `other` - The iterator to compare.
- **Return Value:** `true` if the iterators point to different elements, otherwise `false`.

---

## MatrixColumnIterator Class

### Overview

The `MatrixColumnIterator` class is a custom iterator for column-wise traversal of matrix elements. It supports random access operations.

### Class Template

```cpp
template <typename T>
class MatrixColumnIterator;
```

### Methods

#### `MatrixColumnIterator(T* ptr)`

Constructor that initializes the iterator with a pointer to a matrix element.

- **Parameters:**
  - `ptr` - Pointer to the matrix element.

#### `MatrixColumnIterator<T>& operator++()`

Moves the iterator to the next element in the column (prefix increment).

- **Return Value:** The iterator after incrementing.

#### `T& operator*()`

Dereferences the iterator to access the element it points to.

- **Return Value:** A reference to the matrix element.

#### `bool operator==(const MatrixColumnIterator<T>& other) const`

Compares two iterators for equality.

- **Parameters:**
  - `other` - The iterator to compare.
- **Return Value:** `true` if both iterators point to the same element, otherwise `false`.

#### `bool operator!=(const MatrixColumnIterator<T>& other) const`

Compares two iterators for inequality.

- **Parameters:**
  - `other` - The iterator to compare.
- **Return Value:** `true` if the iterators point to different elements, otherwise `false`.
