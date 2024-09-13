
# Matrix Library

## Overview

This C++ library provides a templated `Matrix` class for performing matrix operations. It supports common operations such as matrix addition, subtraction, scalar multiplication, and division. The class is designed for ease of use and efficiency, with custom iterators for traversing rows and columns of the matrix.

## Features

- **Templated Matrix Class**: The `Matrix` class is templated, allowing you to create matrices with different data types (e.g., `int`, `float`, `double`, etc.).
- **Matrix Arithmetic**: Support for addition, subtraction, multiplication, and division of matrices and scalars.
- **Custom Iterators**: The library provides custom row and column iterators, allowing easy traversal and manipulation of matrix elements.
- **Efficient Access**: Provides direct access to matrix elements with operator overloading for easy indexing.

## Usage

### Creating a Matrix

To create a matrix, simply instantiate the `Matrix` class with the desired data type and dimensions.

```cpp
#include "matrix.h"

int main() {
    // Create a 3x3 matrix of integers
    Matrix<int> mat(3, 3);
    
    // Assign values to the matrix
    mat[0][0] = 1;
    mat[0][1] = 2;
    mat[0][2] = 3;

    // Output the value at (0, 0)
    std::cout << mat[0][0] << std::endl;
    
    return 0;
}
```

### Matrix Arithmetic

You can easily perform arithmetic operations on matrices:

```cpp
#include "matrix.h"

int main() {
    Matrix<int> mat1(3, 3);
    Matrix<int> mat2(3, 3);
    
    // Initialize matrices with some values
    // ...

    // Add two matrices
    Matrix<int> result = mat1 + mat2;
    
    // Subtract two matrices
    result = mat1 - mat2;
    
    // Multiply a matrix by a scalar
    result = mat1 * 2;
    
    return 0;
}
```

### Iterating Over Rows and Columns

The library provides custom iterators for rows and columns, allowing you to iterate through matrix elements:

```cpp
#include "matrix.h"

int main() {
    Matrix<int> mat(3, 3);

    // Fill the matrix with values
    // ...

    // Iterate over rows
    for (auto rowIt = mat.row_begin(0); rowIt != mat.row_end(0); ++rowIt) {
        std::cout << *rowIt << " ";
    }
    
    // Iterate over columns
    for (auto colIt = mat.column_begin(0); colIt != mat.column_end(0); ++colIt) {
        std::cout << *colIt << " ";
    }

    return 0;
}
```

### Error Handling

Ensure that matrices being operated on have compatible dimensions. For instance, attempting to add or subtract matrices of different sizes will result in an error.

## Installation

To use this library in your project, simply include the `matrix.h` header file:

```cpp
#include "matrix.h"
```

Make sure your project is set up to compile C++11 or later to use the iterator functionalities provided.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributions

Contributions to this project are welcome. Feel free to open issues or submit pull requests.
