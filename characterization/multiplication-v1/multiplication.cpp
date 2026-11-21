// Each build includes exactly one legacy definition. Never link these together.
#include LEGACY_HEADER
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>

int passed = 0;
int failed = 0;
void check(bool condition, const std::string& name) {
    std::cout << (condition ? "PASS " : "FAIL ") << name << '\n';
    condition ? ++passed : ++failed;
}

template<class T>
Matrix::Matrix<T> make(int rows, int cols, const T* values) {
    Matrix::Matrix<T> result(rows, cols);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c) result[r][c] = values[r * cols + c];
    return result;
}
template<class T>
bool equals(const Matrix::Matrix<T>& result, int rows, int cols, const T* values) {
    if (result.rows() != static_cast<size_t>(rows) || result.cols() != static_cast<size_t>(cols)) return false;
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
            if (result[r][c] != values[r * cols + c]) return false;
    return true;
}
template<class T>
void ordinary(const std::string& type) {
    // Pinned Neuro rectangular reference, extended to all three numeric types.
    const T av[] = {1,2,3,4,5,6}, bv[] = {7,8,9,10,11,12}, expected[] = {58,64,139,154};
    const auto a = make(2,3,av), b = make(3,2,bv);
    check(equals(a*b,2,2,expected) && equals(a,2,3,av) && equals(b,3,2,bv), type+" rectangular_2x3_3x2");
    // Pinned Bitmap square/int reference and Neuro square/float reference.
    const T sv[] = {1,2,3,4}, tv[] = {5,6,7,8}, se[] = {19,22,43,50};
    check(equals(make(2,2,sv)*make(2,2,tv),2,2,se),type+" square_2x2");
    // Asymmetric column-vector case exercises Neuro's distinct RHS.cols == 1 kernel.
    const T cv[] = {1,-2,3,4,-1,2,2,-3,0,7,-2,1}, dv[] = {2,-1,2,3}, ce[] = {22,-9,-8};
    check(equals(make(3,4,cv)*make(4,1,dv),3,1,ce),type+" column_3x4_4x1");
    bool threw = false;
    try { (void)(a * make(2,2,sv)); } catch (const std::invalid_argument&) { threw = true; }
    check(threw,type+" mismatch_invalid_argument");
}
template<class T>
void cancellation(const std::string& type) {
    static_assert(std::numeric_limits<T>::is_iec559);
    // Adapted from pinned Neuro cancellation reference, including the general kernel.
    const T big = T(1.0e20), small = T(3.25);
    const T av[] = {big,-big,small}, one[] = {1,1,1}, two[] = {1,2,1,2,1,2};
    const T single_expected[] = {small}, multi_expected[] = {small,T(6.5)};
    // Reassociation gives zero; this control ensures the fixture detects that regression.
    volatile T tail = -big + small;
    const bool sensitive = (big + tail != small);
    check(sensitive && equals(make(1,3,av)*make(3,1,one),1,1,single_expected),type+" cancellation_column_left_to_right");
    check(sensitive && equals(make(1,3,av)*make(3,2,two),1,2,multi_expected),type+" cancellation_general_left_to_right");
}
int main() {
    ordinary<float>("float"); ordinary<double>("double"); ordinary<int>("int");
    cancellation<float>("float"); cancellation<double>("double");
    std::cout << "RESULT passed=" << passed << " failed=" << failed << '\n';
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
