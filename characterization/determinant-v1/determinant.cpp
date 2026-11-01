// Each executable includes exactly one immutable legacy definition.
#include LEGACY_HEADER
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
int passed=0,failed=0;
void check(bool ok,const std::string& name) {
    std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n';ok?++passed:++failed;
}
template<class F> bool invalid(F f) {
    try { f(); } catch(const std::invalid_argument&) { return true; }
    return false;
}
template<class T> bool close(T value,T expected) {
    if constexpr(std::is_floating_point<T>::value)
        return std::isfinite(value) && std::abs(value-expected)<=T(32)*std::numeric_limits<T>::epsilon()*std::abs(expected);
    else return value==expected;
}
template<class T> void fixtures(const std::string& type) {
    using M=Matrix::Matrix<T>;
    M e;check(e.Determinant()==T(1),type+" empty_unit");
    M scalar(1,1);scalar[0][0]=T(-7);check(scalar.Determinant()==T(-7),type+" scalar");
    M two(2,2);two[0][0]=4;two[0][1]=7;two[1][0]=2;two[1][1]=6;
    check(two.Determinant()==T(10),type+" two_by_two");
    M three(3,3);const int values[9]={1,2,3,0,1,4,5,6,0};
    for(size_t r=0;r<3;++r)for(size_t c=0;c<3;++c)three[r][c]=T(values[3*r+c]);
    check(three.Determinant()==T(1),type+" three_by_three");
    // Triangular determinant is the product 2*(-3)*4*5=-120.
    M upper(4,4);const int diag[4]={2,-3,4,5};
    for(size_t r=0;r<4;++r)for(size_t c=0;c<4;++c)upper[r][c]=r==c?T(diag[r]):(r<c?T(2*r+c+1):T(0));
    check(close(upper.Determinant(),T(-120)),type+" upper_triangular_four");
    bool unchanged=true;
    for(size_t r=0;r<4;++r)for(size_t c=0;c<4;++c)
        unchanged=unchanged && upper[r][c]==(r==c?T(diag[r]):(r<c?T(2*r+c+1):T(0)));
    check(unchanged,type+" determinant_input_unchanged");
    M one(4,4),both(4,4);const size_t p1[4]={1,0,2,3},p2[4]={1,0,3,2};
    for(size_t r=0;r<4;++r)for(size_t c=0;c<4;++c){one[r][c]=T(c==p1[r]);both[r][c]=T(c==p2[r]);}
    check(one.Determinant()==T(-1),type+" one_swap_zero_initial_pivot");
    check(both.Determinant()==T(1),type+" two_swap_sign");
    M singular=upper;for(size_t c=0;c<4;++c)singular[3][c]=singular[2][c];
    check(singular.Determinant()==T(0),type+" duplicate_row_singular");
    M zero(4,4);check(zero.Determinant()==T(0),type+" all_zero_singular");
    M six(6,6),seven(7,7);
    for(size_t r=0;r<6;++r)for(size_t c=0;c<6;++c)six[r][c]=r==c?T(r+1):T(0);
    for(size_t r=0;r<7;++r)for(size_t c=0;c<7;++c)seven[r][c]=r==c?T(r+1):T(0);
    check(close(six.Determinant(),T(720)),type+" six_dimension_boundary");
#ifndef NEURO
    if constexpr(std::is_integral<T>::value)check(invalid([&]{seven.Determinant();}),type+" seven_integer_rejected");
    else check(close(seven.Determinant(),T(5040)),type+" seven_float_supported");
#else
    check(close(seven.Determinant(),T(5040)),type+" seven_recursive_supported");
#endif
    M wide(2,3),tall(3,2);
    check(invalid([&]{wide.Determinant();}),type+" wide_rejected");
    check(invalid([&]{tall.Determinant();}),type+" tall_rejected");
    if constexpr(std::is_floating_point<T>::value) {
        // Analytic block determinant (delta*0 - 1*1)*2*3=-6.
        // delta is nonzero, below the inverse threshold; determinant has no such cutoff.
        M tiny(4,4);tiny[0][0]=T(1e-12);tiny[0][1]=T(1);tiny[1][0]=T(1);tiny[1][1]=T(0);tiny[2][2]=T(2);tiny[3][3]=T(3);
        check(close(tiny.Determinant(),T(-6)),type+" tiny_pivot_requires_stable_elimination");
    }
}
int main(){fixtures<float>("float");fixtures<double>("double");fixtures<int>("int");
    std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;}
