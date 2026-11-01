// One immutable legacy Matrix definition per executable.
#include LEGACY_HEADER
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
int passed=0,failed=0;
void check(bool ok,const std::string& name) {
    std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n'; ok?++passed:++failed;
}
template<class F> bool invalid(F f) {
    try { f(); } catch(const std::invalid_argument&) { return true; }
    return false;
}
template<class T> void splits(const std::string& type) {
    using M=Matrix::Matrix<T>; M a(4,6);
    for(size_t r=0;r<4;++r) for(size_t c=0;c<6;++c) a[r][c]=T(10*r+c);
    auto vertical=a.SplitVertical(2); bool good=vertical.size()==2;
    for(size_t s=0;s<vertical.size();++s) {
        good=good && vertical[s].rows()==2 && vertical[s].cols()==6;
        for(size_t r=0;r<2;++r) for(size_t c=0;c<6;++c)
            good=good && vertical[s][r][c]==T(10*(2*s+r)+c);
    }
    check(good,type+" vertical_values");
    vertical[0][0][0]=T(99);
    check(a[0][0]==T(0) && vertical[1][0][0]==T(20),type+" vertical_independent");
    auto horizontal=a.SplitHorizontal(3); good=horizontal.size()==3;
    for(size_t s=0;s<horizontal.size();++s) {
        good=good && horizontal[s].rows()==4 && horizontal[s].cols()==2;
        for(size_t r=0;r<4;++r) for(size_t c=0;c<2;++c)
            good=good && horizontal[s][r][c]==T(10*r+2*s+c);
    }
    check(good,type+" horizontal_values");
    horizontal[0][0][0]=T(99);
    check(a[0][0]==T(0) && horizontal[1][0][0]==T(2),type+" horizontal_independent");
    auto vd=a.SplitVertical(); auto hd=a.SplitHorizontal();
    check(vd.size()==2 && vd[1][0][0]==T(20) && hd.size()==2 && hd[1][0][0]==T(3),type+" default_halves");
    check(invalid([&]{a.SplitVertical(0);}),type+" vertical_zero_rejected");
    check(invalid([&]{a.SplitHorizontal(0);}),type+" horizontal_zero_rejected");
    check(invalid([&]{a.SplitVertical(3);}),type+" vertical_nondivisor_rejected");
    check(invalid([&]{a.SplitHorizontal(4);}),type+" horizontal_nondivisor_rejected");
    M odd(3,3);
    check(invalid([&]{odd.SplitVertical();}) && invalid([&]{odd.SplitHorizontal();}),type+" odd_default_rejected");
    M e;
#ifdef NEURO
    check(e.SplitVertical().empty() && e.SplitHorizontal().empty(),type+" empty_default_no_parts");
    auto ev=e.SplitVertical(3); auto eh=e.SplitHorizontal(3);
    good=ev.size()==3 && eh.size()==3;
    for(const auto& m:ev) good=good && m.rows()==0 && m.cols()==0;
    for(const auto& m:eh) good=good && m.rows()==0 && m.cols()==0;
    check(good,type+" empty_explicit_three_empty_parts");
#else
    check(invalid([&]{e.SplitVertical();}) && invalid([&]{e.SplitHorizontal();}),type+" empty_default_rejected");
    check(invalid([&]{e.SplitVertical(3);}) && invalid([&]{e.SplitHorizontal(3);}),type+" empty_explicit_rejected");
#endif
    check(invalid([&]{e.SplitVertical(0);}) && invalid([&]{e.SplitHorizontal(0);}),type+" empty_zero_rejected");
}
template<class T> bool sigmoid_values(const Matrix::Matrix<T>& m) {
    const T expected[6]={T(0.11920292202211755),T(0.2689414213699951),T(0.5),T(0.7310585786300049),T(0.8807970779778823),T(0.9525741268224334)};
    if(m.rows()!=2 || m.cols()!=3) return false;
    for(size_t r=0;r<2;++r) for(size_t c=0;c<3;++c)
        if(!std::isfinite(m[r][c]) || std::abs(m[r][c]-expected[3*r+c])>T(0.000001)) return false;
    return true;
}
template<class T> void sigmoid_checker_regressions(const std::string& type) {
    using M=Matrix::Matrix<T>;
    const T expected[6]={T(0.11920292202211755),T(0.2689414213699951),T(0.5),T(0.7310585786300049),T(0.8807970779778823),T(0.9525741268224334)};
    M a(2,3);
    for(size_t r=0;r<2;++r) for(size_t c=0;c<3;++c) a[r][c]=expected[3*r+c];
    for(size_t i=0;i<6;++i) {
        a[i/3][i%3]=std::numeric_limits<T>::quiet_NaN();
        check(!sigmoid_values(a),type+" checker_rejects_nan_position_"+std::to_string(i));
        a[i/3][i%3]=expected[i];
    }
    a[0][0]=std::numeric_limits<T>::infinity();
    check(!sigmoid_values(a),type+" checker_rejects_positive_infinity");
    a[0][0]=-std::numeric_limits<T>::infinity();
    check(!sigmoid_values(a),type+" checker_rejects_negative_infinity");
    a[0][0]=T(0.25);
    check(!sigmoid_values(a),type+" checker_rejects_finite_mismatch");
    for(size_t r=0;r<2;++r) for(size_t c=0;c<3;++c) a[r][c]=std::numeric_limits<T>::quiet_NaN();
    check(!sigmoid_values(a),type+" checker_rejects_all_nan");
}
template<class T> void sigmoid(const std::string& type) {
    using M=Matrix::Matrix<T>; M a(2,3);
    for(size_t r=0;r<2;++r) for(size_t c=0;c<3;++c) a[r][c]=T(int(3*r+c)-2);
#ifdef NEURO
    static_assert(std::is_same<decltype(std::declval<const M&>().SigmoidMatrix()),M>::value,"Neuro const value return");
    const M& input=a; auto result=input.SigmoidMatrix();
    bool unchanged=true;
    for(size_t r=0;r<2;++r) for(size_t c=0;c<3;++c) unchanged=unchanged && a[r][c]==T(int(3*r+c)-2);
    check(unchanged,type+" sigmoid_input_unchanged");
    check(sigmoid_values(result),type+" sigmoid_values");
    result[0][0]=T(42); check(a[0][0]==T(-2),type+" sigmoid_independent_result");
    const M e; auto result_empty=e.SigmoidMatrix();
    check(result_empty.rows()==0 && result_empty.cols()==0,type+" sigmoid_empty_copy");
#else
    static_assert(std::is_same<decltype(std::declval<M&>().SigmoidMatrix()),M&>::value,"Bitmap reference return");
    auto& result=a.SigmoidMatrix();
    check(&result==&a,type+" sigmoid_returns_input_reference");
    check(sigmoid_values(a),type+" sigmoid_values_mutate_input");
    result[0][0]=T(42); check(a[0][0]==T(42),type+" sigmoid_result_aliases_input");
    M e; auto& result_empty=e.SigmoidMatrix();
    check(&result_empty==&e && e.rows()==0 && e.cols()==0,type+" sigmoid_empty_reference");
#endif
}
int main() {
    splits<float>("float"); splits<double>("double"); splits<int>("int");
    sigmoid<float>("float"); sigmoid<double>("double");
    sigmoid_checker_regressions<float>("float"); sigmoid_checker_regressions<double>("double");
    std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n'; return failed?1:0;
}
