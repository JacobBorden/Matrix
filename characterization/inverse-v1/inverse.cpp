// Exactly one immutable legacy definition per executable.
#include LEGACY_HEADER
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
int passed=0,failed=0;
void check(bool ok,const std::string& name){std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n';ok?++passed:++failed;}
template<class E,class F> bool rejects(F f){try{f();}catch(const E&){return true;}return false;}
template<class T> bool close(T a,T b){
    if constexpr(std::is_floating_point<T>::value)
        return std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=T(64)*std::numeric_limits<T>::epsilon()*std::max(T(1),std::abs(b));
    else return a==b;
}
template<class T> bool values(const Matrix::Matrix<T>& m,const T* expected,size_t rows,size_t cols){
    if(m.rows()!=rows||m.cols()!=cols)return false;
    for(size_t r=0;r<rows;++r)for(size_t c=0;c<cols;++c)if(!close(m[r][c],expected[r*cols+c]))return false;
    return true;
}
template<class T> void fixtures(const std::string& type){
    using M=Matrix::Matrix<T>;
    M empty;auto ei=empty.Inverse();check(ei.rows()==0&&ei.cols()==0,type+" empty");
    M scalar(1,1);scalar[0][0]=T(-4);const T scalar_expected[1]={T(1)/T(-4)};
    check(values(scalar.Inverse(),scalar_expected,1,1),type+" scalar");
    M a(2,2);const T original[4]={4,7,2,6};
    for(size_t r=0;r<2;++r)for(size_t c=0;c<2;++c)a[r][c]=original[2*r+c];
    const M& ca=a;auto inverse=ca.Inverse();
    // Analytic adjugate / 10; legacy integer reciprocal truncates to zero.
    const T expected[4]={T(6)*(T(1)/T(10)),T(-7)*(T(1)/T(10)),T(-2)*(T(1)/T(10)),T(4)*(T(1)/T(10))};
    check(values(inverse,expected,2,2),type+" two_by_two_analytic");
    check(values(a,original,2,2),type+" const_input_unchanged");
    inverse[0][0]=T(99);check(values(a,original,2,2),type+" independent_result");
    M three(3,3);const T v[9]={1,2,3,0,1,4,5,6,0};const T inv3[9]={-24,18,5,20,-15,-4,-5,4,1};
    for(size_t r=0;r<3;++r)for(size_t c=0;c<3;++c)three[r][c]=v[3*r+c];
    check(values(three.Inverse(),inv3,3,3),type+" three_by_three_unimodular");
    M swap(2,2);swap[0][1]=1;swap[1][0]=1;const T permutation[4]={0,1,1,0};
    check(values(swap.Inverse(),permutation,2,2),type+" negative_determinant");
    M wide(2,3),tall(3,2);wide[0][0]=T(9);
    check(rejects<std::invalid_argument>([&]{wide.Inverse();}),type+" wide_rejected");
    check(rejects<std::invalid_argument>([&]{tall.Inverse();}),type+" tall_rejected");
    check(wide.rows()==2&&wide.cols()==3&&wide[0][0]==T(9),type+" nonsquare_input_preserved");
    if constexpr(std::is_floating_point<T>::value){
        M singular(2,2);singular[0][0]=1;singular[0][1]=2;singular[1][0]=2;singular[1][1]=4;
        check(rejects<std::runtime_error>([&]{singular.Inverse();}),type+" singular_rejected");
        const T s[4]={1,2,2,4};check(values(singular,s,2,2),type+" singular_preserved");
        M tiny(1,1);tiny[0][0]=T(1e-12);check(rejects<std::runtime_error>([&]{tiny.Inverse();}),type+" small_nonsingular_rejected");
        const T threshold=T(1e-9);
        tiny[0][0]=std::nextafter(threshold,T(0));check(rejects<std::runtime_error>([&]{tiny.Inverse();}),type+" below_threshold_rejected");
        tiny[0][0]=threshold;
#ifdef NEURO
        const bool exact_reject=false;
#else
        const bool exact_reject=double(threshold)<1e-9;
#endif
        bool rejected=false,correct=true;try{auto result=tiny.Inverse();correct=close(result[0][0],T(1)/threshold);}catch(const std::runtime_error&){rejected=true;}
        check(rejected==exact_reject&&correct,type+" exact_threshold_contract");
        tiny[0][0]=std::nextafter(threshold,std::numeric_limits<T>::infinity());
        auto above=tiny.Inverse();check(close(above[0][0],T(1)/tiny[0][0]),type+" above_threshold_accepted");
        tiny[0][0]=-threshold;check(rejects<std::runtime_error>([&]{tiny.Inverse();})==exact_reject,type+" negative_threshold_contract");
    }
}
int main(int argc,char**){
    if(argc>1){Matrix::Matrix<int> singular(2,2);singular[0][0]=1;singular[0][1]=2;singular[1][0]=2;singular[1][1]=4;
        try{auto result=singular.Inverse();std::cout<<"UNEXPECTED_RETURN "<<result[0][0]<<'\n';return 2;}
        catch(const std::runtime_error&){std::cout<<"INTEGER_SINGULAR_REJECTED\n";return 0;}}
    fixtures<float>("float");fixtures<double>("double");fixtures<int>("int");
    check(!close(std::numeric_limits<double>::quiet_NaN(),1.0),"checker_nan");
    check(!close(std::numeric_limits<double>::infinity(),1.0),"checker_positive_infinity");
    check(!close(-std::numeric_limits<double>::infinity(),1.0),"checker_negative_infinity");
    check(!close(2.0,1.0),"checker_finite_mismatch");check(close(0.0,0.0),"checker_zero");
    std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;
}
