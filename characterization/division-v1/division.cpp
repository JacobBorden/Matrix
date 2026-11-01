// One immutable legacy definition per executable; checks remain active under NDEBUG.
#include LEGACY_HEADER
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
int passed=0,failed=0;
void check(bool ok,const std::string& name){std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n';ok?++passed:++failed;}
template<class F> bool rejects(F f){try{f();}catch(const std::runtime_error&){return true;}return false;}
template<class T> Matrix::Matrix<T> input(){Matrix::Matrix<T> m(2,3);const int a[6]={-9,-4,0,3,8,13};for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)m[r][c]=T(a[3*r+c]);return m;}
template<class T> bool values(const Matrix::Matrix<T>& m,T divisor){const int a[6]={-9,-4,0,3,8,13};for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c){T expected=T(a[3*r+c])/divisor;if(m[r][c]!=expected)return false;}return m.rows()==2&&m.cols()==3;}
template<class T> void basic(const std::string& tag){
    using M=Matrix::Matrix<T>;
    static_assert(std::is_same<decltype(std::declval<const M&>()/T(2)),M>::value,"division returns value");
    static_assert(std::is_same<decltype(std::declval<M&>()/=T(2)),M&>::value,"assignment returns reference");
    M source=input<T>();const M& constant=source;auto result=constant/T(2);
    check(values(result,T(2)),tag+" copy_values_shape");check(values(source,T(1)),tag+" source_unchanged");
    result[0][0]=T(99);check(source[0][0]==T(-9),tag+" copy_independent");
    M assigned=input<T>();M& alias=(assigned/=T(2));check(&alias==&assigned,tag+" assignment_alias");check(values(assigned,T(2)),tag+" assignment_values");
    check(values(constant/T(-2),T(-2)),tag+" negative_divisor");
    M empty;auto e=empty/T(2);check(e.rows()==0&&e.cols()==0,tag+" empty_copy");
    check(&(empty/=T(2))==&empty&&empty.rows()==0&&empty.cols()==0,tag+" empty_assignment");
    if constexpr(std::is_floating_point<T>::value){
        for(T divisor:{T(0),T(-0.0),T(1e-12),T(-1e-12)}){
            M a=input<T>();check(rejects([&]{auto x=a/divisor;(void)x;}),tag+" small_copy_rejected");check(values(a,T(1)),tag+" small_copy_preserved");
            check(rejects([&]{a/=divisor;}),tag+" small_assignment_rejected");check(values(a,T(1)),tag+" small_assignment_preserved");
        }
        const T threshold=T(1e-9);
        const T cases[3]={std::nextafter(threshold,T(0)),threshold,std::nextafter(threshold,std::numeric_limits<T>::infinity())};
        for(T d:cases){
#ifdef NEURO
            bool expected_reject=d<threshold;
#else
            // Bitmap compares promoted float to double 1e-9; rounded float threshold lies below it.
            bool expected_reject=static_cast<double>(d)<1e-9;
#endif
            M a=input<T>();bool copy_reject=false;bool ok=false;try{auto x=a/d;ok=values(x,d);}catch(const std::runtime_error&){copy_reject=true;}
            check(copy_reject==expected_reject&&(copy_reject||ok)&&values(a,T(1)),tag+" threshold_copy");
            bool assignment_reject=false;try{a/=d;}catch(const std::runtime_error&){assignment_reject=true;}
            check(assignment_reject==expected_reject&&values(a,expected_reject?T(1):d),tag+" threshold_assignment");
        }
    }
}
int main(int argc,char** argv){
    if(argc==2){ // Diagnostic-only child. UBSan terminates Neuro before undefined division completes.
        auto m=input<int>();volatile int runtime_zero=0;const int zero=runtime_zero;
        bool rejected=rejects([&]{if(std::string(argv[1])=="--integer-zero-assign")m/=zero;else{auto x=m/zero;(void)x;}});
        std::cout<<(rejected?"INTEGER_ZERO_REJECTED":"INTEGER_ZERO_NOT_REJECTED")<<'\n';return rejected?0:2;
    }
    basic<float>("float");basic<double>("double");basic<int>("int");
    Matrix::Matrix<int> empty;
#ifdef NEURO
    check(!rejects([&]{auto x=empty/0;check(x.rows()==0&&x.cols()==0,"int empty_zero_shape");}),"int empty_zero_copy_bypasses_guard");
    check(!rejects([&]{empty/=0;}),"int empty_zero_assignment_bypasses_guard");
#else
    check(rejects([&]{auto x=empty/0;(void)x;}),"int empty_zero_copy_rejected");
    check(rejects([&]{empty/=0;}),"int empty_zero_assignment_rejected");
    check(empty.rows()==0&&empty.cols()==0,"int empty_zero_shape");
#endif
    std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;
}
