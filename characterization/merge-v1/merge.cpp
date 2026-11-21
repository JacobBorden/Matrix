// Legacy definitions intentionally compile into separate executables.
#include LEGACY_HEADER
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
int passed=0,failed=0;
void check(bool ok,const std::string& label) {
    std::cout<<(ok?"PASS ":"FAIL ")<<label<<'\n';
    ok?++passed:++failed;
}
template<class T> bool equal(const Matrix::Matrix<T>& m,size_t r,size_t c,std::initializer_list<int> values) {
    if(m.rows()!=r || m.cols()!=c || m.size()!=r*c || values.size()!=r*c)return false;
    auto it=values.begin();
    for(size_t i=0;i<r;++i)for(size_t j=0;j<c;++j)if(m[i][j]!=*it++)return false;
    return true;
}
template<class T> void cases(const std::string& label) {
    using M=Matrix::Matrix<T>;
    static_assert(std::is_same<decltype(std::declval<const M&>().MergeVertical(std::declval<const M&>())),M>::value,"vertical value");
    static_assert(std::is_same<decltype(std::declval<const M&>().MergeHorizontal(std::declval<const M&>())),M>::value,"horizontal value");
    M a(2,2),b(1,2),c(2,1);
    a[0][0]=-3;a[0][1]=2;a[1][0]=7;a[1][1]=-1;
    b[0][0]=11;b[0][1]=-9;c[0][0]=4;c[1][0]=8;
    const M& ca=a;const M& cb=b;const M& cc=c;
    auto v=ca.MergeVertical(cb),h=ca.MergeHorizontal(cc);
    check(equal(v,3,2,{-3,2,7,-1,11,-9}),label+" vertical golden");
    check(equal(h,2,3,{-3,2,4,7,-1,8}),label+" horizontal golden");
    check(equal(cb.MergeVertical(ca),3,2,{11,-9,-3,2,7,-1}),label+" reversed vertical");
    check(equal(cc.MergeHorizontal(ca),2,3,{4,-3,2,8,7,-1}),label+" reversed horizontal");
    check(equal(ca,2,2,{-3,2,7,-1}) && equal(cb,1,2,{11,-9}) && equal(cc,2,1,{4,8}),label+" inputs unchanged");
    check(equal(ca.MergeVertical(ca),4,2,{-3,2,7,-1,-3,2,7,-1}),label+" self vertical");
    check(equal(ca.MergeHorizontal(ca),2,4,{-3,2,-3,2,7,-1,7,-1}),label+" self horizontal");
    v[0][0]=99;h[0][0]=88;
    check(a[0][0]==-3,label+" result mutations independent");
    a[1][1]=66;b[0][1]=55;c[1][0]=44;
    check(v[1][1]==-1 && v[2][1]==-9 && h[1][2]==8,label+" operand mutations independent");
    M badv(1,3),badh(3,1);
    bool rejected=false;try{ca.MergeVertical(badv);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,label+" vertical mismatch rejects");
    rejected=false;try{ca.MergeHorizontal(badh);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,label+" horizontal mismatch rejects");
    check(equal(ca,2,2,{-3,2,7,66}),label+" mismatch preserves source");
    const M empty;
    check(equal(empty.MergeVertical(empty),0,0,{}),label+" empty vertical");
    check(equal(empty.MergeHorizontal(empty),0,0,{}),label+" empty horizontal");
    check(equal(ca.MergeVertical(empty),2,2,{-3,2,7,66}),label+" right empty vertical");
    check(equal(empty.MergeVertical(ca),2,2,{-3,2,7,66}),label+" left empty vertical");
    check(equal(ca.MergeHorizontal(empty),2,2,{-3,2,7,66}),label+" right empty horizontal");
    check(equal(empty.MergeHorizontal(ca),2,2,{-3,2,7,66}),label+" left empty horizontal");
    auto copy=ca.MergeVertical(empty);copy[0][0]=77;
    check(a[0][0]==-3,label+" empty identity deep copy");
    const M zeroRows(0,2),zeroCols(2,0);
    check(equal(zeroRows.MergeVertical(ca),2,2,{-3,2,7,66}),label+" zero rows vertical left");
    check(equal(ca.MergeVertical(zeroRows),2,2,{-3,2,7,66}),label+" zero rows vertical right");
    check(equal(zeroCols.MergeHorizontal(ca),2,2,{-3,2,7,66}),label+" zero cols horizontal left");
    check(equal(ca.MergeHorizontal(zeroCols),2,2,{-3,2,7,66}),label+" zero cols horizontal right");
}
int main(){cases<int>("int");cases<float>("float");cases<double>("double");std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;}
