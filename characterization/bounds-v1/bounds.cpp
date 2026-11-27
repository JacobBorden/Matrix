// Each executable includes exactly one immutable legacy definition.
#include LEGACY_HEADER
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
int passed=0,failed=0;
void check(bool ok,const std::string& name) {
    std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n'; ok?++passed:++failed;
}
template<class F> bool rejects(F f) {
    try { f(); } catch(const std::out_of_range&) { return true; } catch(...) { return false; }
    return false;
}
template<class M,class=void> struct has_at:std::false_type {};
template<class M> struct has_at<M,std::void_t<decltype(std::declval<M&>().at(0,0))>>:std::true_type {};
template<class T> void cases(const std::string& name) {
    using M=Matrix::Matrix<T>;
    M m(2,3); for(size_t r=0;r<2;++r) for(size_t c=0;c<3;++c) m[r][c]=T(10*r+c+1);
    const M& cm=m; auto& row=m[1]; const auto& cr=cm[1];
    auto test=[&](bool ok,const char* label){check(ok,name+" "+label);};
    test((std::is_same_v<decltype(m[1][2]),T&>),"mutable_subscript_reference");
    test((std::is_same_v<decltype(cm[1][2]),const T&>),"const_subscript_reference");
    test(cm[0][0]==T(1)&&cm[1][2]==T(13),"asymmetric_valid_values");
#ifdef NEURO
    test((std::is_same_v<decltype(row.at(2)),T>),"mutable_row_at_value");
#else
    test((std::is_same_v<decltype(row.at(2)),T&>),"mutable_row_at_reference");
#endif
    test((std::is_same_v<decltype(cr.at(2)),T>),"const_row_at_value");
    auto copied=cr.at(2);copied=T(90);
    test(cr.at(2)==T(13)&&copied==T(90),"const_row_at_independent_copy");
    decltype(auto) access=row.at(2);access=T(77);
#ifdef NEURO
    test(row[2]==T(13)&&access==T(77),"mutable_row_at_independent_copy");
    test(!has_at<M>::value,"matrix_at_absent");
#else
    test(row[2]==T(77)&&&access==&row[2],"mutable_row_at_alias");
    test(has_at<M>::value,"matrix_at_present");
#endif
    const size_t neg=static_cast<size_t>(-1), max=std::numeric_limits<size_t>::max();
    test(rejects([&]{(void)row.at(3);}),"row_at_mutable_endpoint");
    test(rejects([&]{(void)cr.at(3);}),"row_at_const_endpoint");
    test(rejects([&]{(void)row.at(neg);}),"row_at_mutable_negative_conversion");
    test(rejects([&]{(void)cr.at(neg);}),"row_at_const_negative_conversion");
    test(rejects([&]{(void)row.at(max);}),"row_at_mutable_size_max");
    test(rejects([&]{(void)cr.at(max);}),"row_at_const_size_max");
    Matrix::MatrixRow<T> emptyrow; const auto& cer=emptyrow;
    test(rejects([&]{(void)emptyrow.at(0);}),"empty_row_at_mutable");
    test(rejects([&]{(void)cer.at(0);}),"empty_row_at_const");
    M empty; const M& ce=empty;
#ifdef NEURO
    test(rejects([&]{(void)m[2];}),"matrix_subscript_mutable_endpoint");
    test(rejects([&]{(void)cm[2];}),"matrix_subscript_const_endpoint");
    test(rejects([&]{(void)m[neg];}),"matrix_subscript_mutable_negative_conversion");
    test(rejects([&]{(void)cm[neg];}),"matrix_subscript_const_negative_conversion");
    test(rejects([&]{(void)row[3];}),"row_subscript_mutable_endpoint");
    test(rejects([&]{(void)cr[3];}),"row_subscript_const_endpoint");
    test(rejects([&]{(void)empty[0];}),"empty_matrix_subscript_mutable");
    test(rejects([&]{(void)ce[0];}),"empty_matrix_subscript_const");
#else
    test(rejects([&]{(void)m.at(2,0);}),"matrix_at_mutable_row_endpoint");
    test(rejects([&]{(void)cm.at(2,0);}),"matrix_at_const_row_endpoint");
    test(rejects([&]{(void)m.at(0,3);}),"matrix_at_mutable_column_endpoint");
    test(rejects([&]{(void)cm.at(0,3);}),"matrix_at_const_column_endpoint");
    test(rejects([&]{(void)m.at(neg,0);}),"matrix_at_mutable_negative_row");
    test(rejects([&]{(void)cm.at(neg,0);}),"matrix_at_const_negative_row");
    test(rejects([&]{(void)empty.at(0,0);}),"empty_matrix_at_mutable");
    test(rejects([&]{(void)ce.at(0,0);}),"empty_matrix_at_const");
#endif
}
int main(){cases<float>("float");cases<double>("double");cases<int>("int");
std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;}
