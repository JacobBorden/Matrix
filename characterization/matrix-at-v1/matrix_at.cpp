// Bitmap only: never combine conflicting immutable legacy definitions.
#include LEGACY_HEADER
#include <iostream>
#include <string>
#include <type_traits>
int passed=0,failed=0;
void check(bool ok,const std::string& label){std::cout<<(ok?"PASS ":"FAIL ")<<label<<'\n';ok?++passed:++failed;}
template<class T> void cases(const std::string& name){
    using M=Matrix::Matrix<T>;
    M m(2,3); const M& cm=m;
    auto test=[&](bool ok,const char* label){check(ok,name+" "+label);};
    const int values[2][3]={{2,-7,19},{31,43,-59}};
    for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)m[r][c]=T(values[r][c]);
    test((std::is_same_v<decltype(m.at(0,0)),T&>),"mutable_return_reference");
    test((std::is_same_v<decltype(cm.at(0,0)),T>),"const_return_value");
    bool all_mutable=true,all_const=true,all_addresses=true;
    for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c){
        all_mutable=all_mutable&&m.at(r,c)==T(values[r][c]);
        all_const=all_const&&cm.at(r,c)==T(values[r][c]);
        all_addresses=all_addresses&&&m.at(r,c)==&m[r][c];
    }
    test(all_mutable,"all_asymmetric_mutable_values");
    test(all_const,"all_asymmetric_const_values");
    test(all_addresses,"all_mutable_addresses_match_storage");
    auto snapshot=cm.at(1,2); snapshot=T(71);
    test(snapshot==T(71)&&m[1][2]==T(-59),"const_copy_mutation_independent");
    decltype(auto) alias=m.at(1,2);alias=T(83);
    test(&alias==&m[1][2]&&m[1][2]==T(83),"mutable_alias_writes_storage");
    test(cm.at(1,2)==T(83),"const_read_observes_alias_write");
    m[1][2]=T(-97);
    test(alias==T(-97),"alias_observes_subscript_write");
    test(snapshot==T(71),"const_snapshot_survives_storage_write");
    m.at(0,1)=T(101);
    test(cm[0][1]==T(101),"direct_at_assignment_visible");
    bool neighbors=true;
    for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)
        if(!(r==1&&c==2)&&!(r==0&&c==1))neighbors=neighbors&&cm.at(r,c)==T(values[r][c]);
    test(neighbors,"unrelated_cells_preserved");
    test(m.rows()==2&&m.cols()==3&&m.size()==6,"access_preserves_shape");
}
int main(){cases<float>("float");cases<double>("double");cases<int>("int");
std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;}
