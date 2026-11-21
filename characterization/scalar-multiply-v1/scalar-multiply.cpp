#include LEGACY_HEADER
#include <iostream>
#include <type_traits>
#include <utility>
int passed=0,failed=0;
void check(bool ok,const char* name){std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n';ok?++passed:++failed;}
template<class T> bool values(const Matrix::Matrix<T>& m,const int* golden){
 if(m.rows()!=2||m.cols()!=3||m.size()!=6)return false;
 for(size_t i=0;i<2;++i)for(size_t j=0;j<3;++j)if(m[i][j]!=T(golden[i*3+j]))return false;
 return true;
}
template<class T> void cases(){
 using M=Matrix::Matrix<T>;
 static_assert(std::is_same<decltype(std::declval<const M&>()*T(2)),M>::value,"const multiplication returns value");
 static_assert(std::is_same<decltype(std::declval<M&>()*=T(2)),M&>::value,"assignment returns reference");
 int original[]={-4,2,7,0,-3,9};
 int golden[][6]={{0,0,0,0,0,0},{-4,2,7,0,-3,9},{8,-4,-14,0,6,-18}};
 T factors[]={T(0),T(1),T(-2)};
 M source(2,3);for(size_t i=0;i<2;++i)for(size_t j=0;j<3;++j)source[i][j]=T(original[i*3+j]);
 const M& input=source;
 for(int n=0;n<3;++n){
  auto result=input*factors[n];check(values(result,golden[n]),"scalar golden values");
  check(values(source,original),"const source unchanged");
  result[0][0]=T(99);check(values(source,original),"binary storage independent");
  M compound=source;check(&(compound*=factors[n])==&compound&&values(compound,golden[n]),"compound identity and golden values");
 }
 M alias=source;alias*=alias[0][0];
 int stable[]={16,-8,-28,0,12,-36};int changing[]={16,32,112,0,-48,144};
#ifdef NEURO
 check(values(alias,stable),"aliased scalar captured by value");
 auto left=T(-2)*input;check(values(left,golden[2]),"left scalar overload golden");
#else
 (void)stable;check(values(alias,changing),"aliased scalar changes after first element");
#endif
 for(auto dims:{std::pair<size_t,size_t>{0,0},{0,3},{2,0}}){
  M empty(dims.first,dims.second);auto result=empty*T(-2);
  check(result.rows()==dims.first&&result.cols()==dims.second&&result.size()==0,"empty binary shape retained");
  check(&(empty*=T(-2))==&empty&&empty.rows()==dims.first&&empty.cols()==dims.second&&empty.size()==0,"empty compound shape and identity");
 }
 if constexpr(std::is_floating_point<T>::value){
  M fractional(1,3);fractional[0][0]=T(-1.5);fractional[0][1]=T(.25);fractional[0][2]=T(8);
  auto result=fractional*T(.5);
  check(result[0][0]==T(-.75)&&result[0][1]==T(.125)&&result[0][2]==T(4),"exact binary fractions");
  fractional*=T(.5);check(fractional[0][0]==T(-.75)&&fractional[0][1]==T(.125)&&fractional[0][2]==T(4),"exact compound fractions");
 }
}
int main(){cases<int>();cases<float>();cases<double>();std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;}
