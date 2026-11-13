#include LEGACY_HEADER
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
int passed=0,failed=0;
void check(bool b,const char* label){std::cout<<(b?"PASS ":"FAIL ")<<label<<'\n';b?++passed:++failed;}
template<class T> bool values(const Matrix::Matrix<T>& m, int which){
 if(m.rows()!=2||m.cols()!=3)return false;
 const int golden[][6]={{-4,2,7,0,-3,9},{1,-2,3,4,5,-6},{-3,0,10,4,2,3},{-5,4,4,-4,-8,15},{-2,4,9,2,-1,11},{-6,0,5,-2,-5,7},{-8,4,14,0,-6,18},{0,0,0,0,0,0},{-1,2,-3,-4,-5,6}};
 for(size_t i=0;i<2;++i)for(size_t j=0;j<3;++j)if(m[i][j]!=T(golden[which][i*3+j]))return false;
 return true;
}
template<class T> void cases(){
 using M=Matrix::Matrix<T>;
 static_assert(std::is_same<decltype(std::declval<const M&>()+std::declval<const M&>()),M>::value,"value addition");
 static_assert(std::is_same<decltype(std::declval<M&>()-=std::declval<const M&>()),M&>::value,"compound reference");
 M a(2,3),b(2,3);int av[]={-4,2,7,0,-3,9},bv[]={1,-2,3,4,5,-6};
 for(size_t i=0;i<2;++i)for(size_t j=0;j<3;++j){a[i][j]=T(av[i*3+j]);b[i][j]=T(bv[i*3+j]);}
 const M& ca=a;const M& cb=b;
 auto sum=ca+cb;auto diff=ca-cb;
 check(values(sum,2),"matrix sum golden");check(values(diff,3),"matrix difference golden");
 check(values(ca+T(2),4),"scalar sum golden");check(values(ca-T(2),5),"scalar difference golden");
 check(values(a,0)&&values(b,1),"const operands unchanged");
 sum[0][0]=T(99);check(values(a,0)&&values(b,1),"result storage independent");
 M c=a;check(&(c+=b)==&c&&values(c,2),"compound addition alias and values");
 c=a;check(&(c-=b)==&c&&values(c,3),"compound subtraction alias and values");
 c=a;check(&(c+=T(2))==&c&&values(c,4),"scalar compound addition");
 c=a;check(&(c-=T(2))==&c&&values(c,5),"scalar compound subtraction");
 c=a;c+=c;check(values(c,6),"self addition");c=a;c-=c;check(values(c,7),"self subtraction");
 for(int op=0;op<4;++op){c=a;M wrong(3,2);bool rejected=false;try{if(op==0)(void)(c+wrong);if(op==1)(void)(c-wrong);if(op==2)c+=wrong;if(op==3)c-=wrong;}catch(const std::invalid_argument&){rejected=true;}check(rejected&&values(c,0),"shape error retains operand");}
 M empty;
 for(int op=0;op<4;++op){bool rejected=false;M result;try{if(op==0)result=empty+b;if(op==1)result=empty-b;if(op==2)result=b+empty;if(op==3)result=b-empty;}catch(const std::invalid_argument&){rejected=true;}
#ifdef NEURO
 check(rejected,"empty mismatch rejected");
#else
 check(!rejected&&values(result,op==1?8:1),"empty mismatch identity or negation");
#endif
 }
 c=a;bool rejected=false;try{c+=empty;}catch(const std::invalid_argument&){rejected=true;}check(rejected&&values(c,0),"compound empty addition rejects");
 rejected=false;try{c-=empty;}catch(const std::invalid_argument&){rejected=true;}check(rejected&&values(c,0),"compound empty subtraction rejects");
}
int main(){cases<int>();cases<float>();cases<double>();std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;}
