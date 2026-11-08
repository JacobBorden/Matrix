// Compile each pinned legacy definition separately.
#include LEGACY_HEADER
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
int passed=0,failed=0;
void check(bool ok,const std::string& s){std::cout<<(ok?"PASS ":"FAIL ")<<s<<'\n';ok?++passed:++failed;}
template<class T> void cases(){
 using M=Matrix::Matrix<T>;
 M m(2,3); int golden[]={-3,2,7,11,-9,4};
 for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)m[r][c]=T(golden[r*3+c]);
 auto b=m.begin(),e=m.end();
 check(b!=e,"nonempty range");
 auto old=b++;check((*old)[0]==-3 && (*b)[0]==11,"row postfix increment");
 auto previous=b--;check((*previous)[0]==11 && b==m.begin(),"row postfix decrement");
 check(&(++b)==&b && (*b)[2]==4,"row prefix increment identity");
 check(&(--b)==&b && (*b)[2]==7,"row prefix decrement identity");
 check(b.operator->()==&m[0],"row arrow aliases storage");
 size_t r=0,n=0;bool values=true;
 for(auto row=m.begin();row!=m.end();++row,++r){
  size_t c=0;for(auto it=row->begin();it!=row->end();++it,++c,++n)values=values && *it==golden[r*3+c];
  values=values && c==3;
 }
 check(values && r==2 && n==6,"asymmetric row-major golden traversal");
 // Arithmetic stays within one allocated row, including its one-past endpoint.
 auto first=m[0].begin(),last=m[0].end();
 check(last-first==3 && first-last==-3,"signed row distances");
 auto middle=first+1;
 check(*middle==2 && *first==-3,"addition preserves source");
 check(*(last-1)==7,"subtraction from endpoint");
 check(*(middle+(-1))==-3 && *(middle-(-1))==7,"negative offsets");
 check(&(middle+=1)==&middle && *middle==7,"plus assign identity");
 check(&(middle-=1)==&middle && *middle==2,"minus assign identity");
 check(&(middle+=(-1))==&middle && middle==first,"negative plus assign");
 check(&(middle-=(-1))==&middle && *middle==2,"negative minus assign");
 check(first< middle && middle<last && last>middle && middle>first,"strict ordering");
 check(first<=first && last>=last && first<=last && last>=first,"inclusive ordering");
 check(!(first<first) && !(last>last) && !(last<=first) && !(first>=last),"ordering negatives");
 const auto fixed=middle;
 check(fixed[-1]==-3 && fixed[0]==2 && fixed[1]==7,"const iterator signed indexing");
 check(&fixed[1]==&m[0][2],"index aliases row storage");
 fixed[1]=T(19);check(m[0][2]==19,"indexed mutation observed");
 m[0][2]=T(7);
 check(std::distance(first,last)==3,"standard distance within row");
 auto it=m[0].begin();auto ri=it++;
 check(*ri==-3 && *it==2,"element postfix increment");
 auto rd=it--;check(*rd==2 && *it==-3,"element postfix decrement");
 check(&(++it)==&it && *it==2,"element prefix identity");
 check(&(--it)==&it && *it==-3,"element decrement identity");
 *it=T(33);check(m[0][0]==33,"element mutable alias");
 (*m.begin())[1]=T(44);check(m[0][1]==44,"row mutable alias");
 const auto& row=m[1];auto ci=row.begin();
 static_assert(std::is_same<decltype(*ci),T&>::value,"legacy const row exposes mutable reference");
 *ci=T(55);check(m[1][0]==55,"legacy const row mutation observed");
#ifdef NEURO
 // This source has no const Matrix begin/end; do not cast away const.
#else
 const M& cm=m;auto cmi=cm.begin();(*cmi)[2]=T(66);
 check(m[0][2]==66,"legacy const matrix mutation observed");
 size_t count=0;for(auto x=cm.begin();x!=cm.end();++x)++count;
 check(count==2,"const matrix traverses rows");
#endif
 M empty;check(empty.begin()==empty.end(),"default empty range equality");
}
int main(){cases<int>();cases<float>();cases<double>();std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;}
