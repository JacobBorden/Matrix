// Compile each pinned legacy definition separately.
#include LEGACY_HEADER
#include <algorithm>
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
#ifdef NEURO
 // Row objects are allocated in one array; scalar rows remain separate allocations.
 auto rb=m.begin(),re=m.end(),second=rb+1;
 check(re-rb==2 && rb-re==-2,"matrix signed row distances");
 check((*second)[0]==11 && rb==m.begin(),"row addition preserves source");
 auto tail=re-1;check(tail==second && (*tail)[2]==4,"row subtraction from endpoint");
 check(second+(-1)==rb && rb-(-1)==second,"row negative offsets");
 check(&(second+=(-1))==&second && second==rb,"row negative plus assign identity");
 check(&(second-=(-1))==&second && (*second)[1]==-9,"row negative minus assign identity");
 check(&(second+=1)==&second && second==re,"row plus assign endpoint");
 check(&(second-=1)==&second && second==tail,"row minus assign endpoint");
 check(rb<second && second<re && re>second && second>rb,"row strict ordering");
 check(rb<=rb && re>=re && rb<=re && re>=rb,"row inclusive ordering");
 check(!(rb<rb) && !(re>re) && !(re<=rb) && !(rb>=re),"row ordering negatives");
 check(std::distance(rb,re)==2,"standard distance over row objects");
 (*tail)[0]=T(71);check(m[1][0]==71 && &*tail==&m[1],"row arithmetic storage alias");
 m[1][0]=T(11);
#endif

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
template<class T> void narrow_shapes(){
 using M=Matrix::Matrix<T>;
 for(auto shape : {std::pair<size_t,size_t>{1,1},{1,3},{3,1}}){
  const auto rows=shape.first,cols=shape.second;
  M m(rows,cols);
  for(size_t r=0;r<rows;++r)for(size_t c=0;c<cols;++c)m[r][c]=T(10*r+2*c);
  size_t r=0,n=0;bool ok=true;
  for(auto row=m.begin();row!=m.end();++row,++r){
   size_t c=0;
   for(auto it=row->begin();it!=row->end();++it,++c,++n)
    ok=ok && *it==T(10*r+2*c);
   ok=ok && c==cols;
  }
  check(ok && r==rows && n==rows*cols,"narrow shape forward golden traversal");
  auto row=m.end();r=rows;n=0;ok=true;
  while(row!=m.begin()){
   --row;--r;auto it=row->end();size_t c=cols;
   while(it!=row->begin()){--it;--c;++n;ok=ok && *it==T(10*r+2*c);}
   ok=ok && c==0;
  }
  check(ok && r==0 && n==rows*cols,"narrow shape reverse golden traversal");
  auto tail=m.end();--tail;auto cell=tail->end();--cell;
  check(&*tail==&m[rows-1] && &*cell==&m[rows-1][cols-1],"narrow shape endpoint aliases");
  *cell=T(-17);check(m[rows-1][cols-1]==T(-17),"narrow shape endpoint mutation");
  ++cell;++tail;
  check(cell==m[rows-1].end() && tail==m.end(),"narrow shape endpoint restoration");
 }
}
template<class T> void row_algorithms(){
 using M=Matrix::Matrix<T>;
 M m(2,3);
 const T input[]={T(-7),T(2),T(11)};
 for(size_t c=0;c<3;++c){m[0][c]=input[c];m[1][c]=T(40+c);}
 auto first=m[0].begin(),last=m[0].end();
 T* address=&m[0][0];
 std::reverse(first,last);
 check(m[0][0]==T(11) && m[0][1]==T(2) && m[0][2]==T(-7),"standard reverse golden values");
 check(&*first==address && &m[0][0]==address,"standard reverse preserves storage aliases");
 check(m[1][0]==T(40) && m[1][1]==T(41) && m[1][2]==T(42),"standard reverse preserves other row");
 std::reverse(first,last);
 check(m[0][0]==input[0] && m[0][1]==input[1] && m[0][2]==input[2],"standard reverse twice restores values");
 T copied[3]={};
 auto output=std::copy(first,last,copied);
 check(output==copied+3 && copied[0]==input[0] && copied[1]==input[1] && copied[2]==input[2],"standard copy to array golden values and endpoint");
 const T replacement[]={T(5),T(-13),T(8)};
 auto end=std::copy(replacement,replacement+3,first);
 check(end==last && m[0][0]==T(5) && m[0][1]==T(-13) && m[0][2]==T(8),"standard copy to row golden values and endpoint");
 check(&*first==address && m[1][0]==T(40) && m[1][2]==T(42),"standard copy preserves storage and other row");
}
int main(){row_algorithms<int>();row_algorithms<float>();row_algorithms<double>();cases<int>();cases<float>();cases<double>();narrow_shapes<int>();narrow_shapes<float>();narrow_shapes<double>();std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;}
