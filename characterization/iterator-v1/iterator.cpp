// Compile each pinned legacy definition separately.
#include LEGACY_HEADER
#include <algorithm>
#include <cstdint>
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
template<class T> void row_sort(){
 using M=Matrix::Matrix<T>;
 M m(2,5);
 const T input[]={T(8),T(-3),T(8),T(0),T(2)};
 for(size_t c=0;c<5;++c){m[0][c]=input[c];m[1][c]=T(40+c);}
 auto first=m[0].begin(),last=m[0].end();T* address=&m[0][0];
 std::sort(first,last);
 const T ascending[]={T(-3),T(0),T(2),T(8),T(8)};
 bool ok=true;for(size_t c=0;c<5;++c)ok=ok && m[0][c]==ascending[c];
 check(ok,"standard sort ascending golden values including duplicate");
 check(&*first==address && last==m[0].end(),"standard sort preserves storage and endpoints");
 std::sort(first,last,[](const T& a,const T& b){return a>b;});
 const T descending[]={T(8),T(8),T(2),T(0),T(-3)};
 ok=true;for(size_t c=0;c<5;++c)ok=ok && m[0][c]==descending[c];
 check(ok,"standard sort comparator descending golden values");
 ok=true;for(size_t c=0;c<5;++c)ok=ok && m[1][c]==T(40+c);
 check(ok && &*first==address,"standard sort preserves adjacent row and alias");
}
template<class T> void move_aliases(){
 using M=Matrix::Matrix<T>;
 M source(2,3);
 for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)source[r][c]=T(10*r+2*c-0.0);
 auto row=source.begin(),row_end=source.end();
 auto cell=source[0].begin(),cell_end=source[0].end();
 M moved(std::move(source));
 check(row==moved.begin() && row_end==moved.end(),"move construction preserves row iterators");
 check(cell==moved[0].begin() && cell_end==moved[0].end(),"move construction preserves element iterators");
 check(&*row==&moved[0] && &*cell==&moved[0][0],"move construction storage aliases");
 bool ok=true;for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)ok=ok && moved[r][c]==T(10*r+2*c);
 check(ok && source.size()==0,"move construction values and emptied source");
 M target(1,1);target[0][0]=T(99);
 check(&(target=std::move(moved))==&target,"move assignment return identity");
 check(row==target.begin() && row_end==target.end() && cell==target[0].begin() && cell_end==target[0].end(),"move assignment preserves source iterators");
 *cell=T(-17);(*row)[2]=T(23);
 check(target[0][0]==T(-17) && target[0][2]==T(23) && target[1][2]==T(14),"saved iterators mutate new owner only");
 check(moved.size()==0 && &*row==&target[0] && &*cell==&target[0][0],"move assignment storage aliases and emptied source");
}
template<class T> void copy_aliases(){
 using M=Matrix::Matrix<T>;
 M source(2,3);
 for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)source[r][c]=T(10*r+2*c);
 auto row=source.begin(),row_end=source.end();
 auto cell=source[0].begin(),cell_end=source[0].end();
 M copied(source);
 check(row==source.begin() && row_end==source.end() && cell==source[0].begin() && cell_end==source[0].end(),"copy construction preserves source iterators");
 check(&*copied.begin()!=&*row && &*copied[0].begin()!=&*cell,"copy construction independent row and element storage");
 bool ok=true;for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)ok=ok && copied[r][c]==T(10*r+2*c);
 check(ok,"copy construction golden values");
 *cell=T(-17);(*row)[2]=T(23);
 check(source[0][0]==T(-17) && source[0][2]==T(23) && copied[0][0]==T(0) && copied[0][2]==T(4),"source iterator mutations preserve copy");
 M assigned(1,1);assigned[0][0]=T(99);
 check(&(assigned=source)==&assigned,"copy assignment return identity");
 check(row==source.begin() && row_end==source.end() && cell==source[0].begin() && cell_end==source[0].end(),"copy assignment preserves source iterators");
 check(&*assigned.begin()!=&*row && &*assigned[0].begin()!=&*cell,"copy assignment independent row and element storage");
 ok=true;for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)ok=ok && assigned[r][c]==source[r][c];
 check(ok,"copy assignment golden values");
 *assigned[0].begin()=T(71);(*assigned.begin())[2]=T(73);
 check(*cell==T(-17) && (*row)[2]==T(23) && copied[0][0]==T(0) && copied[0][2]==T(4),"assigned iterator mutations preserve source and earlier copy");
}
template<class T> void self_copy_aliases(){
 using M=Matrix::Matrix<T>;
 M m(2,3);const int golden[]={-3,2,7,11,-9,4};
 for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)m[r][c]=T(golden[r*3+c]);
 auto row=m.begin(),row_end=m.end();
 auto cell=m[0].begin(),cell_end=m[0].end();
 const M& alias=m;
 check(&(m=alias)==&m,"self copy assignment return identity");
 check(row==m.begin() && row_end==m.end(),"self copy preserves row iterators and endpoint");
 check(cell==m[0].begin() && cell_end==m[0].end(),"self copy preserves element iterators and endpoint");
 check(&*row==&m[0] && &*cell==&m[0][0],"self copy preserves storage aliases");
 bool ok=true;for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)ok=ok && m[r][c]==T(golden[r*3+c]);
 check(ok && m.rows()==2 && m.cols()==3,"self copy preserves shape and golden values");
 *cell=T(71);(*row)[2]=T(73);
 check(m[0][0]==T(71) && m[0][2]==T(73) && m[1][2]==T(4),"saved self copy iterators mutate original storage");
 auto last_row=row_end;--last_row;auto last_cell=cell_end;--last_cell;
 check(&*last_row==&m[1] && &*last_cell==&m[0][2],"saved self copy endpoints decrement to valid storage");
}
template<class T> void self_move_aliases() {
 using M = Matrix::Matrix<T>;
 M matrix(2, 3);
 const int golden[2][3] = {{-3, 2, 7}, {11, -9, 4}};
 for (size_t row = 0; row < 2; ++row) {
  for (size_t col = 0; col < 3; ++col) {
   matrix[row][col] = T(golden[row][col]);
  }
 }

 const auto row_begin = matrix.begin();
 const auto row_end = matrix.end();
 const auto cell_begin = matrix[0].begin();
 const auto cell_end = matrix[0].end();
 const auto row_address = reinterpret_cast<std::uintptr_t>(&matrix[0]);
 const auto cell_address = reinterpret_cast<std::uintptr_t>(&matrix[0][0]);

 M& alias = matrix;
 check(&(matrix = std::move(alias)) == &matrix, "self move assignment return identity");
 const bool shape_preserved = matrix.rows() == 2 && matrix.cols() == 3 && matrix.size() == 6;
 check(shape_preserved, "self move preserves shape and size");
 if (!shape_preserved) return;

 bool values_preserved = true;
 for (size_t row = 0; row < 2; ++row) {
  for (size_t col = 0; col < 3; ++col) {
   values_preserved &= matrix[row][col] == T(golden[row][col]);
  }
 }
 check(values_preserved, "self move preserves golden values");

 const bool storage_preserved =
     reinterpret_cast<std::uintptr_t>(&matrix[0]) == row_address &&
     reinterpret_cast<std::uintptr_t>(&matrix[0][0]) == cell_address;
 check(storage_preserved, "self move preserves row and element storage addresses");
 if (!storage_preserved) return;

 auto row = row_begin;
 auto cell = cell_begin;
 check(row == matrix.begin() && row_end == matrix.end(), "self move preserves row iterators and endpoint");
 check(cell == matrix[0].begin() && cell_end == matrix[0].end(), "self move preserves element iterators and endpoint");
 auto last_row = row_end;
 auto last_cell = cell_end;
 --last_row;
 --last_cell;
 check(&*last_row == &matrix[1] && &*last_cell == &matrix[0][2],
       "saved self move endpoints decrement to valid storage");
 *cell = T(71);
 (*last_row)[2] = T(73);
 check(matrix[0][0] == T(71) && matrix[1][2] == T(73) && matrix[0][2] == T(7),
       "saved self move iterators mutate original storage");
}
template<class T> void fixed_width_growth() {
 using M = Matrix::Matrix<T>;
 M matrix(2, 3);
 const int golden[3][3] = {{-3, 2, 7}, {11, -9, 4}, {0, 0, 0}};
 for (size_t row = 0; row < 2; ++row) {
  for (size_t col = 0; col < 3; ++col) {
   matrix[row][col] = T(golden[row][col]);
  }
 }

 // Capture address values while the allocations are live. Use only fresh iterators after resize.
 const auto old_row = reinterpret_cast<std::uintptr_t>(&matrix[0]);
 const auto old_first = reinterpret_cast<std::uintptr_t>(&matrix[0][0]);
 const auto old_second = reinterpret_cast<std::uintptr_t>(&matrix[1][0]);

 matrix.resize(3, 3);
 const bool shape_ok = matrix.rows() == 3 && matrix.cols() == 3 && matrix.size() == 9;
 check(shape_ok, "fixed-width growth updates shape and size");
 if (!shape_ok) return;

 check(reinterpret_cast<std::uintptr_t>(&matrix[0]) != old_row,
       "fixed-width growth replaces the row array");
 const auto first_address = reinterpret_cast<std::uintptr_t>(&matrix[0][0]);
 const auto second_address = reinterpret_cast<std::uintptr_t>(&matrix[1][0]);
#ifdef NEURO
 check(first_address == old_first && second_address == old_second,
       "Neuro fixed-width growth retains both scalar allocations");
#else
 // The second replacement may reuse the freed first allocation.
 check(first_address != old_first && second_address != old_second,
       "Bitmap fixed-width growth replaces each row's scalar allocation");
#endif

 auto row = matrix.begin();
 const auto row_end = matrix.end();
 bool traversal_ok = true;
 size_t r = 0;
 for (; r < 3 && row != row_end; ++r, ++row) {
  traversal_ok &= &*row == &matrix[r];
  auto cell = row->begin();
  const auto cell_end = row->end();
  size_t c = 0;
  for (; c < 3 && cell != cell_end; ++c, ++cell) {
   traversal_ok &= &*cell == &matrix[r][c] && *cell == T(golden[r][c]);
  }
  traversal_ok &= c == 3 && cell == cell_end;
 }
 traversal_ok &= r == 3 && row == row_end;
 check(traversal_ok, "fresh iterators traverse retained values and zero-initialized new row");
 if (!traversal_ok) return;

 auto last_row = row_end;
 --last_row;
 auto last_cell = last_row->end();
 --last_cell;
 check(&*last_row == &matrix[2] && &*last_cell == &matrix[2][2],
       "fresh endpoints decrement to the new row and its final cell");
 *last_cell = T(71);
 check(matrix[2][2] == T(71) && matrix[2][1] == T(0) && matrix[1][2] == T(4),
       "fresh iterator mutation reaches only the selected new cell");
}
int main(){fixed_width_growth<int>();fixed_width_growth<float>();fixed_width_growth<double>();self_move_aliases<int>();self_move_aliases<float>();self_move_aliases<double>();self_copy_aliases<int>();self_copy_aliases<float>();self_copy_aliases<double>();copy_aliases<int>();copy_aliases<float>();copy_aliases<double>();move_aliases<int>();move_aliases<float>();move_aliases<double>();row_sort<int>();row_sort<float>();row_sort<double>();row_algorithms<int>();row_algorithms<float>();row_algorithms<double>();cases<int>();cases<float>();cases<double>();narrow_shapes<int>();narrow_shapes<float>();narrow_shapes<double>();std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;}
