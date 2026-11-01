// One immutable legacy header per executable; no arithmetic APIs instantiated for Value.
#include LEGACY_HEADER
#include <iostream>
#include <string>
#include <utility>
struct Value {
    static int live;
    int value=0;
    Value() { ++live; }
    Value(int n):value(n) { ++live; }
    Value(const Value& v):value(v.value) { ++live; }
    Value& operator=(const Value&)=default;
    ~Value() { --live; }
    bool operator==(const Value& v) const { return value==v.value; }
};
int Value::live=0;
int passed=0,failed=0;
void check(bool ok,const std::string& name) {
    std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n'; ok?++passed:++failed;
}
template<class T> bool original(const Matrix::Matrix<T>& m) {
    if(m.rows()!=2 || m.cols()!=3 || m.size()!=6) return false;
    for(size_t r=0;r<2;++r) for(size_t c=0;c<3;++c)
        if(!(m[r][c]==T(int(r*3+c+1)))) return false;
    return true;
}
template<class T> bool cleared(const Matrix::Matrix<T>& m) {
    return m.rows()==0 && m.cols()==0 && m.size()==0 && m.capacity()==0;
}
template<class T> void cases(const std::string& type) {
    using M=Matrix::Matrix<T>;
    M empty;
    check(cleared(empty),type+" default_empty");
    M a(2,3);
    for(size_t r=0;r<2;++r) for(size_t c=0;c<3;++c) a[r][c]=T(int(r*3+c+1));
    check(original(a),type+" shape_values");
#ifdef NEURO
    check(a.capacity()==6*sizeof(T) && a[0].capacity()==3*sizeof(T),type+" capacity_bytes");
#else
    check(a.capacity()==2 && a[0].capacity()==3,type+" capacity_rows_elements");
#endif
    M copy(a); copy[0][0]=T(99);
    check(original(a) && copy[0][0]==T(99) && &a[0][0]!=&copy[0][0],type+" copy_constructor_independent");
    M assigned(1,1); assigned=a; assigned[1][2]=T(88);
    check(original(a) && assigned[1][2]==T(88),type+" copy_assignment_independent");
    auto* address=&a[0][0]; a=a;
    check(original(a) && address==&a[0][0],type+" self_copy_preserves");
    M moving(a); auto* moved_address=&moving[0][0]; M moved(std::move(moving));
    bool move_ok=original(moved) && cleared(moving) && &moved[0][0]==moved_address;
    moving.resize(1,1); moving[0][0]=T(7);
    check(move_ok && moving[0][0]==T(7),type+" move_constructor_source_reusable");
    M target(1,1); target=std::move(moved);
    check(original(target) && cleared(moved) && &target[0][0]==moved_address,type+" move_assignment_transfers");
    target=std::move(target);
    check(original(target) && &target[0][0]==moved_address,type+" self_move_preserves");
    target.resize(3,4); bool grew=target.rows()==3 && target.cols()==4 && target.size()==12;
    for(size_t r=0;r<3;++r) for(size_t c=0;c<4;++c)
        grew=grew && target[r][c]==T(r<2 && c<3?int(r*3+c+1):0);
    check(grew,type+" grow_preserves_rectangle_defaults");
    target.resize(1,2);
    check(target.rows()==1 && target.cols()==2 && target[0][0]==T(1) && target[0][1]==T(2),type+" shrink_preserves_rectangle");
    target.resize(2,0);
#ifdef NEURO
    check(cleared(target),type+" zero_column_resize_normalizes");
#else
    check(target.rows()==2 && target.cols()==0 && target.size()==0 && target.capacity()==2 && target[0].size()==0,type+" zero_column_resize_retains_rows");
#endif
    bool bounds=false; try { (void)a[0].at(3); } catch(const std::out_of_range&) { bounds=true; }
    check(bounds,type+" row_at_throws");
}
int main() {
    cases<float>("float"); cases<double>("double"); cases<int>("int");
    cases<Value>("nonnumeric_Value");
    check(Value::live==0,"nonnumeric_all_instances_destroyed");
    std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';
    return failed?1:0;
}
