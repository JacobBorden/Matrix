// Exactly one immutable legacy definition per executable.
#include LEGACY_HEADER
#include <iostream>
#include <string>
#include <vector>
struct Fault { int marker; };
struct Value {
    static int live, constructors, assignments, fail_constructor, fail_assignment;
    int value=0;
    Value() { if (++constructors==fail_constructor) throw Fault{101}; ++live; }
    Value(const Value& other):value(other.value) { ++live; }
    Value& operator=(const Value& other) {
        if (++assignments==fail_assignment) throw Fault{202};
        value=other.value; return *this;
    }
    ~Value() { --live; }
    static void reset() { constructors=assignments=0; fail_constructor=fail_assignment=-1; }
};
int Value::live=0,Value::constructors=0,Value::assignments=0;
int Value::fail_constructor=-1,Value::fail_assignment=-1;
int passed=0,failed=0;
void check(bool ok,const std::string& s) {
    std::cout<<(ok?"PASS ":"FAIL ")<<s<<'\n';ok?++passed:++failed;
}
using M=Matrix::Matrix<Value>;
int sentinel(size_t r,size_t c) { return int(17+31*r+7*c); }
// Independent operation schedule for this specific 2x3 -> 3x5 resize.
// Each entry is one default construction or copy assignment, with its row stage.
struct Event { bool construction; int row; };
std::vector<Event> schedule() {
    std::vector<Event> events;
    auto add=[&](bool construction,int row,int n) {
        for(int i=0;i<n;++i)events.push_back({construction,row});
    };
#ifdef NEURO
    add(true,0,5);add(false,0,3);
    add(true,1,5);add(false,1,3);
    add(true,2,5);
#else
    // Bitmap moves both original rows first, then constructs/fills the new row.
    add(true,2,5);add(true,2,1);add(false,2,5);
    add(true,0,5);add(false,0,3);
    add(true,1,5);add(false,1,3);
    add(true,2,5);add(false,2,5); // Even equal-size rows are reallocated.
#endif
    return events;
}
void one(const Event& event,int ctor_count,int assign_count) {
    Value::reset();
    const std::string label=(event.construction?"constructor ":"assignment ")+
        std::to_string(event.construction?ctor_count:assign_count);
    {
        M matrix(2,3);
        for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)matrix[r][c].value=sentinel(r,c);
        auto* second_row_address=&matrix[1][0];const auto capacity=matrix.capacity();
        Value::reset();
        if(event.construction)Value::fail_constructor=ctor_count;
        else Value::fail_assignment=assign_count;
        bool exact=false,wrong=false;
        try { matrix.resize(3,5); }catch(const Fault& e){exact=e.marker==(event.construction?101:202);}
        catch(...){wrong=true;}
        const int actual_ctors=Value::constructors,actual_assigns=Value::assignments;
        Value::reset();
        check(exact&&!wrong,label+" exact_exception");
        check(actual_ctors==ctor_count,label+" exact_constructor_count");
        check(actual_assigns==assign_count,label+" exact_assignment_count");
        check(matrix.rows()==2&&matrix.cols()==3&&matrix.size()==6,label+" unchanged_matrix_metadata");
        check(matrix.capacity()==capacity,label+" unchanged_matrix_capacity");
#ifdef NEURO
        const bool survives=event.row==0;
#else
        const bool survives=false;
#endif
        check(matrix[0].size()==0&&matrix[0].capacity()==0,label+" first_row_moved_out");
        check(matrix[1].size()==size_t(survives?3:0),label+" second_row_size");
        check(matrix[1].capacity()==size_t(survives?
#ifdef NEURO
            3*sizeof(Value)
#else
            3
#endif
            :0),label+" second_row_capacity");
        bool values=true,identity=true;
        // Never index moved-out rows or rely on the stale matrix column count.
        if(survives) {
            if(matrix[1].size()!=3) { values=false;identity=false; }
            else {
                for(size_t c=0;c<3;++c)values=values&&matrix[1][c].value==sentinel(1,c);
                identity=&matrix[1][0]==second_row_address;
            }
        }
        check(values,label+" safely_observable_surviving_values");
        check(identity,label+" surviving_storage_identity");
        check(Value::live==(survives?3:0),label+" temporary_and_moved_storage_destroyed");
        matrix.resize(3,5);
        bool shape=matrix.rows()==3&&matrix.cols()==5&&matrix.size()==15;
        bool retry_values=true;
        for(size_t r=0;r<3;++r) {
            shape=shape&&matrix[r].size()==5;
            if(matrix[r].size()!=5) { retry_values=false;continue; }
            for(size_t c=0;c<5;++c) {
                const int expected=survives&&r==1&&c<3?sentinel(r,c):0;
                retry_values=retry_values&&matrix[r][c].value==expected;
            }
        }
        check(shape,label+" retry_restores_consistent_shape");
        check(retry_values,label+" retry_preserves_only_surviving_values");
        check(Value::live==15,label+" retry_live_balance");
    }
    check(Value::live==0,label+" final_destructor_balance");
}
int main() {
    int ctors=0,assigns=0;
    for(const auto& event:schedule()) {
        if(event.construction)++ctors;else ++assigns;
        one(event,ctors,assigns);
    }
    check(Value::live==0,"all_cases_destroyed");
    std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;
}
