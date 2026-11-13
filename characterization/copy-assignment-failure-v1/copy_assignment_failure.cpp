// Each executable includes exactly one immutable legacy implementation.
#include LEGACY_HEADER
#include <iostream>
#include <string>
struct Fault { int marker; };
struct Value {
    static int live, constructors, assignments, fail_constructor, fail_assignment;
    int value=0;
    Value() { if (++constructors==fail_constructor) throw Fault{101}; ++live; }
    Value(const Value& x):value(x.value) { ++live; }
    Value& operator=(const Value& x) {
        if (++assignments==fail_assignment) throw Fault{202};
        value=x.value; return *this;
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
void fill(M& m,int offset) {
    for(size_t r=0;r<m.rows();++r)for(size_t c=0;c<m.cols();++c)
        m[r][c].value=offset+sentinel(r,c);
}
bool values(const M& m,size_t rows,size_t cols,int offset) {
    if(m.rows()!=rows||m.cols()!=cols||m.size()!=rows*cols)return false;
    for(size_t r=0;r<rows;++r) {
        if(m[r].size()!=cols)return false;
        for(size_t c=0;c<cols;++c)if(m[r][c].value!=offset+sentinel(r,c))return false;
    }
    return true;
}
void one(bool construction,int slot) {
    Value::reset();
    const std::string label=(construction?"constructor ":"assignment ")+std::to_string(slot);
    {
        M source(2,3),dest(3,2);fill(source,0);fill(dest,-200);
        auto* source_address=&source[0][0];
#ifdef NEURO
        auto* dest_address=&dest[0][0];const auto dest_capacity=dest.capacity();
#endif
        const auto source_capacity=source.capacity();
        Value::reset();
        if(construction)Value::fail_constructor=slot;else Value::fail_assignment=slot;
        bool exact=false,wrong=false;
        try { dest=source; }catch(const Fault& e){exact=e.marker==(construction?101:202);}
        catch(...){wrong=true;}
        const int ctor_calls=Value::constructors,assign_calls=Value::assignments;
        Value::reset();
        const int row=(slot-1)/3;
        check(exact&&!wrong,label+" exact_exception");
        check(ctor_calls==(construction?slot:(row+1)*3),label+" constructor_count");
        check(assign_calls==(construction?row*3:slot),label+" assignment_count");
        check(values(source,2,3,0),label+" source_shape_values");
        check(source.capacity()==source_capacity,label+" source_capacity");
        check(&source[0][0]==source_address,label+" source_storage_identity");
#ifdef NEURO
        check(dest.rows()==3&&dest.cols()==2&&dest.size()==6,label+" destination_metadata");
        check(dest.capacity()==dest_capacity,label+" destination_capacity");
        check(values(dest,3,2,-200),label+" destination_values");
        check(&dest[0][0]==dest_address,label+" destination_storage_identity");
        check(Value::live==12,label+" temporary_cleanup_live_balance");
#else
        check(dest.rows()==2&&dest.cols()==3&&dest.size()==6,label+" destination_metadata");
        check(dest.capacity()==source.capacity(),label+" destination_capacity");
        bool sizes=true,observed=true;
        for(int r=0;r<2;++r) {
            sizes=sizes&&dest[r].size()==size_t(r<=row?3:0);
            // Constructor failure leaves the current row size=3 but storage null.
            // Never index it. Only completed rows, or allocated assignment rows, are safe.
            if(r<row||(!construction&&r==row))for(int c=0;c<3;++c) {
                const int expected=(r<row||c<(slot-1)%3)?sentinel(r,c):0;
                observed=observed&&dest[r][c].value==expected;
            }
        }
        check(sizes,label+" partial_destination_row_sizes");
        check(observed,label+" safely_accessible_partial_values");
        check(Value::live==6+3*(row+(construction?0:1)),label+" partial_destination_live_balance");
#endif
        M* returned=&(dest=source);
        check(returned==&dest,label+" retry_returns_destination");
        check(values(dest,2,3,0),label+" retry_shape_values");
        check(&dest[0][0]!=&source[0][0],label+" retry_independent_storage");
        check(Value::live==12,label+" retry_live_balance");
    }
    check(Value::live==0,label+" final_destructor_balance");
}
void self_assignment() {
    Value::reset();
    {
        M m(2,3);fill(m,0);auto* address=&m[0][0];const auto capacity=m.capacity();
        Value::reset();Value::fail_constructor=Value::fail_assignment=1;
        bool threw=false;M* result=nullptr;
        try{result=&(m=m);}catch(...){threw=true;}
        check(!threw&&result==&m,"self return_no_exception");
        check(Value::constructors==0&&Value::assignments==0,"self no_element_operations");
        Value::reset();
        check(values(m,2,3,0)&&m.capacity()==capacity&&&m[0][0]==address,"self state_preserved");
        check(Value::live==6,"self live_balance");
    }
    check(Value::live==0,"self final_destructor_balance");
}
int main() {
    for(int i=1;i<=6;++i)one(true,i);
    for(int i=1;i<=6;++i)one(false,i);
    self_assignment();
    std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;
}
