// Exactly one immutable legacy definition per executable.
#include LEGACY_HEADER
#include <iostream>
#include <string>
struct Fault { int marker; };
struct Value {
    static int live, constructors, assignments, fail_constructor, fail_assignment;
    int value=0;
    Value() {
        if (++constructors==fail_constructor) throw Fault{101};
        ++live;
    }
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
void check(bool ok,const std::string& label) {
    std::cout<<(ok?"PASS ":"FAIL ")<<label<<'\n'; ok?++passed:++failed;
}
using Row=Matrix::MatrixRow<Value>;
bool original(const Row& row) {
    return row.size()==3 && row[0].value==17 && row[1].value==-31 && row[2].value==59;
}
void one(size_t target,bool construction,int slot) {
    Value::reset();
    const std::string label=std::to_string(target)+(construction?" construction ":" assignment ")+std::to_string(slot);
    {
        Row row(3); row[0].value=17;row[1].value=-31;row[2].value=59;
        auto* address=&row[0];const auto capacity=row.capacity();const int live=Value::live;
        Value::reset();
        if(construction)Value::fail_constructor=slot;else Value::fail_assignment=slot;
        bool exact=false,wrong=false;
        try { row.resize(target); } catch(const Fault& e) { exact=e.marker==(construction?101:202); }
        catch(...) { wrong=true; }
        const int calls=construction?Value::constructors:Value::assignments;
        Value::reset();
        check(exact&&!wrong,label+" exact_exception");
        check(calls==slot,label+" injection_reached");
        check(original(row),label+" original_shape_values");
        check(row.capacity()==capacity,label+" original_capacity");
        check(&row[0]==address,label+" original_storage_identity");
        check(Value::live==live,label+" temporary_elements_destroyed");
        row.resize(target);
        bool values=row.size()==target;
        const int expected[]={17,-31,59};
        for(size_t i=0;i<target;++i)values=values&&row[i].value==(i<3?expected[i]:0);
        check(values,label+" retry_succeeds_values");
        check(Value::live==int(target),label+" retry_live_balance");
    }
    check(Value::live==0,label+" final_destructor_balance");
}
int main() {
    // Grow: five constructor sites and three assignment sites.
    // Shrink: two constructor sites and two assignment sites.
    for(int i=1;i<=5;++i)one(5,true,i);
    for(int i=1;i<=3;++i)one(5,false,i);
    for(int i=1;i<=2;++i)one(2,true,i);
    for(int i=1;i<=2;++i)one(2,false,i);
    check(Value::live==0,"all_cases_destroyed");
    std::cout<<"RESULT passed="<<passed<<" failed="<<failed<<'\n';return failed?1:0;
}
