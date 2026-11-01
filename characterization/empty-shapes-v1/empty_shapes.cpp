// Each pinned legacy implementation is compiled independently.
#include LEGACY_HEADER
#include <iostream>
#include <string>
#include <utility>
int passed=0, failed=0;
void check(bool ok, const std::string& label) {
    std::cout << (ok ? "PASS " : "FAIL ") << label << '\n';
    ok ? ++passed : ++failed;
}
template<class T> bool shape(const Matrix::Matrix<T>& m, size_t r, size_t c) {
    return m.rows()==r && m.cols()==c && m.size()==r*c;
}
template<class T> void cases(const std::string& type) {
    using M=Matrix::Matrix<T>;
    for (auto dims : {std::pair<size_t,size_t>{0,0}, {0,3}, {2,0}}) {
        auto r=dims.first, c=dims.second;
        const auto label=type+" "+std::to_string(r)+"x"+std::to_string(c);
        M source(r,c);
        check(shape(source,r,c),label+" constructor_preserves_dimensions");
        M copy(source);
        check(shape(copy,r,c),label+" copy_preserves_dimensions");
        M assigned(1,1); assigned[0][0]=T(97); assigned=source;
        check(shape(assigned,r,c),label+" assignment_preserves_dimensions");
        M moved(std::move(copy));
        check(shape(moved,r,c) && shape(copy,0,0),label+" move_transfers_dimensions");
        bool rows_empty=true;
        for(size_t i=0;i<r;++i) rows_empty=rows_empty && source[i].size()==0;
        check(rows_empty,label+" allocated_rows_are_empty");
        moved.resize(2,3);
        bool defaults=shape(moved,2,3);
        for(size_t i=0;i<2;++i)for(size_t j=0;j<3;++j)defaults=defaults && moved[i][j]==T(0);
        check(defaults,label+" growth_default_initializes");
        moved[1][2]=T(71);
        check(shape(source,r,c),label+" growth_independent_of_source");
        moved.resize(r,c);
#ifdef NEURO
        check(shape(moved,0,0),label+" resize_normalizes_empty");
#else
        check(shape(moved,r,c),label+" resize_preserves_dimensions");
#endif
        moved.resize(1,2);
        check(shape(moved,1,2) && moved[0][0]==T(0) && moved[0][1]==T(0),label+" regrowth_discards_old_values");
    }
}
int main() {
    cases<int>("int"); cases<float>("float"); cases<double>("double");
    std::cout << "RESULT passed=" << passed << " failed=" << failed << '\n';
    return failed ? 1 : 0;
}
