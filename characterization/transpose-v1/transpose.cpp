// Compile each immutable legacy definition in its own executable.
#include LEGACY_HEADER
#include <iostream>
#include <string>
#include <type_traits>
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
    static_assert(std::is_same<decltype(std::declval<const M&>().Transpose()),M>::value,"const transpose returns value");
    for(auto dims : {std::pair<size_t,size_t>{2,3},{3,2},{1,4},{4,1},{1,1},{3,3}}) {
        auto r=dims.first,c=dims.second;
        M source(r,c);
        for(size_t i=0;i<r;++i)for(size_t j=0;j<c;++j)source[i][j]=T(int(i*17+j*3)-11);
        const M& input=source;
        auto result=input.Transpose();
        auto label=type+" "+std::to_string(r)+"x"+std::to_string(c);
        check(shape(result,c,r),label+" dimensions");
        bool mapped=true,unchanged=true;
        for(size_t i=0;i<r;++i)for(size_t j=0;j<c;++j) {
            mapped=mapped && result[j][i]==T(int(i*17+j*3)-11);
            unchanged=unchanged && source[i][j]==T(int(i*17+j*3)-11);
        }
        check(mapped,label+" independent_golden_mapping");
        check(unchanged,label+" const_source_unchanged");
        auto twice=result.Transpose(); bool restored=shape(twice,r,c);
        for(size_t i=0;i<r;++i)for(size_t j=0;j<c;++j)restored=restored && twice[i][j]==T(int(i*17+j*3)-11);
        check(restored,label+" double_transpose");
        result[0][0]=T(99); source[r-1][c-1]=T(88);
        check(twice[0][0]==T(-11) && result[c-1][r-1]==(r*c==1?T(99):T(int((r-1)*17+(c-1)*3)-11)),label+" independent_storage");
    }
    for(auto dims : {std::pair<size_t,size_t>{0,0},{0,3},{2,0}}) {
        const M source(dims.first,dims.second);auto result=source.Transpose();
#ifdef NEURO
        check(shape(result,dims.second,dims.first),type+" empty_swaps_dimensions");
#else
        check(shape(result,0,0),type+" empty_normalizes_dimensions");
#endif
    }
}
int main() {
    cases<int>("int");cases<float>("float");cases<double>("double");
    std::cout << "RESULT passed=" << passed << " failed=" << failed << '\n';
    return failed ? 1 : 0;
}
