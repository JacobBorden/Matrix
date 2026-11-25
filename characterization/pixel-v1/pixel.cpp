#include "bitmap/bitmap.h"

#include <array>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
using Image = Matrix::Matrix<Pixel>;

static_assert(std::is_same_v<decltype(std::declval<Image&>().at(0, 0)), Pixel&>);
static_assert(std::is_same_v<decltype(std::declval<const Image&>().at(0, 0)), Pixel>);

constexpr Pixel clear{};
constexpr Pixel a{1, 23, 45, 67};
constexpr Pixel b{89, 101, 123, 145};
constexpr Pixel c{167, 189, 201, 223};
constexpr Pixel d{245, 254, 0, 255};
constexpr Pixel fill{31, 62, 93, 124};
constexpr Pixel changed{17, 34, 51, 68};

int passed = 0;
int failed = 0;

void check(bool condition, std::string_view label)
{
    if (condition) {
        ++passed;
    } else {
        ++failed;
        std::cerr << "FAIL " << label << '\n';
    }
}

bool equal(Pixel lhs, Pixel rhs)
{
    return lhs.blue == rhs.blue && lhs.green == rhs.green
        && lhs.red == rhs.red && lhs.alpha == rhs.alpha;
}

template <std::size_t Rows, std::size_t Cols>
bool matches(const Image& image, const Pixel (&expected)[Rows][Cols])
{
    if (image.rows() != Rows || image.cols() != Cols
        || image.size() != Rows * Cols || image.capacity() != Rows || image.empty()) {
        return false;
    }
    for (std::size_t row = 0; row < Rows; ++row) {
        for (std::size_t col = 0; col < Cols; ++col) {
            if (!equal(image.at(row, col), expected[row][col])) {
                return false;
            }
        }
    }
    return true;
}

bool is_empty(const Image& image)
{
    return image.empty() && image.rows() == 0 && image.cols() == 0
        && image.size() == 0 && image.capacity() == 0;
}

auto row_storage(Image& image)
{
    return std::array{&image.at(0, 0), &image.at(1, 0)};
}

void check_bitmap_conversion()
{
    constexpr Pixel e{13, 29, 53, 255};
    constexpr Pixel f{71, 97, 113, 255};
    const Pixel first_row[][3]{{{1, 23, 45, 255}, {89, 101, 123, 255}, e},
                               {{167, 189, 201, 255}, {245, 254, 0, 255}, f}};
    const Pixel last_row[][3]{{{167, 189, 201, 255}, {245, 254, 0, 255}, f},
                              {{1, 23, 45, 255}, {89, 101, 123, 255}, e}};

    Bitmap::File input;
    input.bitmapInfoHeader.biWidth = 3;
    input.bitmapInfoHeader.biHeight = -2;
    input.bitmapInfoHeader.biBitCount = 24;
    input.bitmapData = {1, 23, 45, 89, 101, 123, 13, 29, 53, 240, 241, 242,
                        167, 189, 201, 245, 254, 0, 71, 97, 113, 243, 244, 245};
    check(matches(CreateMatrixFromBitmap(input), first_row),
          "24-bit top-down input retains stored rows and skips padding");

    input.bitmapInfoHeader.biHeight = 2;
    input.bitmapData = {167, 189, 201, 245, 254, 0, 71, 97, 113, 243, 244, 245,
                        1, 23, 45, 89, 101, 123, 13, 29, 53, 240, 241, 242};
    check(matches(CreateMatrixFromBitmap(input), last_row),
          "24-bit bottom-up input retains stored rows and skips padding");
    input.bitmapData.pop_back();
    check(is_empty(CreateMatrixFromBitmap(input)),
          "24-bit input missing a byte of final row padding yields an empty matrix");

    Image image(2, 3);
    const Pixel expected[][3]{{a, b, e}, {c, d, f}};
    for (std::size_t row = 0; row < 2; ++row) {
        for (std::size_t col = 0; col < 3; ++col) {
            image.at(row, col) = expected[row][col];
        }
    }
    Bitmap::File output = CreateBitmapFromMatrix(image);
    const auto& info = output.bitmapInfoHeader;
    check(output.IsValid() && info.biSize == 40 && info.biWidth == 3 && info.biHeight == 2
              && info.biPlanes == 1 && info.biBitCount == 32 && info.biCompression == 0
              && info.biSizeImage == 24,
          "output has a positive-height uncompressed 32-bit info header");
    const auto& file = output.bitmapFileHeader;
    check(file.bfType == 0x4D42 && file.bfOffBits == 54 && file.bfSize == 78
              && file.bfReserved1 == 0 && file.bfReserved2 == 0,
          "output file header has the expected signature, size and pixel offset");
    const std::vector<std::uint8_t> bytes{1, 23, 45, 67, 89, 101, 123, 145, 13, 29, 53, 255,
                                          167, 189, 201, 223, 245, 254, 0, 255, 71, 97, 113, 255};
    check(output.bitmapData == bytes, "output stores exact BGRA bytes in matrix row order");
    check(matches(CreateMatrixFromBitmap(output), expected) && matches(image, expected),
          "32-bit round trip preserves every channel and leaves source pixels intact");
}

void check_copy(Image& copy, Image& source)
{
    const Pixel original[][2]{{a, b}, {c, d}};
    check(matches(copy, original), "copy retains shape and channels");
    const auto copy_rows = row_storage(copy);
    const auto source_rows = row_storage(source);
    check(copy_rows[0] != source_rows[0] && copy_rows[1] != source_rows[1],
          "copy owns both rows independently");
    copy.at(0, 1) = changed;
    const Pixel updated[][2]{{a, changed}, {c, d}};
    check(matches(copy, updated) && matches(source, original), "copy mutation stays local");
}
}

int main()
{
    Image source(2, 2);
    const Pixel blank[][2]{{clear, clear}, {clear, clear}};
    check(matches(source, blank), "default pixels have zeroed BGRA channels");

    source.at(0, 0) = a;
    source[0][1] = b;
    source.at(1, 0) = c;
    source[1][1] = d;
    const Pixel original[][2]{{a, b}, {c, d}};
    check(matches(source, original), "mutable indexing and const reads retain all channels");
    check(&source.at(1, 1) == &source[1][1], "checked and unchecked access alias");

    const Image& const_source = source;
    Pixel value = const_source.at(0, 0);
    value.alpha = changed.alpha;
    check(value.alpha != source.at(0, 0).alpha && matches(source, original),
          "const access returns an independent value");

    bool rejected = false;
    try {
        (void)source.at(0, 2);
    } catch (const std::out_of_range&) {
        rejected = true;
    }
    check(rejected && matches(source, original), "checked access rejects a column past the end");

    Image copied(source);
    check_copy(copied, source);
    Image assigned(1, 3, fill);
    check(&(assigned = source) == &assigned, "copy assignment returns the destination");
    check_copy(assigned, source);

    const auto storage = row_storage(source);
    Image moved(std::move(source));
    check(is_empty(source), "move construction empties the source");
    check(matches(moved, original) && row_storage(moved) == storage,
          "move construction transfers both rows");

    Image destination(1, 1, fill);
    check(&(destination = std::move(moved)) == &destination,
          "move assignment returns the destination");
    check(is_empty(moved), "move assignment empties the source");
    check(matches(destination, original) && row_storage(destination) == storage,
          "move assignment transfers both rows");
    source.resize(1, 1, fill);
    const Pixel reused[][1]{{fill}};
    check(matches(source, reused) && matches(destination, original),
          "moved-from matrix can be reused independently");

    destination.resize(3, 3);
    const Pixel grown[][3]{{a, b, clear}, {c, d, clear}, {clear, clear, clear}};
    check(matches(destination, grown), "growth retains overlap and zeroes new pixels");
    destination.resize(2, 1);
    const Pixel shrunk[][1]{{a}, {c}};
    check(matches(destination, shrunk), "shrink retains the exact upper-left region");
    destination.resize(3, 2, fill);
    const Pixel regrown[][2]{{a, fill}, {c, fill}, {fill, fill}};
    check(matches(destination, regrown), "fill resize initializes only newly exposed pixels");

    check_bitmap_conversion();

    std::cout << "RESULT passed=" << passed << " failed=" << failed << '\n';
    return failed != 0;
}
