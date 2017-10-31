#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"
#include "../src/SliceableArrays.h"

using namespace SliceableArrays;

TEST_CASE("Specific instantiations of arrays can be declared directly") {
    Array2D<int> arr(3,4);
    Array3D<double> arr2(5, 7, 100);
}

TEST_CASE("Specific instantiations of arrays can be declared with external memory specified") {
    int mem[] = { 5, 9, 1, 6, 10, -5 };
    Array2D<int> arr(mem, 3, 2);
    arr.ind(2,1);
}

TEST_CASE("Super-dimensional arrays can be created and indexed") {
    ArrayND<int,10> arr(5,1,5,1,1,1,1,1,2,4);
    arr.ind(3,0,2,0,0,0,0,0,1,2);
}

TEST_CASE("Arrays fill function changes all elements") {
    constexpr size_t sz = 10;
    int fillVal = 4;
    Array1D<int> m(sz);
    m.fill(fillVal);

    for (size_t i = 0; i < sz; ++i) {
        REQUIRE(m[i] == fillVal);
    }
}

TEST_CASE("Arrays are indexed with the first index changing most quickly") {
    Array2D<int> mat(4, 4);
    for (size_t i = 0; i < 16; ++i) {
        mat[i] = i;
    }

    REQUIRE(mat.ind(2, 0) == 2);
    REQUIRE(mat.ind(0, 2) == 8);
}

TEST_CASE("Array copy constructor and copy assignment op create a deep copy") {
    Array1D<int> mat(4);
    mat.fill(0);
    Array1D<int> mat_copied(mat);
    mat_copied.ind(1) = 42;

    REQUIRE(mat.ind(1) == 0);

    Array1D<int> mat_copyAssigned(mat);
    mat_copyAssigned.ind(1) = 42;

    REQUIRE(mat.ind(1) == 0);
}
