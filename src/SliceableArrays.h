/*
 * MIT License
 *
 * Copyright (c) 2017 Daniel Politte
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files(the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions :
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <iostream>
#include <cassert>
#include <numeric>
#include <vector>
#include <array>

#ifndef __H_SLICEABLEARRAYS
#define __H_SLICEABLEARRAYS

namespace SliceableArrays {

   /*
    * Sliceable Arrays
    *
    * Features:
    * - 0-indexed
    * - Column-major
    * - Arbitrarily high dimensionality, known at compile time
    * - Fixed-size
    * - Dimensions need not be known until runtime
    *
    * descended from class described at:
    * http://www.cplusplus.com/forum/articles/17108/
    *
    * AN IMPORTANT DETAIL different from usual C arrays:
    * The elements are stored internally in 2nd dimension-major order. That is, a
    * 2D array with dimensions n by m has a memory layout of m series n elements,
    * in which each series shares a common 2nd index. The first index changes the
    * fastest as memory is traversed. This also applies to the higher-dimensioned
    * arrays; the first dimension changes the fastest, and the last dimension
    * changes most slowly.
    *
    */

    template <typename T, size_t NDIMS>
    class ArrayND {
    private:
        const std::array<size_t, NDIMS> dims_;
        const size_t numEls_;
        bool dataAutoDestroy_;
        T *arr_;

        size_t getFlattenedIndex(std::vector<size_t>& position) const {
            // make sure we have the proper number of indices
            if (position.size() != NDIMS) {
                std::cerr << "Illegally-sized index into ArrayND!" << std::endl;
            }
            position.resize(NDIMS, 0); // tack on additional zeros as necc.

            // If we're debugging, check out validity of each position vector item (i.e., less than its dim's size)
            assert(indexValid(position, dims_));

            size_t index = 0;
            for (size_t i = dims_.size() - 1; i > 0; --i) {
                index += position[i];
                index *= dims_[i - 1];
            }
            index += position[0];
            return index;
        }

        size_t buildEndOfIndex(size_t posThis) {
            // The tail
            return posThis;
        }

        template<typename... Args>
        size_t buildEndOfIndex(size_t posThis, Args... posRest) {
            const size_t correspondingDimensionIndex = NDIMS - (1 + sizeof...(posRest));
            // If we could have the full position list as a vector, posThis would be pos[correspondingDimensionIndex]
            size_t tailResult = buildEndOfIndex(posRest...);
            return posThis + (dims_[correspondingDimensionIndex] * tailResult);
        }

        static bool indexValid(const std::vector<size_t>& pos, const std::array<size_t, NDIMS>& dims) {
            for (size_t i = 0; i < pos.size(); ++i) {
                if (pos[i] >= dims[i]) {
                    return false;
                }
            }
            return true;
        }

        static void compactArgsToArrayInner(std::array<size_t, NDIMS>& resultArray, size_t only) {
            resultArray[NDIMS - 1] = only;
        }

        template<typename... Args>
        static void compactArgsToArrayInner(std::array<size_t, NDIMS>& resultArray, size_t first, Args... vals) {
            resultArray[NDIMS - (sizeof...(vals)+1)] = first; // the +1 is to account for first not being counted among the args
            compactArgsToArrayInner(resultArray, vals...);
        }

        template<typename... Args>
        static std::array<size_t, NDIMS> compactArgsToArray(Args... vals) {
            static_assert(sizeof...(vals) == NDIMS, "Requires NDIMS arguments exactly");
            std::array<size_t, NDIMS> resultArray; // an array we'll gradually fill with copies of the relavent values
            compactArgsToArrayInner(resultArray, vals...);
            return resultArray;
        }

        // allocates space on heap for this matrix
        void initializeData() {
            if (numEls_ > 0) {
                arr_ = new T[numEls_];
            }
        }

        // Deallocate space we used, if applicable
        void destroyDataIfNecessary() {
            if (dataAutoDestroy_) {
                delete[] arr_;
            }
        }

    public:
        /*
        * Constructor & Deconstructor: The ArrayND class is responsible for the
        * management of the space allocated/adopted for the matrix's data
        */

        // Constructor in which user doesn't provide a data pointer, and space is automatically alloc'd and marked for destruction
        template<typename... Args>
        ArrayND(size_t dim1, Args... remainingDims)
            : dims_(compactArgsToArray(dim1, remainingDims...)),
            numEls_(accumulate(dims_.begin(), dims_.end(), 1, [](size_t a, size_t b) -> size_t {return a*b; })),
            dataAutoDestroy_(true),
            arr_(0)
        {
            static_assert(NDIMS == sizeof...(remainingDims)+1, "ArrayND constructor requires exactly as many size arguments as the array has dimensions");
            initializeData();
        }

        // Constructor in which user provides a data pointer, and that data is not marked for destruction
        template<typename... Args>
        ArrayND(T* dataPtr, size_t dim1, Args... remainingDims)
            : dims_(compactArgsToArray(dim1, remainingDims...)),
            numEls_(accumulate(dims_.begin(), dims_.end(), 1, [](size_t a, size_t b) -> size_t {return a*b; })),
            dataAutoDestroy_(false),
            arr_(dataPtr)
        {
            static_assert(NDIMS == sizeof...(remainingDims)+1, "ArrayND constructor requires exactly as many size arguments as the array has dimensions");
        }

        // copy constructor: deep copy of data
        ArrayND(ArrayND& other)
            : arr_(0),
            dataAutoDestroy_(true), // since we always allocate space here
            numEls_(other.numEls_), dims_(other.dims_)
        {
            initializeData(); // allocating space automatically
            for (size_t i = 0; i < numEls_; ++i) {
                arr_[i] = other.arr_[i];
            }
        }

        // copy assignment operator: deep copy of data
        ArrayND& operator=(ArrayND& other) {
            dataAutoDestroy_ = true; // since we always allocate space here
            numEls_ = other.numEls_;
            dims_ = other.dims_;
            initializeData(); // allocating space automatically
            for (size_t i = 0; i < numEls_; ++i) {
                arr_[i] = other.arr_[i];
            }
            return *this;
        }

        // move constructor
        ArrayND(ArrayND&& other) noexcept
            : arr_(other.arr_), dataAutoDestroy_(other.dataAutoDestroy_),
            numEls_(other.numEls_), dims_(other.dims_)
        {
            // get old one to state where running destructor is safe
            other.arr_ = nullptr;
        }

        // move assignment operator
        ArrayND& operator=(ArrayND&& other) {
            if (this != &rhs) {
                destroyDataIfNecessary(); // flush data that was at new location, if any

                arr_ = other.arr_;
                dataAutoDestroy_ = other.dataAutoDestroy_;
                numEls_ = other.numEls_;
                dims_ = other.dims_;

                // get old one to state where running destructor is safe
                other.arr_ = nullptr;
            }
            return *this;
        }

        // destructor
        ~ArrayND()
        {
            destroyDataIfNecessary();
        }

        // get dims
        size_t numEls() const { return numEls_; }

        // Get the number of elements along the nth dimension (counting from 0, of course)
        size_t getDim(size_t n) const {
            assert(n < dims_.size());
            return dims_[n];
        }

        // 1-D indexing, which is common to all dimensionalities
        T ind(size_t i) const {
            assert(i >= 0 && i < numEls_);
            return arr_[i];
        }
        T& ind(size_t i) {
            assert(i >= 0 && i < numEls_);
            return arr_[i];
        }
        T operator[](size_t i) const {
            return ind(i);
        }
        T& operator[](size_t i) {
            return ind(i);
        }

        // non-general multi-D indexing, only available to certain dimensionalities
        T ind(size_t p0, size_t p1) const {
            static_assert(NDIMS == 2, "Only 2D arrays can be indexed with 2 dimensions");
            return ind(getOffsetAtIndex(p0, p1));
        }
        T& ind(size_t p0, size_t p1) {
            static_assert(NDIMS == 2, "Only 2D arrays can be indexed with 2 dimensions");
            return ind(getOffsetAtIndex(p0, p1));
        }
        size_t getOffsetAtIndex(size_t p0, size_t p1) const {
            static_assert(NDIMS == 2, "Only 2D arrays can be indexed with 2 dimensions");
            size_t m = getDim(0);
            size_t n = getDim(1);
            assert(p0 < m && p1 < n);
            return p1*m + p0;
        }

        T ind(size_t p0, size_t p1, size_t p2) const {
            static_assert(NDIMS == 3, "Only 3D arrays can be indexed with 3 dimensions");
            return ind(getOffsetAtIndex(p0, p1, p2));
        }
        T& ind(size_t p0, size_t p1, size_t p2) {
            static_assert(NDIMS == 3, "Only 3D arrays can be indexed with 3 dimensions");
            return ind(getOffsetAtIndex(p0, p1, p2));
        }
        size_t getOffsetAtIndex(size_t p0, size_t p1, size_t p2) const {
            static_assert(NDIMS == 3, "Only 3D arrays can be indexed with 3 dimensions");
            size_t m = getDim(0);
            size_t n = getDim(1);
            size_t p = getDim(2);
            assert(p0 < m && p1 < n && p2 < p);
            return (p2*n + p1)*m + p0;
        }

        T ind(size_t p0, size_t p1, size_t p2, size_t p3) const {
            static_assert(NDIMS == 4, "Only 4D arrays can be indexed with 4 dimensions");
            return ind(getOffsetAtIndex(p0, p1, p2, p3));
        }
        T& ind(size_t p0, size_t p1, size_t p2, size_t p3) {
            static_assert(NDIMS == 4, "Only 4D arrays can be indexed with 4 dimensions");
            return ind(getOffsetAtIndex(p0, p1, p2, p3));
        }
        size_t getOffsetAtIndex(size_t p0, size_t p1, size_t p2, size_t p3) const {
            static_assert(NDIMS == 4, "Only 4D arrays can be indexed with 4 dimensions");
            size_t m = getDim(0);
            size_t n = getDim(1);
            size_t p = getDim(2);
            size_t q = getDim(3);
            assert(p0 < m && p1 < n && p2 < p && p3 < q);
            return ((p3*p + p2)*n + p1)*m + p0;
        }

        T ind(size_t p0, size_t p1, size_t p2, size_t p3, size_t p4) const {
            static_assert(NDIMS == 5, "Only 5D arrays can be indexed with 5 dimensions");
            return ind(getOffsetAtIndex(p0, p1, p2, p3, p4));
        }
        T& ind(size_t p0, size_t p1, size_t p2, size_t p3, size_t p4) {
            static_assert(NDIMS == 5, "Only 5D arrays can be indexed with 5 dimensions");
            return ind(getOffsetAtIndex(p0, p1, p2, p3, p4));
        }
        size_t getOffsetAtIndex(size_t p0, size_t p1, size_t p2, size_t p3, size_t p4) const {
            static_assert(NDIMS == 5, "Only 5D arrays can be indexed with 5 dimensions");
            size_t m = getDim(0);
            size_t n = getDim(1);
            size_t p = getDim(2);
            size_t q = getDim(3);
            size_t r = getDim(4);
            assert(p0 < m && p1 < n && p2 < p && p3 < q && p4 < r);
            return (((p4*q + p3)*p + p2)*n + p1)*m + p0;
        }

        T ind(size_t p0, size_t p1, size_t p2, size_t p3, size_t p4, size_t p5) const {
            static_assert(NDIMS == 6, "Only 6D arrays can be indexed with 6 dimensions");
            return ind(getOffsetAtIndex(p0, p1, p2, p3, p4, p5));
        }
        T& ind(size_t p0, size_t p1, size_t p2, size_t p3, size_t p4, size_t p5) {
            static_assert(NDIMS == 6, "Only 6D arrays can be indexed with 6 dimensions");
            return ind(getOffsetAtIndex(p0, p1, p2, p3, p4, p5));
        }
        size_t getOffsetAtIndex(size_t p0, size_t p1, size_t p2, size_t p3, size_t p4, size_t p5) const {
            static_assert(NDIMS == 6, "Only 6D arrays can be indexed with 6 dimensions");
            size_t m = getDim(0);
            size_t n = getDim(1);
            size_t p = getDim(2);
            size_t q = getDim(3);
            size_t r = getDim(4);
            size_t s = getDim(5);
            assert(p0 < m && p1 < n && p2 < p && p3 < q && p4 < r && p5 < s);
            return ((((p5*r + p4)*q + p3)*p + p2)*n + p1)*m + p0;
        }

        T ind(size_t p0, size_t p1, size_t p2, size_t p3, size_t p4, size_t p5, size_t p6) const {
            static_assert(NDIMS == 7, "Only 7D arrays can be indexed with 7 dimensions");
            return ind(getOffsetAtIndex(p0, p1, p2, p3, p4, p5, p6));
        }
        T& ind(size_t p0, size_t p1, size_t p2, size_t p3, size_t p4, size_t p5, size_t p6) {
            static_assert(NDIMS == 7, "Only 7D arrays can be indexed with 7 dimensions");
            return ind(getOffsetAtIndex(p0, p1, p2, p3, p4, p5, p6));
        }
        size_t getOffsetAtIndex(size_t p0, size_t p1, size_t p2, size_t p3, size_t p4, size_t p5, size_t p6) const {
            static_assert(NDIMS == 7, "Only 7D arrays can be indexed with 7 dimensions");
            size_t m = getDim(0);
            size_t n = getDim(1);
            size_t p = getDim(2);
            size_t q = getDim(3);
            size_t r = getDim(4);
            size_t s = getDim(5);
            size_t t = getDim(6);
            assert(p0 < m && p1 < n && p2 < p && p3 < q && p4 < r && p5 < s && p6 < t);
            return (((((p6*s + p5)*r + p4)*q + p3)*p + p2)*n + p1)*m + p0;
        }


        // TODO: this should be destroyed if general version of 'ind' is good
        T indexGeneral(const std::vector<size_t>& position) const {
            // A slow but general indexing function, so larger dimensions can kinda stand on its own

            size_t index = getFlattenedIndex(position);
            return ind(index);
        }
        T& indexGeneral(const std::vector<size_t>& position) {
            // A slow but general indexing function, so larger dimensions can kinda stand on its own

            size_t index = getFlattenedIndex(position);
            return ind(index);
        }

        
        template<typename... Args>
        T ind(size_t posFirst, Args... posRest) const {
            // Generic offset computation for arbitrarily high dimensions.
            // This offset building is equiv to using getFlattenedIndex if we
            // had all the position elements in a container
            size_t offset = buildEndOfIndex(posFirst, posRest...);
            return ind(offset);
        }
        template<typename... Args>
        T& ind(size_t posFirst, Args... posRest) {
            // Generic offset computation for arbitrarily high dimensions.
            // This offset building is equiv to using getFlattenedIndex if we
            // had all the position elements in a container
            size_t offset = buildEndOfIndex(posFirst, posRest...);
            return ind(offset);
        }

        // reset entire matrix to a value
        void fill(T val) {
#pragma omp parallel for
            for (int i_tmp = 0; i_tmp < (int)numEls_; ++i_tmp) {
                size_t i = i_tmp; // To provide OpenMP 2.0 compatibility
                arr_[i] = val;
            }
        }

        // Retrieves the location of internal data
        T* getData() const {
            return arr_;
        }

        // Modifies the location of internal data. Will not autodestroy this new location by default
        void setData(T *data) {
            arr_ = data;
            dataAutoDestroy_ = true;
        }

        // choose whether the internal data of this array will be deleted when it is destructed.
        // The default value when this has not been called is true.
        void setDataAutoDestroy(bool isDataForfeit) {
            dataAutoDestroy_ = isDataForfeit;
        }

        bool getDataAutoDestroy() {
            return dataAutoDestroy_;
        }

    };

    template <typename T> using Array1D = ArrayND<T, 1>;
    template <typename T> using Array2D = ArrayND<T, 2>;
    template <typename T> using Array3D = ArrayND<T, 3>;
    template <typename T> using Array4D = ArrayND<T, 4>;
    template <typename T> using Array5D = ArrayND<T, 5>;
    template <typename T> using Array6D = ArrayND<T, 6>;
    template <typename T> using Array7D = ArrayND<T, 7>;

};

#endif /* __H_SLICEABLEARRAYS */
