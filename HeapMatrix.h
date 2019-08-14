#ifndef HEAP_MATRIX_H
#define HEAP_MATRIX_H
#include "MatrixInterface.h"

template<typename T> class HeapMatrix : public MatrixInterface<HeapMatrix<T>> {
    using Interface = MatrixInterface<HeapMatrix<T>>;
    friend Interface;

public:
    using ValueType = T;

    static const mtl::size_t Rows = Interface::DynamicValue;
    static const mtl::size_t Cols = Interface::DynamicValue;
    static const mtl::size_t EffCols = Interface::DynamicValue;

private:
    T* data_;
    mtl::size_t rowCount_;
    mtl::size_t columnCount_;

public:
    HeapMatrix() :
        data_(nullptr),
        rowCount_(0),
        columnCount_(0)
    {}
    HeapMatrix(mtl::size_t rowCount, mtl::size_t columnCount) :
        data_(new T[rowCount * columnCount]),
        rowCount_(rowCount),
        columnCount_(columnCount)
    {}
    HeapMatrix(T* newedPtr, mtl::size_t rowCount, mtl::size_t columnCount) :
        data_(newedPtr),
        rowCount_(rowCount),
        columnCount_(columnCount)
    {}
    template<typename MT> HeapMatrix(const MatrixInterface<MT>& m) :
        data_(new T[m.rowCount() * m.columnCount()]),
        columnCount_(m.columnCount()),
        rowCount_(m.rowCount())
    {
        this->copy(m);
    }
    template<typename MT> HeapMatrix(const MatrixInterface<MT>& m, int rCount, int cCount) :
        data_(new T[rCount * cCount]),
        columnCount_(cCount),
        rowCount_(rCount)
    {
        this->copy(m);
    }
    HeapMatrix(const HeapMatrix<T>& m) :
        data_(new T[m.rowCount() * m.columnCount()]),
        columnCount_(m.columnCount()),
        rowCount_(m.rowCount())
    {
        this->copy(m);
    }
    HeapMatrix(HeapMatrix<T>&& m) : 
        data_(m.data_),
        columnCount_(m.columnCount()),
        rowCount_(m.rowCount())
    {
        m.data_ = nullptr;
    }
    ~HeapMatrix() {
        delete[] data_;
    }

    template<typename MT> HeapMatrix<T>& operator=(const MatrixInterface<MT>& m) {
        if (reinterpret_cast<const HeapMatrix<T>*>(&m) == this) {
            return *this;
        }

        delete[] data_;

        data_ = new T[m.rowCount() * m.columnCount()];
        rowCount_ = m.rowCount();
        columnCount_ = m.columnCount();
        this->copy(m);

        return *this;
    }
    HeapMatrix<T>& operator=(const HeapMatrix<T>& m) {
        return operator=(Interface(m));
    }
    HeapMatrix<T>& operator=(HeapMatrix<T>&& m) {
        if (&m == this) {
            return *this;
        }
            
        delete[] data_;

        data_ = m.data_;
        rowCount_ = m.rowCount_;
        columnCount_ = m.columnCount_;
        m.data_ = nullptr;

        return *this;
    }
};
#endif
