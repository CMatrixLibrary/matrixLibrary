#ifndef CACHE_PADDING_HEAP_MATRIX_H
#define CACHE_PADDING_HEAP_MATRIX_H
#include "MatrixInterface.h"
#include <cstdlib>

template<typename T> class CachePaddedHeapMatrix : public MatrixInterface<CachePaddedHeapMatrix<T>> {
    using Interface = MatrixInterface<CachePaddedHeapMatrix<T>>;
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
    mtl::size_t effectiveColumnCount_;

public:
    CachePaddedHeapMatrix() :
        data_(nullptr),
        rowCount_(0),
        columnCount_(0),
        effectiveColumnCount_(0)
    {}
    CachePaddedHeapMatrix(mtl::size_t rowCount, mtl::size_t columnCount) :
        rowCount_(rowCount),
        columnCount_(columnCount),
        effectiveColumnCount_(columnCount)
    {
        if (columnCount_ % 64 == 0) {
            effectiveColumnCount_ += 1;
        }
        arrayNew();
    }
    template<typename MT> CachePaddedHeapMatrix(const MatrixInterface<MT>& m) :
        columnCount_(m.columnCount()),
        rowCount_(m.rowCount()),
        effectiveColumnCount_(m.effectiveColumnCount())
    {
        arrayNew();
        this->copy(m);
    }
    template<typename MT> CachePaddedHeapMatrix(const MatrixInterface<MT>& m, mtl::size_t rCount, mtl::size_t cCount) :
        columnCount_(cCount),
        rowCount_(rCount),
        effectiveColumnCount_(cCount)
    {
        if (columnCount_ % 64 == 0) {
            effectiveColumnCount_ += 1;
        }
        arrayNew();
        this->copy(m);
    }
    CachePaddedHeapMatrix(const CachePaddedHeapMatrix<T>& m) :
        columnCount_(m.columnCount()),
        rowCount_(m.rowCount()),
        effectiveColumnCount_(m.effectiveColumnCount())
    {
        arrayNew();
        this->copy(m);
    }
    CachePaddedHeapMatrix(CachePaddedHeapMatrix<T>&& m) :
        data_(m.data_),
        columnCount_(m.columnCount()),
        rowCount_(m.rowCount()),
        effectiveColumnCount_(m.effectiveColumnCount())
    {
        m.data_ = nullptr;
    }
    ~CachePaddedHeapMatrix() {
        arrayDelete();
    }

    template<typename MT> CachePaddedHeapMatrix<T>& operator=(const MatrixInterface<MT>& m) {
        if (reinterpret_cast<const CachePaddedHeapMatrix<T>*>(&m) == this) {
            return *this;
        }

        arrayDelete();

        rowCount_ = m.rowCount();
        columnCount_ = m.columnCount();
        effectiveColumnCount_ = m.effectiveColumnCount();
        arrayNew();

        this->copy(m);

        return *this;
    }
    CachePaddedHeapMatrix<T>& operator=(const CachePaddedHeapMatrix<T>& m) {
        return operator=(Interface(m));
    }
    CachePaddedHeapMatrix<T>& operator=(CachePaddedHeapMatrix<T>&& m) {
        if (&m == this) {
            return *this;
        }
            
        arrayDelete();

        data_ = m.data_;
        rowCount_ = m.rowCount_;
        columnCount_ = m.columnCount_;
        effectiveColumnCount_ = m.effectiveColumnCount_;
        m.data_ = nullptr;

        return *this;
    }

    CachePaddedHeapMatrix _createNew() const {
        return CachePaddedHeapMatrix<T>(rowCount_, columnCount_);
    }
    CachePaddedHeapMatrix _createNew(mtl::size_t rowCount, mtl::size_t columnCount) const {
        return CachePaddedHeapMatrix<T>(rowCount, columnCount);
    }

    mtl::size_t _effectiveColumnCount() const { return effectiveColumnCount_; }

private:
    void arrayNew() {
        data_ = static_cast<T*>(std::malloc(rowCount_ * effectiveColumnCount_ * sizeof(T)));
        if constexpr (std::is_class_v<T>) {
            for (std::size_t row = 0; row < rowCount_; ++row) {
                for (std::size_t column = 0; column < columnCount_; ++column) {
                    new(&data_[column + row*effectiveColumnCount_]) T();
                }
            }
        }
    }
    void arrayDelete() {
        if constexpr (std::is_class_v<T>) {
            for (std::size_t row = 0; row < rowCount_; ++row) {
                for (std::size_t column = 0; column < columnCount_; ++column) {
                    data_[column + row * effectiveColumnCount_].~T();
                }
            }
        }
        std::free(data_);
    }
};
#endif
