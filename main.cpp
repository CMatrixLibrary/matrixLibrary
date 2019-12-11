#include <iostream>
#if __has_include("mkl.h")
    #include "mkl.h"
    #define USE_BLAS
#endif
#include "matrixLibrary.h"

void example() {
    HeapMatrix<int> dh1(5, 5);        // HeapMatrix<int> = Matrix<int>
    StaticHeapMatrix<int, 5, 5> sh1;  // 
    StackMatrix<int, 5, 5> ss1;       // StackMatrix<int, 5, 5> = Matrix<int, 5, 5>

    int initValue = 0;

    // iterate over all rows of each matrix
    for (auto[rowDh, rowSh, rowSS] : rangeZip(dh1, sh1, ss1)) {

        // iterate over all values of each matrix
        for (auto[valueDh, valueSh, valueSs] : rangeZip(rowDh, rowSh, rowSS)) { 
            valueDh = initValue;
            valueSh = initValue;
            valueSs = initValue;
            initValue += 1;
        }
    }

    for (auto rowSh : sh1) { // iterate over all rows of 'sh1' matrix
        for (auto index : indicies(rowSh)) { // iterate over all indexes in rowSh
            rowSh[index] *= 2;
        }
    }

    // traditional indexing using .at(row, column)
    for (int i = 0; i < ss1.rowCount(); ++i) {
        for (int j = 0; j < ss1.columnCount(); ++j) {
            ss1.at(i, j) = ss1.size() - i*j;
            // ss1[i][j] = ... would be logically the same 
            // (with good compiler code should also be the same)
        }
    }

    // transposition
    dh1 = transpose(dh1);
    sh1 = transpose(sh1);
    ss1 = transpose(ss1);

    // subMatrixView with 2 rows and 5 columns starting from 0th row and 0th column
    // best possible type (one with most compile time information for optimization) is automaticly chosen
    auto v1 = dh1.subMatrixView(0, 0, 2, 5); // type(v1) = MatrixDynamicView<int>
    auto v2 = sh1.subMatrixView(0, 0, 2, 5); // type(v2) = MatrixStaticEffectiveColumnView<int, 5>
    auto v3 = ss1.subMatrixView<2, 5>(0, 0); // type(v3) = MatrixStaticView<int, 2, 2, 5>

    // overloaded +, -, +=, -=, *, *=, << operators 
    // (note how different matrix types seamlessly work with each other)
    v1 *= 10;
    v1 += v2;
    dh1 = v1 + v2 + v3;
    sh1 -= ss1;
    sh1 = sh1 - ss1;
    auto c = dh1 * sh1;
    std::cout << c << '\n';
}



// ------------------------------------------------------
// passing matrices to functions: 
// ------------------------------------------------------

// function which can be with any matrix type arguments.
// In case of MatrixView<> the matrix argument has to be editable so non-const specified
// and cannot be other const view.
// The only problem with this function is it will generate sub-optimal code if argument
// had compile time known size of matrix which could be used to optimize the code by compiler
void simpleGenericFunction(MatrixView<int> editable, MatrixConstView<int> nonEditable) {
    editable.at(0, 0) = nonEditable.at(0, 0);
    //nonEditable.at(0, 0) = 0; // compile error "you cannot assign to a variable that is const"
}


// same as above but with generic type
template<typename T> void simpleGenericFunction2(MatrixView<T> matrix) {
    std::cout << matrix.data()[0] << '\n'; // directly access pointer to data with .data()
}


// compared to the previous functions this one besides working with any matrix type
// will also use compile time constant size if such matrix has them provided.
// This is similar to dynamic polimorphism using pointer to base class, but unlike
// traditional class polymorphism this one will evaluate everything staticly (at compile-time)
template<typename MatrixType> 
void fullyGenericFunction(const MatrixInterface<MatrixType>& matrix) {
    // to access value type of matrix you can use "typename MatrixType::ValueType"
    typename MatrixType::ValueType value = matrix[0][0];
}

void exampleFunctionCalls() {
    HeapMatrix<int>                   a(3, 5);
    StaticHeapMatrix<int, 7, 1>       b;
    StackMatrix<double, 1, 2>         c;
    MatrixDynamicView<int>            d = a;
    MatrixStaticView<double, 1, 2, 2> e = c;

    simpleGenericFunction(a, a);
    simpleGenericFunction(a, b);
    simpleGenericFunction(a, d);
    simpleGenericFunction(d, b);

    // unfortunatelly in case of using generic type MatrixView<> you need 
    // to provide the type explicitly on call
    simpleGenericFunction2<int>(a);
    simpleGenericFunction2<int>(b);
    simpleGenericFunction2<double>(c);
    simpleGenericFunction2<int>(d);
    simpleGenericFunction2<double>(e);

    fullyGenericFunction(a);
    fullyGenericFunction(b);
    fullyGenericFunction(c);
    fullyGenericFunction(d);
    fullyGenericFunction(e);
}

int main() {
    return 0;
}

