# matrixLibrary

## Interface macierzy:

Wszystkie macierze dziedziczą po klasie ```MatrixInterface<MatrixType>```, gdzie ```MatrixType``` jest typem macierzy, która dziedziczy po ```MatrixInterface<>```.

To znaczy, że przykładowa deklaracja konkretnej macierzy może wyglądać w ten sposób:
```C++
template<typename T> class SomeMatrix : public MatrixInterface<SomeMatrix<T>>
```
Czyli jest to zastosowanie wzorca projektowego znanego jako "Curiously recurring template pattern": https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern

&nbsp;  

### Wymagania co do dziedziczącej macierzy:
- **Publiczna deklaracja```ValueType``` jako aliasu dla typu trzymanych wartości w macierzy**.

  np.: jeżeli macierz przechowuje wartości typu ```int``` to deklaracja może wyglądać następująco:
  ```C++
  using ValueType = int;
  ```
- **Publiczna deklaracja statycznych stałych typu ```mtl::size_t```: ```Rows```, ```Cols``` i ```EffCols```**.

  Wartości te powinny przedstawiać liczbę wierszy (```Rows```), liczbę kolumn (```Cols```) oraz efektywną liczbę kolumn (```EffCols```).
  
  Wartości te muszą być znane w trakcie kompilacji, a że dana macierz nie musi wymagać, aby jej rozmiar był znany w trakcie kompilacji, to istnieje specjalna statyczna zmienna ```MatrixInterface<>::DynamicValue```, którą można przypisać aby zaznaczyć ten fakt.
  
  Stąd przykładowo dla macierzy, która ma znaną liczbę kolumn, liczbę efektywną kolumn taką samą jak liczbę kolumn oraz nieznaną liczbę wierszy może zadeklarować zmienne w taki sposób:
  ```C++
  static const mtl::size_t Rows = Interface::DynamicValue;
  static const mtl::size_t Cols = 10;
  static const mtl::size_t EffCols = 10;
  ```
  (gdzie ```using Interface = MatrixInterface<TYPE_OF_THIS_MATRIX>;```)
  
- **Istnienie dostępnej dla klasy interface-u zmiennej ```data_``` typu ```ValueType*``` lub funkcji ```ValueType* _data()```**.

  Innymi słowy wymagany jest dostęp do wskaźnika na pierwszy element macierzy.
  
- **Istnienie dostępnej dla klasy interface-u zmiennej ```data_``` typu ```ValueType*``` lub ```const ValueType*``` albo funkcji ```const ValueType* _data() const```**.

  Czyli to samo co wyżej, ale wersja potrzebna w przypadku, gdy macierz jest stała (```const```)
  
- **Istnienie dostępnej dla klasy interface-u zmiennej ```rowCount_``` typu ```mtl::size_t``` lub funkcji ```mtl::size_t _rowCount() const```  lub zadeklarowaną wartość statyczną ```Rows``` o wartości innej niż ```DynamicValue```**.
  
  Czyli istnieje sposób na pobranie liczby wierszy w czasie wykonywania programu.
  
- **Istnienie dostępnej dla klasy interface-u zmiennej ```columnCount_``` typu ```mtl::size_t``` lub funkcji ```mtl::size_t _columnCount() const```  lub zadeklarowaną wartość statyczną ```Cols``` o wartości innej niż ```DynamicValue```**.
  
  Czyli to samo co wyżej, ale dla kolumn.
  
  Dodatkowo to samo można powiedzieć o ```effectiveColumnCount```, ale nie jest to wymagane - w przypadku jeżeli nie ma sposóbu na pobranie wartości ```effectiveColumnCount``` wykorzystywana jest liczba kolumn.

&nbsp;  

### Publiczny interface klasy ```MatrixInterface```:
  ```C++
  static constexpr mtl::size_t ConstexprRowCount();
  static constexpr mtl::size_t ConstexprColumnCount();
  static constexpr mtl::size_t ConstexprEffectiveColumnCount();

  // same as the ones above but with shorter names
  static constexpr mtl::size_t CRow();
  static constexpr mtl::size_t CCol();
  static constexpr mtl::size_t CEffCol();

  static constexpr mtl::size_t DynamicValue;

  static constexpr bool HasConstexprRowAndColumnCount();

  auto createNew() const;
  auto createNew(int rowCount, int columnCount) const;
  template<mtl::size_t RowCount, mtl::size_t ColumnCount> auto  createNew() const;

  auto effectiveColumnCount() const;
  auto columnCount() const;
  auto rowCount() const;
  auto data();
  auto data() const;
  auto operator[](mtl::size_t i);
  auto operator[](mtl::size_t i) const;
  auto& at(mtl::size_t row, mtl::size_t column);
  const auto& at(mtl::size_t row, mtl::size_t column) const;
  auto size() const;
  template<typename MatrixType2> void copy(const MatrixInterface<MatrixType2>& other);
  auto subMatrixView(mtl::size_t startRow, mtl::size_t startColumn, mtl::size_t rowCount, mtl::size_t columnCount);
  auto subMatrixView(mtl::size_t startRow, mtl::size_t startColumn, mtl::size_t rowCount, mtl::size_t columnCount) const;
  template<mtl::size_t RowCount, mtl::size_t ColumnCount> auto subMatrixView(mtl::size_t startRow, mtl::size_t startColumn);
  template<mtl::size_t RowCount, mtl::size_t ColumnCount> auto subMatrixView(mtl::size_t startRow, mtl::size_t startColumn) const;
  auto begin();
  auto end();
  auto begin() const;
  auto end() const;
  ```

  Metody te mogą być przeciążone przez klasy dziedziczące poprzez stworzenie dostępnej dla klasy interface-u metody wewnątrz klasy dziedziczącej o nazwie poprzedzonej znakiem ```_``` (np.: aby przeciążyć metodę ```auto begin()``` trzeba zadeklarować ```T _begin()```, gdzie `T` może być dowolnym typem. Warto zaznaczyć, że odpowiedzialność za dobrą integrację tej funkcji z resztą metod spoczywa na klasie dziedziczącej - np.: przy przeciązeniu ```begin()``` musi ona współpracować z ```end()```).

  W przypadku ```operator[](mtl::size_t)``` metoda przeciążająca musi mieć formę ```_index(mtl::size_t)```.

&nbsp;  

### Stworzone macierze dziedziczące po ```MatrixInterface```:
  #### Przechowujące dane:
  - **```HeapMatrix<typename T>```**

    Macierz o nieznanym w trakcie kompilacji rozmiarze, przechowująca elementy w pamięci alokowanej dynamicznie 

  - **```StaticHeapMatrix<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_>```**

    Macierz o znanym w trakcie kompilacji rozmiarze, przechowująca elementy w pamięci alokowanej dynamicznie 

  - **```StackMatrix<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_>```**

    Macierz o znanym w trakcie kompilacji rozmiarze, przechowująca elementy na stosie 

  #### Widoki:
  - **```MatrixDynamicView<typename T>``` i ```MatrixConstDynamicView<typename T>```**
  - **```MatrixStaticEffectiveColumnCountView<typename T, mtl::size_t EffCols_>``` i ```MatrixConstStaticEffectiveColumnCountView<typename T, mtl::size_t EffCols_>```**
  - **```MatrixStaticSizeView<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_>``` i ```MatrixConstStaticSizeView<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_>```**
  - **```MatrixStaticView<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_, mtl::size_t EffCols_>``` i ```MatrixConstStaticView<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_, mtl::size_t EffCols_>```**

&nbsp;  

### Przykład kodu wykorzystującego macierze:
```C++ 
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
    nonEditable.at(0, 0) = 0; // compile error "you cannot assign to a variable that is const"
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
```
