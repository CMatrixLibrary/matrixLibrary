# matrixLibrary

## Spis treści

- [Interface macierzy](#interface-macierzy)
  * [Wymagania co do dziedziczacej macierzy](#wymagania-co-do-dziedziczacej-macierzy)
  * [Publiczny interface klasy MatrixInterface](#publiczny-interface-klasy-matrixinterface)
  * [Stworzone macierze dziedziczace po MatrixInterface](#stworzone-macierze-dziedziczace-po-matrixinterface)
    + [Przechowujace dane](#przechowujace-dane)
    + [Widoki](#widoki)
  * [Wygodne aliasy dla typow macierzowych](#wygodne-aliasy-dla-typow-macierzowych)
  * [MatrixRowIterator](#matrixrowiterator)
  * [Przyklad kodu wykorzystujacego macierze](#przyklad-kodu-wykorzystujacego-macierze)
- [ThreadPool](#threadpool)
- [Alligned allocation](#alligned-allocation)
- [StackAllocator](#stackallocator)
- [Generyczne bazowe implementacje mnozenia](#generyczne-bazowe-implementacje-mnozenia)
- [AVX i FMA](#avx-i-fma)
- [BLAS](#blas)
- [Generyczne operacje arytmetyczne](#generyczne-operacje-arytmetyczne)
- [Algorytmy szybkiego mnozenia macierzy FMM Fast Matrix Multiplication](#algorytmy-szybkiego-mnozenia-macierzy-fmm-fast-matrix-multiplication)

## Interface macierzy
[MatrixInterface.h](MatrixInterface.h)

Wszystkie macierze dziedziczą po klasie ```MatrixInterface<MatrixType>```, gdzie ```MatrixType``` jest typem macierzy, która dziedziczy po ```MatrixInterface<>```.

To znaczy, że przykładowa deklaracja konkretnej macierzy może wyglądać w ten sposób:
```C++
template<typename T> class SomeMatrix : public MatrixInterface<SomeMatrix<T>>
```
Czyli jest to zastosowanie wzorca projektowego znanego jako "Curiously recurring template pattern": https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern

&nbsp;  

### Wymagania co do dziedziczacej macierzy
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

### Publiczny interface klasy MatrixInterface
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

### Stworzone macierze dziedziczace po MatrixInterface
  #### Przechowujace dane
  - **```HeapMatrix<typename T>```** 
    
    [HeapMatrix.h](HeapMatrix.h)

    Macierz o nieznanym w trakcie kompilacji rozmiarze, przechowująca elementy w pamięci alokowanej dynamicznie 

  - **```StaticHeapMatrix<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_>```** 
    
    [StaticHeapMatrix.h](StaticHeapMatrix.h)

    Macierz o znanym w trakcie kompilacji rozmiarze, przechowująca elementy w pamięci alokowanej dynamicznie 

  - **```StackMatrix<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_>```** 
    
    [StackMatrix.h](StackMatrix.h)

    Macierz o znanym w trakcie kompilacji rozmiarze, przechowująca elementy na stosie 

  #### Widoki
  - **```MatrixDynamicView<typename T>``` i ```MatrixConstDynamicView<typename T>```**
  
    [MatrixDynamicView.h](MatrixDynamicView.h) 
    
    [MatrixConstDynamicView.h](MatrixConstDynamicView.h)
    
  - **```MatrixStaticEffectiveColumnCountView<typename T, mtl::size_t EffCols_>``` i ```MatrixConstStaticEffectiveColumnCountView<typename T, mtl::size_t EffCols_>```**
  
    [MatrixStaticEffectiveColumnCountView.h](MatrixStaticEffectiveColumnCountView.h) 
    
    [MatrixConstStaticEffectiveColumnCountView.h](MatrixConstStaticEffectiveColumnCountView.h)
  
  - **```MatrixStaticSizeView<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_>``` i ```MatrixConstStaticSizeView<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_>```**
  
    [MatrixStaticSizeView.h](MatrixStaticSizeView.h) 
    
    [MatrixConstStaticSizeView.h](MatrixConstStaticSizeView.h)
  
  - **```MatrixStaticView<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_, mtl::size_t EffCols_>``` i ```MatrixConstStaticView<typename T, mtl::size_t RowCount_, mtl::size_t ColumnCount_, mtl::size_t EffCols_>```**

    [MatrixStaticView.h](MatrixStaticView.h) 
    
    [MatrixConstStaticView.h](MatrixConstStaticView.h)

&nbsp;  

### Wygodne aliasy dla typow macierzowych
  - **```Matrix```** 
    
    [Matrix.h](Matrix.h)
    ```C++
    Matrix<typename T, mtl::size_t Rows, mtl::size_t Columns> = 
    - if      (Rows == 0 && Columns == 0)      -> HeapMatrix<T>
    - else if (Rows*Columns <= SOME_THRESHOLD) -> StackMatrix<T, Rows, Columns>
    - else                                     -> StaticHeapMatrix<T, Rows, Columns>
    ```
  - **```MatrixView```** 
    
    [MatrixView.h](MatrixView.h)
    ```C++
    MatrixView<typename T, mtl::size_t Rows, mtl::size_t Columns, mtl::size_t EffColumns> = 
    - if      (Rows != 0 && Columns != 0 && EffColumns != 0) -> MatrixStaticView<T, Rows, Columns, EffColumns>
    - else if (Rows != 0 && Columns != 0)                    -> MatrixStaticSizeView<T, Rows, Columns>
    - else if (Rows != 0)                                    -> MatrixStaticEffectiveColumnView<T, Rows>
    - else                                                   -> MatrixDynamicView<T>
    ```
    
  - **```MatrixConstView```** 
    
    [MatrixView.h](MatrixView.h)
    ```C++
    MatrixConstView<typename T, mtl::size_t Rows, mtl::size_t Columns, mtl::size_t EffColumns> = 
    - if      (Rows != 0 && Columns != 0 && EffColumns != 0) -> MatrixConstStaticView<T, Rows, Columns, EffColumns>
    - else if (Rows != 0 && Columns != 0)                    -> MatrixConstStaticSizeView<T, Rows, Columns>
    - else if (Rows != 0)                                    -> MatrixConstStaticEffectiveColumnView<T, Rows>
    - else                                                   -> MatrixConstDynamicView<T>
    ```

### MatrixRowIterator
[MatrixRowIterator.h](MatrixRowIterator.h); [MatrixConstRowIterator.h](MatrixConstRowIterator.h); [MatrixStaticRowIterator.h](MatrixStaticRowIterator.h); [MatrixConstStaticRowIterator.h](MatrixConstStaticRowIterator.h); 

Iterator po wszystkich rzędach macierzy. Jest typem zwracamym przez metody ```begin()``` i ```end()``` macierzy.

Pozwala to na iteracje po kolejnych rzędach macierzy np. przy użyciu pętli ```for (auto row : matrix)```.

Przy takim zastosowaniu użytkownik nie ma nigdy styczności z tymi typami.

### Przyklad kodu wykorzystujacego macierze
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

## ThreadPool
[ThreadPool.h](ThreadPool.h)

Klasa przeznaczona do przyspieszenia wykonywania danego algorytmu poprzez podzielenie go na podzadania, które są wykonywane równolegle na różnych rdzeniach. Istotne jest, aby nie występowały żadne hazardy (np 2 zadania modyfikujące ten sam obszar pamięci mogą skutkować w błędnym wyniku).

Aby pozwolić na zrównoleglenie, obiekt ThreadPool tworzy i zarządza pulą wątków, które wykonują kolejno zadania (funkcje) dodane przez metodę ```void addTask(std::function<void(void)>)```. Tzn, że jeżeli wszystkie wątki są zajęte, to dodane zadanie czeka w kolejce aż jeden z wątków skończy poprzednio przydzielone zadanie - to oznacza również, że metoda ```addTask``` nie zatrzymuje działania programu.

Aby zatrzymać działanie programu do czasu zakończenia wszystkich zadań, wykorzystywana jest metoda ```void completeTasksAndStop()``` lub ewentualnie automatyczne wywołanie destruktora obiektu ThreadPool.

Przykład użycia:
```C++
void normalVectorSum(int* a, const int* b, const int* c, int size) {
    for (int i = 0; i < size; ++i) {
        a[i] = b[i] + c[i];
    }
}

void parallelVectorSum(int* a, const int* b, const int* c, int size) {
    constexpr int chunkSize = 1'000'000;

    ThreadPool pool;
    for (int i = 0; i < size; i += chunkSize) {
        int partSize = std::min(chunkSize, size - i);
        pool.addTask([=]() { normalVectorSum(a + i, b + i, c + i, partSize); });
    }
}
```

Warto zauważyć, że wzrost szybkosci dzięki paralelizacji jest widoczny dopiero przy bardzo długich zadaniach, a dla krótkich zadań dodatkowy koszt zarządzania wątkami może nawet spowolnić pracę.

## Alligned allocation
[alignedAllocation.h](alignedAllocation.h)

Odpowiedniki funkcji ```malloc```/```free``` i operatorów ```new```, ```new[]```, ```delete```, ```delete[]``` alokujące pamięć o adresach wyrównanych do podanej wartości.

Odpowiedniki operatorów ```new``` i ```new[]``` mają dodatkowo wersję ```Construct```, która w przeciwieństwie do standardowej wersji wymusza inicjalizacje. Powodem tego jest fakt, że operatory ```new``` i ```new[]``` nie inicjalizują obiektów, które nie są klasami, chyba, że doda się nawiasy - ```new()```, ```new[]()```. Odpowiednikiem dodania tych nawiasów jest wywołanie funkcji w wersji z ```Construct``` w nazwie.

Przyklad porównania ```new``` i ```new[]``` dla typów nie klasowych (```Align``` to wybrany alignment): 
``` C++
new int       == alignedNew<int>(Align)
new int()     == alignedNew<int>(Align, 0)  {OR} alignedConstructNew<int>(Align)
new int(123)  == alignedNew<int>(Align, 123) {OR} alignedConstructNew<int>(Align, 123)
new int[10]   == alignedArrayNew<int>(10, Align)
new int[10]() == alignedArrayNew<int>(10, Align, 0) {OR} alignedArrayConstructNew<int>(10, Align)
```
W przypadku alokacji obiektów klas, wersje ```Construct``` i nie-```Construct``` są sobie równe (bo obiekty klas zawsze są inicjalizowane przez konstruktor).

Pamieć zaalokowana przez ```alignedArrayNew``` oraz ```alignedArrayConstructNew``` musi być zwolniona przez ```alignedArrayDelete``` (a nie ```alignedDelete```). W przeciwieństwie do operatora ```delete[]``` poza wskaźnikiem wymagane jest podanie liczby elementów tablicy.

Istnieje również druga wersja alokacji tablicy ```alignedArrayNewWithSize```, która nie wymaga podania rozmiaru przy zwolnieniu za pomocą ```alignedArrayDeleteWithSize```, ale ta wersja jest mniej wydajna, ponieważ wymaga trzymania rozmiaru w alokowanym bloku pamięci.

## StackAllocator
[StackAllocator.h](StackAllocator.h)

Alokator pamięci o kompozycji stosu.

Wszystkie alokacje i dealokacje pamięci mają gwarancję być praktycznie bezpłatne.

W odróżnieniu jednak do typowego alokatora ogólnego przeznaczenia wymaga on, aby zawsze dealokowany był ostatnio zaalokowany blok pamięci (tzn. jeżeli zaalokuje obiekt1 a potem zaalokuje obiekt2 to muszę zdealokować najpierw obiekt2, a potem obiekt1 - dealokacja polega na zdjęciu elementu ze stosu).

```StackAllocator``` poza alokacją zajmuje się również wywołaniem konstruktorów i destruktorów.

## Generyczne bazowe implementacje mnozenia
[naiveMul.h](naiveMul.h); [blockMul.h](blockMul.h); [parallelMul.h](parallelMul.h); [parallelBlockMul.h](parallelBlockMul.h);

- **naiveMul**:

  Najprostsza implementacja (3 pętle for)
  
- **parallelMul**

  To samo co **naiveMul**, ale najbardziej zewnętrzna pętla jest zrównoleglona (wykorzystując ```ThreadPool```)
  
- **blockMul**

  Rozwiącuje problem nieefektywnego wykorzystywania z cache'a (cache thrashing).
  
  Naiwny algorytm oblicza wynik rząd po rzędzie, a żeby obliczyć jeden rząd wynikowy potrzebujemy wczytać każdy element drugiej macierzy. To oznacza, że w momencie gdy skończymy obliczać jeden rząd i przeskakujemy na kolejny, to od nowa zaczynamy czytać elementy z drugiej macierzy, które przy większych macierzach już nie będą w cache'u (bo ostatnie elementy drugiej macierzy zajeły miejsce w cachu pierwszym elementom) - trzeba będzie je znowu wczytać z pamięci.
  
  Rozwiązaniem tego problemu jest operowanie na blokach zamiast na pojedynczych elementach.
  
  Wizualizacja: https://users.ics.aalto.fi/suomela/cache-blocking-demo/ 
  
  Najbardziej na lewo widać naiwny algorytm, a najbardziej na prawo algorytm blokowy.
  - czarny kolor to aktualnie używany element
  - niebieski kolor to L1 cache
  - czerwony kolor to L2 cache
  - ciemno szary kolor to L3 cache
  - jasno szary kolor to pamiec nie w cache.
  Oczywiscie chodzi o to, aby zmaksymalizowac wykorzystanie cache. (animacja jest troche przekłamana, bo każda macierz traktowana jest oddzielnie, tzn każda macierz ma swoją własną niezależną pamięć cache, a w rzeczywistości jest ona współdzielona)

- **parallelBlockMul**

  To samo co **blockMul**, ale najbardziej zewnętrzna pętla jest zrównoleglona (wykorzystując ```ThreadPool```)

## AVX i FMA
[avxSimd.h](avxSimd.h)

Wysokopoziomowy interface ułatwiający wykorzystanie wybranych instrukcji AVX/AVX2/FMA (256 bitowych).

Wspierane typy to:
- 32 i 64 bitowe liczby zmiennoprzecinkowe zgodne ze standardem IEEE 754
- 8,16,32,64 bitowe liczby całkowite (ze znakiem i bez znaku)

256 bitowe rejestry zamodelowane są przez klasę ```AvxType<T>``` (gdzie T jest jednym ze wspieranych typów)

Dostępne funkcje (gdzie T jest jednym ze wspieranych typów):

```C++
template<typename T> constexpr int packedCount(); // liczba elementów mieszczących się w jednym rejestrze
template<typename T> AvxType<T> loadAligned(const T* ptr);
template<typename T> AvxType<T> loadUnaligned(const T* ptr);
template<typename T> void storeAligned(T* dst, AvxType<T> src);
template<typename T> void storeUnaligned(T* dst, AvxType<T> src);
template<typename T> AvxType<T> add(AvxType<T> a, AvxType<T> b);
template<typename T> AvxType<T> sub(AvxType<T> a, AvxType<T> b);
template<typename T> AvxType<T> mul(AvxType<T> a, AvxType<T> b);
template<typename T> AvxType<T> zero();
template<typename T> AvxType<T> setAllElements(T value);
template<typename T> AvxType<T> fma(AvxType<T> a, AvxType<T> b, AvxType<T> c); // a * b + c
```
Dodatkowo wspierane są operatory ```+, -, *, +=, -=, *=```

Interface ten pozwala wygodnie pisać funkcje, które bedą poprawnie działały dla dowolnego wspieranego typu.
Przykład generycznej funkcji, która dodaje do siebie 2 wektory ```a``` i ```b``` zachowując wynik w ```c```:
``` C++
template<typename T> void normalVersion(T* a, T* b, T* c, int size) {
    for (int i = 0; i < size; ++i) {
        c[i] = a[i] + b[i];
    }
}
template<typename T> void avxVersion(T* a, T* b, T* c, int size) {
    int i = 0;
    
    // increment by how many elements fit inside register
    for (; i <= size - avx::packedCount<T>(); i += avx::packedCount<T>()) {
        auto aVector = avx::loadUnaligned(a);
        auto bVector = avx::loadUnaligned(b);
        avx::storeUnaligned(c, aVector + bVector);
    }

    // calculate last few elements when size is not divisible by packedCount<T>.
    for (; i < size; ++i) {
        c[i] = a[i] + b[i];
    }
}
```

## BLAS
[blasMul.h](blasMul.h)

Wrapper na funkcje mnożenia macierzy zgodne z interfacem blas (czyli dla typów ```float```, ```double```, ```std::complex<float>```, ```std::complex<double>```).

Biblioteka sama w sobie nie ma żadnej implementacji interfacu BLAS. Jeżeli użytkownik chce wykorzystać dowolną implementację (np.: MKL czy openBLAS) to musi przed #includem naszej biblioteki:
- #includowac implementację BLAS
- dodać ```#define USE_BLAS```, sygnalizując w ten sposób bibliotece, aby (kiedy to możliwe) wykorzystywała funkcje BLAS do mnożenia macierzy.

To oznacza również, że jeżeli użytkownik nie wykona powyższych kroków, a wywoła bezpośrednio funkcje mnożenia BLAS (tzn jej wrapper należący do naszej biblioteki), to otrzyma błąd kompilacji z odpowiednim komunikatem.

## Generyczne operacje arytmetyczne
[genericArithmeticOperations.h](genericArithmeticOperations.h)

metaprogramowe funkcje do tworzenia wyrażeń arytmetycznych o dowolnej długości. Przykłady:
```C++
calculate<Add>(a)                            >> a              // just return value
calculate<Sub>(a)                            >> -a             // negation
calculate<Add, Add>(a, b)                    >> a + b
calculate<Assign, Add, Add>(a, b, c)         >> a = b + c      
calculate<Assign, Add, Sub, Add>(a, b, c, d) >> a = b - c + d
```

Te pozornie zbędne argumenty ```Add``` są potrzebne, aby było możliwe wszystko zamodelować.

Na przykład musi istnieć możliwość rozróżnienia między ```a = b``` i ```a = -b``` jako odpowiednio ```calculate<Assign, Add>(a, b)``` i ```calculate<Assign, Sub>(a, b)```.

Ze względów wygody zostały jednak dodane skróty, dzięki którym operacje ```Assign``` i ```Add``` dodawane są automatycznie.

Przykłady wykorzystania skróconej notacji:

```C++
calculate<Assign>(a, b) >> a = b
calculate<Add>(a, b, c) >> a = b + c
calculate<Sub>(a, b, c) >> a = b - c
```

Oczywiście te funkcje są same w sobie bezużyteczne, bo zamiast wykorzystywać ```calculate<>``` można poprostu napisać dane wyrażenie.

Funkcja ```calculate<>``` jest za to użyteczna, gdy wykorzystamy ją w większej funkcji. Na przykład można stworzyć funkcję, która wykona operacje arytmetyczne na wszystkich elementach 3 wektorów:

```C++
template<OpType op1, OpType op2, OpType op3> operateOnVectors(int* a, int* b, int* C, int size) {
    for (int i = 0; i < size; ++i) {
        calculate<op1, op2, op3>(a[i], b[i], c[i]);
    }
}
```
Teraz możemy wywołać tą funkcję z dowolną permutacją operacji. 

Np. jeżeli chcemy dodać do siebie wektory ```b``` i ```c``` i zapisać wynik w ```a``` to wywołujemy funkcję ```operateOnVectors<Assign, Add, Add>(...)```. 

Albo jeżeli chcemy skopiować wektor ```c``` do ```a``` i ```b``` (a\[i\] = b\[i\] = c\[i\]) To możemy wywołać funkcję ```operateOnVectors<Assign, Assign, Add>(...)```

## Algorytmy szybkiego mnozenia macierzy FMM Fast Matrix Multiplication 
[fmmUtility.h](fmmUtility.h)

[strassen.h](strassen.h)

Mnożenie macierzy polegające na rekurencyjnym dzieleniu macierzy na mniejsze bloki, które umiemy pomnożyć wykorzystując mniejszą liczbę operacji niż w przypadku zwykłego mnożenia (z mniejszą złożonością obliczeniową niż ```O(n^3)```).

Ze względu na znacznie większe stałe (w sensie złożoności obliczeniowej) występujące w implementacji fmm w porównaniu do zwykłego podejścia, najlepiej jest zależnie od typu danych i rozmiarów macierzy "przerzucić się" w pewnym momencie na zwykłe mnożenie. W implementacji jest to zamodelowane przez liczbę kroków rekurencyjnych (```steps```), które chcemy wykonać przed wywołaniem zwkłego algorytmu mnożenia.

Algorytmy fmm można zaimplementować na wiele sposobów i mamy kilka różnych implementacji:
- Podejście **High-Level**:
  
  Najbardziej naiwna implementacja wykorzystująca wysokopoziomowy interface macierzy. 
  
  Jej głównymi problemami jest bardzo dużo dynamicznych alokacji nowych macierzy oraz wykonywanie zbyt dużej liczby kopii macierzy. W przypadku niewielu kroków rekurencyjnych nie ma to większego znaczenia, ale im więcej kroków wykonujemy tym bardziej te problemy są widoczne w szybkości wykonania.
  
- Podejście **Low-Level**:

  Implementacja Low-Level poza wykorzystywaniem surowych wskaźników + int-ów do trzymania rozmiarów (które same w sobie nie dają praktycznie żadnej przewagi w porównaniu z **High-Level**) korzysta z faktu, że znając na start liczbę kroków rekurencyjnych oraz rozmiary macierzy możemy dokładnie obliczyć ile potrzeba zarezerwować pamięci, aby wystarczyło na zapis wszystkich pośrednich macierzy. 
  
  Oznacza to, że mamy tylko jedną alokację pamięci i wielokrotnie wykorzystujemy ją do różnych celów (do implementacji tego wykorzystany jest ```StackAllocator<T>```).
  
  Przykład obliczenia rozmiaru wymaganego dla strassena ```C = A * B``` (mnożenie o rozmiarach ```<N, M, P>``` czyli ```Cnp = Anm * Bmp```):
  - Na każdy krok rekurencyjny potrzebujemy 1 macierz wynikową o rozmiarze ```(N/2) * (P/2)```. Mamy 7 kroków rekurencyjnych stąd mamy ```7 * (N/2)*(P/2)```. 
  - Potrzebujemy pomocnicze macierze do trzymanie wyniku dodawania/odejmowania pod-macierzy A i B. Stąd mamy ```(N/2)*(M/2) + (M/2)*(P/2)```.
  
  Łącznie wychodzi ```(N/2)*(M/2) + (M/2)*(P/2) + 7*(N/2)*(P/2) = (NM + MP + 7NP)/4```.
  
  Nie jest to jeszcze ostateczny wynik. To jest dopiero łączna suma pamięci potrzebna dla danego kroku rekurencyjnego, a my wkonujemy ```steps``` kroków rekurencyjnych.
  
  Ostateczne równanie wygląda więc następująco:
  
  ```
                      { (NM + MP + 7NP)/4 + f(N/2, M/2, P/2, steps - 1); dla steps > 1
  f(N, M, P, steps) = { (NM + MP + 7NP)/4;                               dla steps = 1
                      { 0;                                               dla steps < 1
  ```
  
  Aby zrozumieć skalę wymaganej pamięci najlepiej popatrzeć na mnożenie macierzy kwadratowych (```N = M = P```). Dla kolejnych kroków rekurencyjnych mamy wtedy wymagane rozmiary:
  ```
  - 1 krok:  2.25     * N^2
  - 2 kroki: 2.8125   * N^2
  - 3 kroki: 2.953125 * N^2
  - limit:   3        * N^2
  ```
  W takim wypadku macierze ```A```, ```B``` i ```C``` mają łączny rozmiar ```3 * N^2```, czyli wymagana dodatkowa pamięć operacyjna w porównaniu do pamięci na macierze wejściowe wynosi aż ```100%```.
  
  Poza alokacją pamięci, aby zwiększyć szybkość wykonywania, wszystkie dłuższe ciągi obliczeń arytmetycznych, np. jak obliczenie sumy macierzy ```r = x + y + z``` są wykonywane bezpośrednio, zamiast wykonywać ```temp = x + y; r = temp + z```
  
- Podejście **Min-Space**:
  
  Implementacja bazująca na **Low-Level**, która zmniejsza wymaganą dodatkową pamięć do minimum.
  
  W tej implementacji zamiast trzymać wszystkie macierze wynikowe z wywołań rekurencyjnych (np.: w przypadku strassena trzymanie 7 macierzy), po każdym wywołaniu rekurencyjnym odpowiednie podmacierze wynikowe wpisują/dodawają/odejmują wynikową macierz z wywołania rekurencyjnego. 
  
  Dzięki temu potrzebujemy jedynie jednej macierzy trzymającej pośredni wynik z wywołania rekurencyjnego.
  
  Funkcja rozmiaru wygląda więc następująco (dla strassena):
  ```
                      { (NM + MP + NP)/4 + f(N/2, M/2, P/2, steps - 1); dla steps > 1
  f(N, M, P, steps) = { (NM + MP + NP)/4;                               dla steps = 1
                      { 0;                                              dla steps < 1
  ```
  A więc rozmiary dla kolejnych kroków (dla strassena):
  ```
  - 1 krok:  0.75     * N^2
  - 2 kroki: 0.9375   * N^2
  - 3 kroki: 0.984375 * N^2
  - limit:   1        * N^2
  ```
  Czyli wymagana dodatkowa pamięć operacyjna w porównaniu do pamięci na macierze wejściowe wynosi ```33.(3)%```.

  W przeciwieństwie do **Low-Level** ze względu na to, że na bierząco wykorzystujemy wyniki rekurencyjnych wywołań, nie ma możliwości wykonywania dłuższych ciągów arytmetycznych naraz. W przypadku macierzy o nieznanym rozmiarze w trakcie kompilacji to praktycznie nie ma wpływu na szybkość wykonywania, jednak w przypadku macierzy o znanym rozmiarze w trakcie kompilacji różnica jest znacząca.

- Podejście **Parallel Low-Level**:
  
  Poprzednie podejścia jedynie pozwalały na zrównoleglanie mnożenia bazowego (dla ```steps = 0```).
  
  Zarówno **Low-Level** jak i **Min-Space** nie dają bezpośredniej możliwości zrównoleglenia. Jest to spowodowane tym, że wielokrotnie wykorzystują one ten sam obszar pamięci dla kolejnych wywołań rekurencyjnych. 
  
  Aby umożliwić zrównoleglenie kilku wywołań rekurencyjnych naraz, musimy oddzielnie dla każdego z nich zaalokować pamięć. Sprowadza się to do tego, że nie możemy mieć jedynie po 1 macierzy pomocniczej do trzymanie wyniku dodawania/odejmowania pod-macierzy A i B. Potrzebujemy ich tyle, ile jest kroków rekurencyjnych (lub ewentualnie przynajmniej tyle, ile mamy wątków).
  
  Znowu patrząc na strassena dla porównania, przy alokacji wszystkich dodatkowych pomocniczych macierzy mamy rozmiar ```4.25 * N^2``` dla pierwszegeo kroku. 
  
  Kolejne kroki rekurencyjne są wykonywane wykorzystując **Min-Space** (możnaby również wykorzystać **Low-Level**, ale narazie w implementacji wykorzystywany jest **Min-Space**).
  
  Przy takich założeniach pełna funkcja rozmiaru dla strassena wygląda następująco (```f```):
  ```
  f(N, M, P, steps) = steps > 0 ? 5NM + 5MP + 7NP + g(N, M, P, steps - 1) : 0
  
                      { (NM + MP + NP)/4 + f(N/2, M/2, P/2, steps - 1); dla steps > 1
  g(N, M, P, steps) = { (NM + MP + NP)/4;                               dla steps = 1
                      { 0;                                              dla steps < 1
  ```
  A więc rozmiary dla kolejnych kroków mamy:
  ```
  - 1 krok:  4.25     * N^2
  - 2 kroki: 4.4375   * N^2
  - 3 kroki: 4.484375 * N^2
  - limit:   4.5      * N^2
  ```
  Czyli wymagana dodatkowa pamięć operacyjna w porównaniu do pamięci na macierze wejściowe wynosi aż ```150%```.
  
  Wymagany rozmiar jest więc spory, zwłaszcza w porównaniu do **Min-Space**, jednak jest to najszybsza implementacja dla dużych macierzy.

W przypadku każdej z tych implementacji możemy jeszcze dodatkowo wprowadzić modyfikację, w której przed samym bazowym mnożeniem sprawdzamy czy efektywna liczba kolumn macierzy równa się ich liczbie kolumn i w przypadku, jeżeli nie, to skopiować je do tymczasowych zmiennych. Dzięki temu algorytm bazowego mnożenia może potencjalnie działać troche szybciej co poprawi szybkość całego mnożenia. Jest to zamodelowane przez wybór ```BaseMulSize::Normal``` lub ```BaseMulSize::Effective```.

W algorytmach fmm występuje również problem w przypadku, gdy rozmiary macierzy nie są podzielne przez dane dla danego algorytmu wartości (np.: w przypadku strassena mamy podział na 2). Istnieją w tym wypadku 4 rozwiązania:
  - static padding
  - dynamic padding
  - static peeling
  - dynamic peeling
  
W naszej bibliotece zaimplementowane są static padding i dynamic peeling, zamodelowane jako wybór ```ResizeStrategy::StaticPadding``` lub ```ResizeStrategy::DynamicPeeling```.
