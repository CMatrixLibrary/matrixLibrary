#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>
#include <tuple>
#include <optional>
#include <map>
#include <numeric>

namespace fs = std::filesystem;

template<typename T> struct Matrix {
    Matrix() {}
    Matrix(int rowCount, int colCount, int mulCount) :
        rowCount(rowCount),
        colCount(colCount),
        mulCount(mulCount),
        data(rowCount * colCount * mulCount)
    {}

    T& operator[](int index) {
        return data[index];
    }
    const T& operator[](int index) const {
        return data[index];
    }

    T& at(int column, int row, int mul) {
        return data[column + row * colCount + mul * rowCount * colCount];
    }
    const T& at(int column, int row, int mul) const {
        return data[column + row * colCount + mul * rowCount * colCount];
    }

    int rowCount;
    int colCount;
    int mulCount;
    std::vector<T> data;
};

void skipLine(std::istream& in) {
    std::string s;
    std::getline(in, s);
}
std::optional<std::string> readString(std::istream& in) {
    std::string str;
    while (in >> str) {
        if (str[0] == '#') {
            skipLine(in);
        } else {
            break;
        }
    }
    if (in) return str;
    else    return std::nullopt;
}
template<typename T> std::optional<T> read(std::istream& in) {
    static_assert(false, "no implementation for read<T>(std::istream&) for given type");
}
template<> std::optional<int> read(std::istream& in) {
    if (auto str = readString(in); str) {
        return std::stoi(*str);
    } else {
        return std::nullopt;
    }
}
template<> std::optional<double> read(std::istream& in) {
    if (auto str = readString(in); str) {
        return std::stod(*str);
    } else {
        return std::nullopt;
    }
}
template<typename T> std::optional<Matrix<T>> readMatrix(std::istream& in, int rowCount, int colCount, int mulCount) {
    Matrix<T> matrix(rowCount, colCount, mulCount);
    for (int row = 0; row < rowCount; ++row) {
        for (int col = 0; col < colCount; ++col) {
            for (int mul = 0; mul < mulCount; ++mul) {
                auto value = read<T>(in);
                if (!value) {
                    return std::nullopt;
                }
                matrix.at(col, row, mul) = *value;
            }
        }
    }
    return matrix;
}

struct Input {
    int xSize;
    int ySize;
    int zSize;
    int mulCount;
    Matrix<double> a;
    Matrix<double> b;
    Matrix<double> c;
};

std::optional<Input> parseInputFile(std::ifstream& inputFile) {
    Input input;
    
    auto xOpt = read<int>(inputFile);
    auto yOpt = read<int>(inputFile);
    auto zOpt = read<int>(inputFile);
    auto mOpt = read<int>(inputFile);

    if (!xOpt || !yOpt || !zOpt || !mOpt) {
        return std::nullopt;
    }
    input.xSize = *xOpt;
    input.ySize = *yOpt;
    input.zSize = *zOpt;
    input.mulCount = *mOpt;

    auto a = readMatrix<double>(inputFile, input.xSize, input.ySize, input.mulCount);
    auto b = readMatrix<double>(inputFile, input.ySize, input.zSize, input.mulCount);
    auto c = readMatrix<double>(inputFile, input.xSize, input.zSize, input.mulCount);

    if (a && b && c) {
        input.a = *a;
        input.b = *b;
        input.c = *c;
        return input;
    } else {
        return std::nullopt;
    }
}

struct MatrixValue {
    MatrixValue(int firstIndex, int secondIndex, double coeff) :
        firstIndex(firstIndex),
        secondIndex(secondIndex),
        coeff(coeff)
    {}
    MatrixValue(int firstIndex, int secondIndex) :
        firstIndex(firstIndex),
        secondIndex(secondIndex),
        coeff(0)
    {}

    std::string getName(const std::string& matrixName) {
        if (coeff == 0) {
            return matrixName + "_" + std::to_string(firstIndex) + "_" + std::to_string(secondIndex);
        } else {
            return matrixName + "[" + std::to_string(firstIndex) + "][" + std::to_string(secondIndex) + "]";
        }
    }

    bool isBaseValue() {
        return coeff != 0;
    }

    int firstIndex;
    int secondIndex;
    double coeff;
};

void printRecursiveCallMultiplyFactor(
    std::ofstream& out, 
    int mul, 
    const std::string& matrixName, 
    const Matrix<double>& m,
    std::vector<MatrixValue> matrixValues,
    std::vector<std::vector<int>> matrixIndexes
) {
    bool isFirstElement = true;
    for (int index : matrixIndexes[mul]) {
        if (isFirstElement && matrixValues[index].coeff < 0) {
            out << "-";
        } else if (!isFirstElement) {
            out << " " << (matrixValues[index].coeff > 0 ? '+' : '-') << " ";
        }
        isFirstElement = false;
        if (abs(matrixValues[index].coeff) != 1 && matrixValues[index].coeff != 0) {
            out << abs(matrixValues[index].coeff) << "*";
        }
        out << matrixValues[index].getName(matrixName);
    }
}
std::pair<std::vector<MatrixValue>, std::vector<std::vector<int>>> createMatrixIndexes(std::ofstream& out, const std::string& matrixName, const Matrix<double>& m) {
    std::vector<MatrixValue> matrixValues;
    std::vector<std::vector<int>> matrixIndexes(m.mulCount);
    
    for (int mul = 0; mul < m.mulCount; ++mul) {
        for (int row = 0; row < m.rowCount; ++row) {
            for (int col = 0; col < m.colCount; ++col) {
                auto coeff = m.at(col, row, mul);
                if (coeff != 0) {
                    auto foundValue = std::find_if(matrixValues.begin(), matrixValues.end(), [&](auto& value) {
                        return value.firstIndex == row && value.secondIndex == col && value.coeff == coeff;
                    });
                    if (foundValue == matrixValues.end()) {
                        matrixValues.emplace_back(row, col, coeff);
                        matrixIndexes[mul].emplace_back(matrixValues.size() - 1);
                    } else {
                        matrixIndexes[mul].emplace_back(std::distance(matrixValues.begin(), foundValue));
                    }
                }
            }
        }
    }

    while (true) {
        std::vector<int> indexCounts(matrixValues.size() * matrixValues.size(), 0);
        for (int i = 0; i < matrixIndexes.size(); ++i) {
            for (int j = 0; j < matrixIndexes[i].size(); ++j) {
                for (int k = j + 1; k < matrixIndexes[i].size(); ++k) {
                    indexCounts[matrixIndexes[i][j] + matrixIndexes[i][k] * matrixValues.size()] += 1;
                }
            }
        }
        std::vector<int> indexCountsIndexes(indexCounts.size());
        std::iota(indexCountsIndexes.begin(), indexCountsIndexes.end(), 0);
        std::sort(indexCountsIndexes.begin(), indexCountsIndexes.end(), [&indexCounts](auto i1, auto i2) {
            return indexCounts[i1] < indexCounts[i2];
        });
        auto maxCount = indexCounts[indexCountsIndexes.back()];
        if (maxCount > 1) {
            int firstValueIndex = indexCountsIndexes.back() % matrixValues.size();
            int secondValueIndex = indexCountsIndexes.back() / matrixValues.size();
            matrixValues.emplace_back(firstValueIndex, secondValueIndex);
            for (int i = 0; i < matrixIndexes.size(); ++i) {
                for (int j = 0; j < matrixIndexes[i].size(); ++j) {
                    for (int k = j + 1; k < matrixIndexes[i].size(); ++k) {
                        if (firstValueIndex == matrixIndexes[i][j] && secondValueIndex == matrixIndexes[i][k]) {
                            matrixIndexes[i].erase(matrixIndexes[i].begin() + k);
                            matrixIndexes[i].erase(matrixIndexes[i].begin() + j);
                            matrixIndexes[i].emplace_back(matrixValues.size() - 1);
                        }
                    }
                }
            }
        } else {
            break;
        }
    }

    for (auto& value : matrixValues) {
        if (!value.isBaseValue()) {
            out << "        auto ";
            out << value.getName(matrixName);
            out << " = ";
            if (matrixValues[value.firstIndex].coeff > 0) {
                out << matrixValues[value.firstIndex].getName(matrixName);
                out << (matrixValues[value.secondIndex].coeff > 0 ? " + " : " - ");
                if (abs(matrixValues[value.secondIndex].coeff) != 1 && matrixValues[value.secondIndex].coeff != 0) {
                    out << abs(matrixValues[value.secondIndex].coeff) << "*";
                }
                out << matrixValues[value.secondIndex].getName(matrixName);
            } else if (matrixValues[value.secondIndex].coeff > 0) {
                out << matrixValues[value.secondIndex].getName(matrixName);
                out << (matrixValues[value.firstIndex].coeff > 0 ? " + " : " - ");
                if (abs(matrixValues[value.firstIndex].coeff) != 1 && matrixValues[value.firstIndex].coeff != 0) {
                    out << abs(matrixValues[value.firstIndex].coeff) << "*";
                }
                out << matrixValues[value.firstIndex].getName(matrixName);
            } else {
                out << "-";
                if (abs(matrixValues[value.firstIndex].coeff) != 1 && matrixValues[value.firstIndex].coeff != 0) {
                    out << abs(matrixValues[value.firstIndex].coeff) << "*";
                }
                out << matrixValues[value.firstIndex].getName(matrixName);
                out << " - ";
                if (abs(matrixValues[value.secondIndex].coeff) != 1 && matrixValues[value.secondIndex].coeff != 0) {
                    out << abs(matrixValues[value.secondIndex].coeff) << "*";
                }
                out << matrixValues[value.secondIndex].getName(matrixName);
            }
            
            out << ";\n";
        }
    }

    return std::pair(matrixValues, matrixIndexes);
}
void printRecursiveCalls(std::ofstream& out, const Input& input, const std::string& algorithmName) {
    out << "        std::array<FullMatrix<T>, " << input.mulCount << "> m;\n";
    auto[aMatrixValues, aMatrixIndexes] = createMatrixIndexes(out, "a", input.a);
    auto[bMatrixValues, bMatrixIndexes] = createMatrixIndexes(out, "b", input.b);

    for (int i = 0; i < input.mulCount; ++i) {
        out << "        m[" << i << "] = " << algorithmName << "("; 
        printRecursiveCallMultiplyFactor(out, i, "a", input.a, aMatrixValues, aMatrixIndexes);
        out << ", ";
        printRecursiveCallMultiplyFactor(out, i, "b", input.b, bMatrixValues, bMatrixIndexes);
        out << ", steps - 1);\n";
    }
}

void printCMatrixConstruction(std::ofstream& out, const Input& input) {
    auto& c = input.c;
    for (int row = 0; row < c.rowCount; ++row) {
        for (int col = 0; col < c.colCount; ++col) {
            out << "        c[" << row << "][" << col << "].copy(";
            bool isFirstElement = true;
            for (int mul = 0; mul < c.mulCount; ++mul) {
                auto elementValue = c.at(col, row, mul);
                if (elementValue != 0) {
                    if (isFirstElement && elementValue < 0) {
                        out << "-1*";
                    }
                    else if (!isFirstElement) {
                        out << " " << (elementValue > 0 ? '+' : '-') << " ";
                    }
                    isFirstElement = false;

                    if (abs(elementValue) != 1) {
                        out << abs(elementValue) << "*";
                    }
                    out << "m[" << mul << "]";
                }
            }
            out << ");\n";
        }
    }
}

void generateAlgorithm(const fs::path& inputPath, const fs::path& outputPath) {    
    std::ifstream inputFile(inputPath);
    auto inputOpt = parseInputFile(inputFile);
    if (!inputOpt) {
        std::cerr << inputPath << ": wrong file format\n";
        return;
    }
    auto input = *inputOpt;
    auto algorithmName = inputPath.stem().string();

    std::ofstream out(outputPath);
    out << "// this file was generated using fastMatrixMultiplyAlgorithms/generator\n";
    out << "#pragma once\n";
    out << "#include <array>\n";
    out << "#include \"../FullMatrix.h\"\n";
    out << "#include \"../FullMatrixView.h\"\n";
    out << "#include \"../FullMatrixConstView.h\"\n";
    out << "#include \"../naiveOperations.h\"\n";
    out << "#include \"../operators.h\"\n";
    out << "#include \"../utilityDetails.h\"\n";
    out << '\n';
    out << "namespace details {\n";
    out << "    template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>\n";
    out << "    FullMatrix<T> " << algorithmName << "(const MatrixA<T>& A, const MatrixB<T>& B, int steps) {\n";
    out << "        if (steps <= 0) {\n";
    out << "            return naiveMul(A, B);\n";
    out << "        }\n";
    out << '\n';
    out << "        auto a = matrixDivide<" << input.xSize << ", " << input.ySize << ">(A);\n";
    out << "        auto b = matrixDivide<" << input.ySize << ", " << input.zSize << ">(B);\n";
    out << '\n';
    printRecursiveCalls(out, input, algorithmName);
    out << '\n';
    out << "        FullMatrix<T> C(A.rowCount(), B.columnCount());\n";
    out << "        auto c = matrixDivide<" << input.xSize << ", " << input.zSize << ">(C);\n";
    out << '\n';
    printCMatrixConstruction(out, input);
    out << '\n';
    out << "        return C;\n";
    out << "    }\n";
    out << "}\n";
    out << '\n';
    out << "template<typename T, template<typename> typename MatrixA, template<typename> typename MatrixB>\n";
    out << "FullMatrix<T> " << algorithmName << "(const MatrixA<T>& A, const MatrixB<T>& B) {\n";
    out << "    return details::fastMul<" << input.xSize << ", " << input.ySize << ", " << input.zSize << ">(A, B, 50, details::" << algorithmName << "<T, MatrixA, MatrixB>);\n";
    out << "}\n";
}

int main(int argc, char** argv) {
    if (argc <= 1) {
        std::cerr << "missing required first argument: input directory path\n";
        return 1;
    }
    if (argc <= 2) {
        std::cerr << "missing required second argument: output directory path\n";
        return 2;
    }

    fs::path coefficientsFolder = argv[1];
    if (!fs::exists(coefficientsFolder)) {
        std::cerr << "provided input directory path does not exist";
        return 3;
    }

    fs::path algorithmsFolder = argv[2];
    if (!fs::exists(coefficientsFolder)) {
        std::cerr << "provided output directory path does not exist";
        return 4;
    }

    for (auto& dirEntry : fs::recursive_directory_iterator(coefficientsFolder)) {
        if (dirEntry.is_regular_file()) {
            generateAlgorithm(dirEntry.path(), algorithmsFolder / fs::path(dirEntry.path().stem().string() + ".h"));
        }
    }

    return 0;
}
