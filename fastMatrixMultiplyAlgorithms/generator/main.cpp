#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>
#include <tuple>
#include <optional>

namespace fs = std::filesystem;

template<typename T> struct Matrix {
    Matrix() {}
    Matrix(int rowCount, int colCount) :
        rowCount(rowCount),
        colCount(colCount),
        data(rowCount * colCount)
    {}

    T& operator[](int index) {
        return data[index];
    }
    const T& operator[](int index) const {
        return data[index];
    }

    T& at(int column, int row) {
        return data[column + row * colCount];
    }
    const T& at(int column, int row) const {
        return data[column + row * colCount];
    }

    int rowCount;
    int colCount;
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
template<typename T> std::optional<Matrix<T>> readMatrix(std::istream& in, int rowCount, int colCount) {
    Matrix<T> matrix(rowCount, colCount);
    for (int i = 0; i < rowCount*colCount; ++i) {
        auto value = read<T>(in);
        if (!value) {
            return std::nullopt;
        }
        matrix[i] = *value;
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

    auto a = readMatrix<double>(inputFile, input.xSize * input.ySize, input.mulCount);
    auto b = readMatrix<double>(inputFile, input.ySize * input.zSize, input.mulCount);
    auto c = readMatrix<double>(inputFile, input.xSize * input.zSize, input.mulCount);

    if (a && b && c) {
        input.a = *a;
        input.b = *b;
        input.c = *c;
        return input;
    } else {
        return std::nullopt;
    }
}

void printRecursiveCallMultiplyFactor(std::ofstream& out, int column, int matrixColSize, const std::string& matrixName, const Matrix<double>& m) {
    bool isFirstElement = true;
    for (int row = 0; row < m.rowCount; ++row) {
        if (m.at(column, row) != 0) {
            if (isFirstElement && m.at(column, row) < 0) {
                out << "-1*";
            } else if (!isFirstElement) {
                out << " " << (m.at(column, row) > 0 ? '+' : '-') << " ";
            }
            isFirstElement = false;

            if (abs(m.at(column, row)) != 1) {
                out << abs(m.at(column, row)) << "*";
            }
            out << matrixName << "[" << row / matrixColSize << "][" << row % matrixColSize << "]";
        }
    }
}
void printRecursiveCalls(std::ofstream& out, const Input& input, const std::string& algorithmName) {
    out << "        std::array<FullMatrix<T>, " << input.mulCount << "> m;\n";
    for (int i = 0; i < input.mulCount; ++i) {
        out << "        m[" << i << "] = " << algorithmName << "("; 
        printRecursiveCallMultiplyFactor(out, i, input.ySize, "a", input.a);
        out << ", ";
        printRecursiveCallMultiplyFactor(out, i, input.zSize, "b", input.b);
        out << ", steps - 1);\n";
    }
}

void printCMatrixConstruction(std::ofstream& out, const Input& input) {
    for (int row = 0; row < input.xSize; ++row) {
        for (int column = 0; column < input.zSize; ++column) {
            out << "        c[" << row << "][" << column << "].copy(";
            bool isFirstElement = true;
            for (int i = 0; i < input.mulCount; ++i) {
                auto elementValue = input.c.at(i, column + row * input.zSize);
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
                    out << "m[" << i << "]";
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
