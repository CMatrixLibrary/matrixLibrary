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

std::pair<std::vector<std::pair<std::string, std::string>>, int> 
getOperationEffValues(std::string matrixName, int n, int m, int mulIndex, const Matrix<double>& matrix) {
    std::vector<std::pair<std::string, std::string>> result;
    int coeff = 0;
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < m; ++col) {
            auto value = matrix.at(col, row, mulIndex);
            if (value != 0 && value != 1 && value != -1) {
                coeff = value;
            }
            if (value > 0) {
                result.emplace_back("Add", matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
            } else if (value < 0) {
                result.emplace_back("Sub", matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
            }
        }
    }
    return { result, coeff };
}
std::string printOperationEff(std::ofstream& out, std::string matrixLetter, int n, int m, int mulIndex, const Matrix<double>& matrix) {
    auto[operationValues, coeff] = getOperationEffValues("d" + matrixLetter, n, m, mulIndex, matrix);
    bool performOperation = coeff || operationValues.size() > 1 || operationValues[0].first == "Sub";
    if (performOperation) {
        if (coeff) out << "            operationEffWithCoeff<Assign";
        else       out << "            operationEff<Assign";
        for (auto& operationValue : operationValues) {
            out << ", " << operationValue.first;
        }
        out << ">(dn, dm, dm, eff" << matrixLetter << ", temp" << matrixLetter;
        if (coeff) out << ", " << coeff;
        for (auto& operationValue : operationValues) {
            out << ", " << operationValue.second;
        }
        out << ");\n";
        return "";
    } else {
        return operationValues[0].second;
    }
}
std::vector<std::pair<std::string, std::string>> 
getOperationsOnFirstArgValues(std::string matrixName, int n, int m, int mulIndex, const Matrix<double>& matrix, std::vector<bool>& assignedResultMatrices) {
    std::vector<std::pair<std::string, std::string>> result;
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < m; ++col) {
            auto value = matrix.at(col, row, mulIndex);
            if (value != 0) {
                if (!assignedResultMatrices[col + row * m]) {
                    result.emplace_back("Assign", matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
                } else if (value > 0) {
                    result.emplace_back("AddAssign", matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
                } else if (value < 0) {
                    result.emplace_back("SubAssign", matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
                }
                assignedResultMatrices[col + row * m] = true;
            }
        }
    }
    return result;
}

void printRecursiveCalls(std::ofstream& out, const Input& input, const std::string& functionName) {
    std::vector<bool> assignedResultMatrices(input.xSize*input.zSize, false);

    for (int i = 0; i < input.mulCount; ++i) {
        std::string aMatrix = printOperationEff(out, "A", input.xSize, input.ySize, i, input.a);
        std::string bMatrix = printOperationEff(out, "B", input.ySize, input.zSize, i, input.b);
        out << "            nextStep<Method, BaseN, BaseM, BaseP, " << functionName << ">(tempC, ";
        if (aMatrix.empty()) out << "tempA";
        else                 out << aMatrix;
        out << ", ";
        if (bMatrix.empty()) out << "tempB";
        else                 out << bMatrix;
        out << ", dn, dm, dp, dp";
        if (aMatrix.empty()) out << ", dm";
        else                 out << ", effA";
        if (bMatrix.empty()) out << ", dp";
        else                 out << ", effB";
        out << ", steps - 1, allocator);\n";
        auto operationValues = getOperationsOnFirstArgValues("dC", input.xSize, input.zSize, i, input.c, assignedResultMatrices);
        out << "            operationsOnFirstArg<" << operationValues[0].first;
        for (int i = 1; i < operationValues.size(); ++i) {
            out << ", " << operationValues[i].first;
        }
        out << ">(dn, dp, dp, effC, tempC";
        for (auto& operationValue : operationValues) {
            out << ", " << operationValue.second;
        }
        out << ");\n";
        out << "\n";
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
    auto recursiveFunctionName = algorithmName + "Recursive";
    recursiveFunctionName[0] = toupper(recursiveFunctionName[0]);

    std::ofstream out(outputPath);
    out << "// this file was generated using fastMatrixMultiplyAlgorithms/generator\n";

    std::string defineName = "GEN_FMM_"
        + std::to_string(input.xSize) + "_"
        + std::to_string(input.ySize) + "_"
        + std::to_string(input.zSize) + "_"
        + std::to_string(input.mulCount) + "_H";

    out << "#ifndef " << defineName << '\n';
    out << "#define " << defineName << '\n';
    out << "#include \"../fmmUtility.h\"\n";
    out << '\n';
    out << "namespace fmm::detail {\n";
    out << "    struct " << recursiveFunctionName << " {\n";
    out << "        template<int Method, typename T>\n";
    out << "        static void Run(T* c, T* a, T* b, int n, int m, int p, int effC, int effA, int effB, int steps, StackAllocator<T>& allocator) {\n";
    out << "            using namespace ArithmeticOperation;\n";
    out << "\n";
    out << "            constexpr int BaseN = " << input.xSize << ";\n";
    out << "            constexpr int BaseM = " << input.ySize << ";\n";
    out << "            constexpr int BaseP = " << input.zSize << ";\n";
    out << "\n";
    out << "            auto[dn, dm, dp] = divideSizes<BaseN, BaseM, BaseP>(n, m, p);\n";
    out << "\n";
    out << "            auto dA = divideView<BaseN, BaseM>(a, n, m, effA);\n";
    out << "            auto dB = divideView<BaseM, BaseP>(b, m, p, effB);\n";
    out << "            auto dC = divideView<BaseN, BaseP>(c, n, p, effC);\n";
    out << "\n";
    out << "            auto tempA = allocator.alloc(dn * dm);\n";
    out << "            auto tempB = allocator.alloc(dm * dp);\n";
    out << "            auto tempC = allocator.alloc(dn * dp);\n";
    out << "\n";
    printRecursiveCalls(out, input, recursiveFunctionName);
    out << "            allocator.dealloc(tempC, dn*dp);\n";
    out << "            allocator.dealloc(tempB, dm*dp);\n";
    out << "            allocator.dealloc(tempA, dn*dm);\n";
    out << "        }\n";
    out << "    };\n";
    out << "}\n";
    out << '\n';
    out << "namespace fmm {\n";
    out << "    template<int Method = 0, typename M1, typename M2>\n";
    out << "    auto "
        << algorithmName
        << "MinSpace(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {\n";
    out << "        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::MinSpace>, "
        << input.xSize << ", " << input.ySize << ", " << input.zSize << ", " << input.mulCount << ", " 
        << "detail::" << recursiveFunctionName << ">(a, b, steps);\n";
    out << "    }\n";
    out << "}\n";
    out << '\n';
    out << "#endif\n";
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
