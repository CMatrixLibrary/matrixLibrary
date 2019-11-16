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

struct OpMatrix {
    std::string letter;
    std::string operation;
    int row;
    int col;
    double coeff;
    
    OpMatrix(const std::string& letter, const std::string& operation, int row, int col, double coeff) : letter(letter), operation(operation), row(row), col(col), coeff(coeff) {}

    std::string name() const {
        return "d" + letter + "[" + std::to_string(row) + "][" + std::to_string(col) + "]";
    }
};
struct TempMatrix {
    std::string letter;
    std::string specialName;
    std::vector<OpMatrix> opMatrices;
    
    bool isView() const {
        return !(opMatrices.size() > 1 || opMatrices[0].operation == "Sub" || hasCoeff());
    }
    bool hasCoeff() const {
        for (int i = 0; i < opMatrices.size(); ++i) {
            if (opMatrices[i].coeff != 1) {
                return true;
            }
        }
        return false;
    }
    std::string eff() const {
        if (isView()) return "eff" + letter;
        else          return (letter == "A") ? "dm" : "dp";
    }
    std::string name() const {
        if (isView())                  return opMatrices[0].name();
        else if (!specialName.empty()) return specialName + letter;
        else                           return "temp" + letter;
    }
};
struct ABMatrix {
    TempMatrix a;
    TempMatrix b;
};
void setTemporaryMatrix(const std::string& letter, const Matrix<double>& matrix, int n, int m, int mulIndex, TempMatrix& tempMatrix) {
    tempMatrix.letter = letter;
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < m; ++col) {
            double value = matrix.at(col, row, mulIndex);
            if (value > 0) {
                tempMatrix.opMatrices.emplace_back(letter, "Add", row, col, std::abs(value));
            } else if (value < 0) {
                tempMatrix.opMatrices.emplace_back(letter, "Sub", row, col, std::abs(value));
            }
        }
    }
}
std::vector<ABMatrix> getTemporaryABMatrices(const Input& input) {
    std::vector<ABMatrix> result(input.mulCount);
    for (int i = 0; i < input.mulCount; ++i) {
        setTemporaryMatrix("A", input.a, input.xSize, input.ySize, i, result[i].a);
        setTemporaryMatrix("B", input.b, input.ySize, input.zSize, i, result[i].b);
    }
    return result;
}

void printCommonAlgorithmStart(std::ofstream& out, const std::string& algorithmRecursiveName, const Input& input) {
    out << "namespace fmm::detail {\n";
    out << "    struct " << algorithmRecursiveName << " {\n";
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
}
void printCommonAlgorithmEnd(std::ofstream& out) {
    out << "        }\n";
    out << "    };\n";
    out << "}\n";
    out << '\n';
}

void printOperationEff(const std::string& indent, std::ofstream& out, int n, int m, int mulIndex, const Matrix<double>& matrix, const TempMatrix& tempMatrix) {
    if (!tempMatrix.isView()) {
        bool hasCoeff = tempMatrix.hasCoeff();
        out << indent;
        if (hasCoeff) out << "operationEffWithCoeffs<Assign";
        else          out << "operationEff<Assign";
        for (auto& opMatrix : tempMatrix.opMatrices) {
            out << ", " << opMatrix.operation;
        }
        if (tempMatrix.letter == "A") out << ">(dn, dm, dm, effA, " << tempMatrix.name();
        if (tempMatrix.letter == "B") out << ">(dm, dp, dp, effB, " << tempMatrix.name();
        for (auto& opMatrix : tempMatrix.opMatrices) {
            if (hasCoeff) out << ", std::pair(" << opMatrix.coeff << ", " << opMatrix.name() << ")";
            else          out << ", " << opMatrix.name();
        }
        out << ");\n";
    }
}

std::vector<std::pair<std::string, std::string>> 
getOperationsOnFirstArgValues(std::string matrixName, int n, int m, int mulIndex, const Matrix<double>& matrix, std::vector<bool>& assignedResultMatrices) {
    std::vector<std::pair<std::string, std::string>> result;
    bool hasCoeff = false;
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < m; ++col) {
            auto value = matrix.at(col, row, mulIndex);
            if (std::abs(value) != 1 && value != 0) {
                hasCoeff = true;
            }
            if (!assignedResultMatrices[col + row * m] && value == -1) {
                hasCoeff = true;
            }
        }
    }
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < m; ++col) {
            auto value = matrix.at(col, row, mulIndex);
            if (value != 0) {
                if (hasCoeff) {
                    if (!assignedResultMatrices[col + row * m]) {
                        result.emplace_back("Assign", "std::pair(" + std::to_string(value) + ", " + matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "])");
                    } else if (value > 0) {
                        result.emplace_back("AddAssign", "std::pair(" + std::to_string(std::abs(value)) + ", " + matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "])");
                    } else if (value < 0) {
                        result.emplace_back("SubAssign", "std::pair(" + std::to_string(std::abs(value)) + ", " + matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "])");
                    }
                } else {
                    if (!assignedResultMatrices[col + row * m]) {
                        result.emplace_back("Assign", matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
                    } else if (value > 0) {
                        result.emplace_back("AddAssign", matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
                    } else if (value < 0) {
                        result.emplace_back("SubAssign", matrixName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
                    }
                }
                assignedResultMatrices[col + row * m] = true;
            }
        }
    }
    return result;
}
void generateMinSpaceAlgorithm(std::ofstream& out, const std::string& algorithmName, const std::string& algorithmRecursiveName, const Input& input, const std::vector<ABMatrix>& abMatrices) {
    printCommonAlgorithmStart(out, algorithmRecursiveName, input);

    // temporary matrix allocation
    out << "            auto tempA = allocator.alloc(dn * dm);\n";
    out << "            auto tempB = allocator.alloc(dm * dp);\n";
    out << "            auto tempC = allocator.alloc(dn * dp);\n\n";

    // recursive calls
    std::vector<bool> assignedResultMatrices(input.xSize*input.zSize, false);
    for (int i = 0; i < input.mulCount; ++i) {
        printOperationEff("            ", out, input.xSize, input.ySize, i, input.a, abMatrices[i].a);
        printOperationEff("            ", out, input.ySize, input.zSize, i, input.b, abMatrices[i].b);
        out << "            nextStep<Method, BaseN, BaseM, BaseP, " << algorithmRecursiveName << ">(tempC, "
            << abMatrices[i].a.name() << ", " << abMatrices[i].b.name() << ", dn, dm, dp, dp, "
            << abMatrices[i].a.eff() << ", " << abMatrices[i].b.eff() << ", steps - 1, allocator);\n";
        auto operationValues = getOperationsOnFirstArgValues("dC", input.xSize, input.zSize, i, input.c, assignedResultMatrices);
        if (operationValues[0].second[0] == 's') out << "            operationsOnFirstArgWithCoeffs<" << operationValues[0].first;
        else out << "            operationsOnFirstArg<" << operationValues[0].first;
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

    // temporary matrix deallocation
    out << "            allocator.dealloc(tempC, dn*dp);\n";
    out << "            allocator.dealloc(tempB, dm*dp);\n";
    out << "            allocator.dealloc(tempA, dn*dm);\n";

    printCommonAlgorithmEnd(out);
}

void printOperationEffCMatrix(std::ofstream& out, const Input& input) {
    bool hasCoeff = false;
    for (int row = 0; row < input.xSize; ++row) {
        for (int col = 0; col < input.zSize; ++col) {
            std::vector<std::pair<std::string, std::string>> usedMatrices;

            bool hasCoeff = false;
            for (int mulIndex = 0; mulIndex < input.mulCount; ++mulIndex) {
                auto value = input.c.at(col, row, mulIndex);
                if (std::abs(value) != 1 && value != 0) {
                    hasCoeff = true;
                }
            }

            for (int mulIndex = 0; mulIndex < input.mulCount; ++mulIndex) {
                auto value = input.c.at(col, row, mulIndex);
                if (hasCoeff) {
                    if (value > 0) {
                        usedMatrices.emplace_back("Add", "std::pair(" + std::to_string(std::abs(value)) + ", m" + std::to_string(mulIndex + 1) + ")");
                    } else if (value < 0) {
                        usedMatrices.emplace_back("Sub", "std::pair(" + std::to_string(std::abs(value)) + ", m" + std::to_string(mulIndex + 1) + ")");
                    }
                } else {
                    if (value > 0) {
                        usedMatrices.emplace_back("Add", "m" + std::to_string(mulIndex + 1));
                    } else if (value < 0) {
                        usedMatrices.emplace_back("Sub", "m" + std::to_string(mulIndex + 1));
                    }
                }
            }

            if (hasCoeff) out << "            operationEffWithCoeffs<Assign";
            else          out << "            operationEff<Assign";
            for (auto& m : usedMatrices) {
                out << ", " << m.first;
            }
            out << ">(dn, dp, effC, dp, dC[" << std::to_string(row) << "][" << std::to_string(col) << "]";
            for (auto& m : usedMatrices) {
                out << ", " << m.second;
            }
            out << ");\n";
        }
    }
    out << "\n";
}
void generateLowLevelAlgorithm(std::ofstream& out, const std::string& algorithmName, const std::string& algorithmRecursiveName, const Input& input, std::vector<ABMatrix>& abMatrices) {
    printCommonAlgorithmStart(out, algorithmName, input);

    // recursive calls
    for (int i = 0; i < input.mulCount; ++i) {
        out << "            auto m" << i+1 << " = allocator.alloc(dn * dp);\n";

        auto& abMatrix = abMatrices[i];
        if (!abMatrix.a.isView()) {
            abMatrix.a.specialName = "m" + std::to_string(i+1);
            out << "            auto " << abMatrix.a.name() << " = allocator.alloc(dn*dm);\n";
        }
        if (!abMatrix.b.isView()) {
            abMatrix.b.specialName = "m" + std::to_string(i+1);
            out << "            auto " << abMatrix.b.name() << " = allocator.alloc(dm*dp);\n";
        }

        printOperationEff("            ", out, input.xSize, input.ySize, i, input.a, abMatrices[i].a);
        printOperationEff("            ", out, input.ySize, input.zSize, i, input.b, abMatrices[i].b);
        out << "            nextStep<Method, BaseN, BaseM, BaseP, " << algorithmRecursiveName << ">(m"
            << std::to_string(i+1) << ", "
            << abMatrices[i].a.name() << ", " << abMatrices[i].b.name() << ", dn, dm, dp, dp, "
            << abMatrices[i].a.eff() << ", " << abMatrices[i].b.eff() << ", steps - 1, allocator);\n";

        if (!abMatrix.b.isView()) {
            out << "            allocator.dealloc(" << abMatrix.b.name() << ", dm*dp);\n";
        }
        if (!abMatrix.a.isView()) {
            out << "            allocator.dealloc(" << abMatrix.a.name() << ", dn*dm);\n";
        }
        out << "\n";
    }

    printOperationEffCMatrix(out, input);

    // deallocation of result matrices
    for (int i = input.mulCount; i >= 1; --i) {
        out << "            allocator.dealloc(m" << i << ", dn*dp);\n";
    }
    printCommonAlgorithmEnd(out);
}

void generateLowLevelParallelAlgorithm(std::ofstream& out, const std::string& algorithmName, const std::string& algorithmRecursiveName, const Input& input, std::vector<ABMatrix>& abMatrices) {
    printCommonAlgorithmStart(out, algorithmName, input);

    // allocation of result matrices
    for (int i = 1; i <= input.mulCount; ++i) {
        out << "            auto m" << i << " = allocator.alloc(dn * dp);\n";
    } 

    // allocation of temporary A and B matrices
    for (int i = 0; i < abMatrices.size(); ++i) {
        auto& abMatrix = abMatrices[i];
        if (!abMatrix.a.isView()) {
            abMatrix.a.specialName = "m" + std::to_string(i+1);
            out << "            auto " << abMatrix.a.name() << " = allocator.alloc(dn*dm);\n";
        }
        if (!abMatrix.b.isView()) {
            abMatrix.b.specialName = "m" + std::to_string(i+1);
            out << "            auto " << abMatrix.b.name() << " = allocator.alloc(dm*dp);\n";
        }
    }

    // recursive calls
    out << "\n            ThreadPool pool;\n\n";
    for (int i = 0; i < input.mulCount; ++i) {
        out << "            pool.addTask([=, dn = dn, dm = dm, dp = dp]() {\n";
        printOperationEff("                ", out, input.xSize, input.ySize, i, input.a, abMatrices[i].a);
        printOperationEff("                ", out, input.ySize, input.zSize, i, input.b, abMatrices[i].b);
        out << "                minSpaceRun<Method, " << input.xSize << ", " << input.ySize << ", " << input.zSize << ", "
            << input.mulCount << ", " << algorithmRecursiveName << ">(m" << std::to_string(i + 1) << ", "
            << abMatrices[i].a.name() << ", " << abMatrices[i].b.name() << ", dn, dm, dp, dp, "
            << abMatrices[i].a.eff() << ", " << abMatrices[i].b.eff() << ", steps - 1);\n";
        out << "            });\n";
    }
    out << "\n            pool.completeTasksAndStop();\n\n";

    // saving result of recursive calls to C matrix
    printOperationEffCMatrix(out, input);

    // deallocation of temporary A and B matrices
    for (int i = abMatrices.size() - 1; i >= 0; --i) {
        auto& abMatrix = abMatrices[i];
        if (!abMatrix.b.isView()) {
            out << "            allocator.dealloc(" << abMatrix.b.name() << ", dm*dp);\n";
        }
        if (!abMatrix.a.isView()) {
            out << "            allocator.dealloc(" << abMatrix.a.name() << ", dn*dm);\n";
        }
    }

    // deallocation of result matrices
    for (int i = input.mulCount; i >= 1; --i) {
        out << "            allocator.dealloc(m" << i << ", dn*dp);\n";
    }

    printCommonAlgorithmEnd(out);
}
void generateCall(std::ofstream& out, const std::string& algorithmName, const std::string& algorithmRecursiveName, const std::string& postfix, const Input& input, int tempACount=0, int tempBCount=0) {
    out << "    template<int Method = 0, typename M1, typename M2>\n";
    out << "    auto "
        << algorithmName << postfix << "(const MatrixInterface<M1>& a, const MatrixInterface<M2>& b, int steps) {\n";
    out << "        return detail::runAlgorithm<getNewWithAlgorithm<Method, Algorithm::" << postfix << ">, "
        << input.xSize << ", " << input.ySize << ", " << input.zSize << ", " << input.mulCount << ", "
        << "detail::" << algorithmRecursiveName << postfix << ">(a, b, steps, " 
        << std::to_string(tempACount) << ", " << std::to_string(tempBCount) << ");\n";
    out << "    }\n";
}

void generateAlgorithms(const fs::path& inputPath, const fs::path& outputPath) {    
    std::ifstream inputFile(inputPath);
    auto inputOpt = parseInputFile(inputFile);
    if (!inputOpt) {
        std::cerr << inputPath << ": wrong file format\n";
        return;
    }
    auto input = *inputOpt;
    auto algorithmName = inputPath.stem().string();
    auto algorithmRecursiveName = algorithmName + "Recursive";
    algorithmRecursiveName[0] = toupper(algorithmRecursiveName[0]);

    auto temporaryABMatrices = getTemporaryABMatrices(input);

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
    out << "#include \"../ThreadPool.h\"\n";
    out << '\n';

    generateMinSpaceAlgorithm(out, algorithmRecursiveName + "MinSpace", algorithmRecursiveName + "MinSpace", input, temporaryABMatrices);
    generateLowLevelAlgorithm(out, algorithmRecursiveName + "LowLevel", algorithmRecursiveName + "LowLevel", input, temporaryABMatrices);
    generateLowLevelParallelAlgorithm(out, algorithmRecursiveName + "LowLevelParallel", algorithmRecursiveName + "MinSpace", input, temporaryABMatrices);

    out << "namespace fmm {\n";
    generateCall(out, algorithmName, algorithmRecursiveName, "MinSpace", input);
    generateCall(out, algorithmName, algorithmRecursiveName, "LowLevel", input);
    int tempACount = 0;
    int tempBCount = 0;
    for (auto& abMatrix : temporaryABMatrices) {
        tempACount += !abMatrix.a.isView();
        tempBCount += !abMatrix.b.isView();
    }
    generateCall(out, algorithmName, algorithmRecursiveName, "LowLevelParallel", input, tempACount, tempBCount);
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
            generateAlgorithms(dirEntry.path(), algorithmsFolder / fs::path(dirEntry.path().stem().string() + ".h"));
        }
    }

    return 0;
}
