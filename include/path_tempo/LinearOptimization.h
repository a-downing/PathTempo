#pragma once

#include <cstddef>
#include <expected>
#include <initializer_list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace path_tempo {
    inline constexpr double SMALL_MATRIX_VALUE = 1e-12;

    constexpr bool retainsMatrixCoefficient(const double value) {
        return value > SMALL_MATRIX_VALUE || value < -SMALL_MATRIX_VALUE;
    }

    double linearProgramInfinity() noexcept;

    class SparseLinearProgram {
        struct Implementation;
        std::unique_ptr<Implementation> m_implementation;

    public:
        explicit SparseLinearProgram(std::size_t columns);
        ~SparseLinearProgram();
        SparseLinearProgram(SparseLinearProgram &&) noexcept;
        SparseLinearProgram &operator=(SparseLinearProgram &&) noexcept;
        SparseLinearProgram(const SparseLinearProgram &) = delete;
        SparseLinearProgram &operator=(const SparseLinearProgram &) = delete;

        void addRow(double lower, double upper, std::initializer_list<std::pair<std::size_t, double>> entries);

        [[nodiscard]] std::size_t columnCount() const noexcept;
        [[nodiscard]] std::size_t rowCount() const noexcept;
        [[nodiscard]] std::size_t nonzeroCount() const noexcept;
        [[nodiscard]] double &columnCost(std::size_t column);
        [[nodiscard]] double &columnLower(std::size_t column);
        [[nodiscard]] double &columnUpper(std::size_t column);
        [[nodiscard]] double rowLower(std::size_t row) const;
        [[nodiscard]] double rowUpper(std::size_t row) const;
        [[nodiscard]] std::size_t rowBegin(std::size_t row) const;
        [[nodiscard]] std::size_t rowEnd(std::size_t row) const;
        [[nodiscard]] std::size_t entryColumn(std::size_t entry) const;
        [[nodiscard]] double entryValue(std::size_t entry) const;

    private:
        friend class PersistentLinearSolver;
    };

    enum class LinearSolveStatus {
        Optimal,
        TimeLimit,
        IterationLimit,
    };

    struct LinearSolveSettings {
        std::size_t simplexIterationLimit = 0;
        bool reuseBasis = true;
    };

    struct LinearSolveDiagnostics {
        std::size_t simplexIterations = 0;
        bool basisReuseAttempted = false;
        bool basisReuseApplied = false;
        bool basisDimensionMismatch = false;
        bool modelUpdateAttempted = false;
        bool modelUpdateApplied = false;
        bool modelStructureMismatch = false;
    };

    struct LinearSolveResult {
        LinearSolveStatus status = LinearSolveStatus::Optimal;
        std::vector<double> values;
        LinearSolveDiagnostics diagnostics;
    };

    class PersistentLinearSolver {
        struct Implementation;
        std::unique_ptr<Implementation> m_implementation;

    public:
        PersistentLinearSolver();
        ~PersistentLinearSolver();
        PersistentLinearSolver(PersistentLinearSolver &&) noexcept;
        PersistentLinearSolver &operator=(PersistentLinearSolver &&) noexcept;
        PersistentLinearSolver(const PersistentLinearSolver &) = delete;
        PersistentLinearSolver &operator=(const PersistentLinearSolver &) = delete;

        std::expected<void, std::string> configure(double timeLimit);
        std::expected<LinearSolveResult, std::string> solve(const SparseLinearProgram &program, const LinearSolveSettings &settings);
    };
}
