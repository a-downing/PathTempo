#include "path_tempo/LinearOptimization.h"

#include <cassert>
#include <cmath>
#include <format>
#include <limits>
#include <optional>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#pragma clang diagnostic ignored "-Wunused-parameter"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include <highs/Highs.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace path_tempo {
    struct SparseLinearProgram::Implementation {
        HighsLp model;
    };

    struct PersistentLinearSolver::Implementation {
        Highs solver;
        std::optional<HighsLp> model;
        std::optional<HighsBasis> basis;
        bool configured = false;
    };

    double linearProgramInfinity() noexcept {
        return kHighsInf;
    }

    SparseLinearProgram::SparseLinearProgram(const std::size_t columns) : m_implementation(std::make_unique<Implementation>()) {
        auto &model = m_implementation->model;
        model.num_col_ = static_cast<HighsInt>(columns);
        model.num_row_ = 0;
        model.sense_ = ObjSense::kMinimize;
        model.offset_ = 0.0;
        model.col_cost_.assign(columns, 0.0);
        model.col_lower_.assign(columns, -kHighsInf);
        model.col_upper_.assign(columns, kHighsInf);
        model.a_matrix_.format_ = MatrixFormat::kRowwise;
        model.a_matrix_.start_.assign(1, 0);
    }

    SparseLinearProgram::~SparseLinearProgram() = default;

    SparseLinearProgram::SparseLinearProgram(SparseLinearProgram &&) noexcept = default;

    SparseLinearProgram &SparseLinearProgram::operator=(SparseLinearProgram &&) noexcept = default;

    void SparseLinearProgram::addRow(const double lower, const double upper, const std::initializer_list<std::pair<std::size_t, double>> entries) {
        auto &model = m_implementation->model;

        for (const auto &[column, value] : entries) {
            assert(std::isfinite(value));
            assert(column < static_cast<std::size_t>(model.num_col_));

            if (!retainsMatrixCoefficient(value)) {
                continue;
            }

            model.a_matrix_.index_.push_back(static_cast<HighsInt>(column));
            model.a_matrix_.value_.push_back(value);
        }

        model.a_matrix_.start_.push_back(static_cast<HighsInt>(model.a_matrix_.index_.size()));
        model.row_lower_.push_back(lower);
        model.row_upper_.push_back(upper);
        ++model.num_row_;
    }

    std::size_t SparseLinearProgram::columnCount() const noexcept {
        return static_cast<std::size_t>(m_implementation->model.num_col_);
    }

    std::size_t SparseLinearProgram::rowCount() const noexcept {
        return static_cast<std::size_t>(m_implementation->model.num_row_);
    }

    std::size_t SparseLinearProgram::nonzeroCount() const noexcept {
        return m_implementation->model.a_matrix_.index_.size();
    }

    double &SparseLinearProgram::columnCost(const std::size_t column) {
        assert(column < columnCount());

        return m_implementation->model.col_cost_[column];
    }

    double &SparseLinearProgram::columnLower(const std::size_t column) {
        assert(column < columnCount());

        return m_implementation->model.col_lower_[column];
    }

    double &SparseLinearProgram::columnUpper(const std::size_t column) {
        assert(column < columnCount());

        return m_implementation->model.col_upper_[column];
    }

    double SparseLinearProgram::rowLower(const std::size_t row) const {
        assert(row < rowCount());

        return m_implementation->model.row_lower_[row];
    }

    double SparseLinearProgram::rowUpper(const std::size_t row) const {
        assert(row < rowCount());

        return m_implementation->model.row_upper_[row];
    }

    std::size_t SparseLinearProgram::rowBegin(const std::size_t row) const {
        assert(row < rowCount());

        return static_cast<std::size_t>(m_implementation->model.a_matrix_.start_[row]);
    }

    std::size_t SparseLinearProgram::rowEnd(const std::size_t row) const {
        assert(row < rowCount());

        return static_cast<std::size_t>(m_implementation->model.a_matrix_.start_[row + 1]);
    }

    std::size_t SparseLinearProgram::entryColumn(const std::size_t entry) const {
        assert(entry < nonzeroCount());

        return static_cast<std::size_t>(m_implementation->model.a_matrix_.index_[entry]);
    }

    double SparseLinearProgram::entryValue(const std::size_t entry) const {
        assert(entry < nonzeroCount());

        return m_implementation->model.a_matrix_.value_[entry];
    }

    PersistentLinearSolver::PersistentLinearSolver() : m_implementation(std::make_unique<Implementation>()) {
    }

    PersistentLinearSolver::~PersistentLinearSolver() = default;

    PersistentLinearSolver::PersistentLinearSolver(PersistentLinearSolver &&) noexcept = default;

    PersistentLinearSolver &PersistentLinearSolver::operator=(PersistentLinearSolver &&) noexcept = default;

    std::expected<void, std::string> PersistentLinearSolver::configure(const double timeLimit) {
        if (!std::isfinite(timeLimit) || timeLimit <= 0.0) {
            return std::unexpected("linear solver time limit must be finite and positive");
        }

        auto &solver = m_implementation->solver;

        if (solver.setOptionValue("output_flag", false) != HighsStatus::kOk || solver.setOptionValue("threads", HighsInt {1}) != HighsStatus::kOk || solver.setOptionValue("solver", std::string {"simplex"}) != HighsStatus::kOk || solver.setOptionValue("small_matrix_value", SMALL_MATRIX_VALUE) != HighsStatus::kOk || solver.setOptionValue("time_limit", timeLimit) != HighsStatus::kOk) {
            return std::unexpected("could not configure HiGHS");
        }

        m_implementation->configured = true;

        return {};
    }

    std::expected<LinearSolveResult, std::string> PersistentLinearSolver::solve(const SparseLinearProgram &program, const LinearSolveSettings &settings) {
        if (!m_implementation->configured) {
            return std::unexpected("linear solver has not been configured");
        }

        if (settings.simplexIterationLimit == 0 || settings.simplexIterationLimit > static_cast<std::size_t>(std::numeric_limits<HighsInt>::max())) {
            return std::unexpected("linear solver simplex iteration limit is invalid");
        }

        const auto &model = program.m_implementation->model;
        auto &solver = m_implementation->solver;
        LinearSolveResult result;

        if (solver.setOptionValue("simplex_iteration_limit", static_cast<HighsInt>(settings.simplexIterationLimit)) != HighsStatus::kOk) {
            return std::unexpected("could not configure the HiGHS simplex iteration limit");
        }

        const auto sameStructure = settings.reuseBasis && m_implementation->model && m_implementation->model->num_col_ == model.num_col_ && m_implementation->model->num_row_ == model.num_row_ && m_implementation->model->a_matrix_.format_ == model.a_matrix_.format_ && m_implementation->model->a_matrix_.start_ == model.a_matrix_.start_ && m_implementation->model->a_matrix_.index_ == model.a_matrix_.index_;

        if (settings.reuseBasis && m_implementation->model) {
            result.diagnostics.modelUpdateAttempted = true;
            result.diagnostics.modelUpdateApplied = sameStructure;
            result.diagnostics.modelStructureMismatch = !sameStructure;
        }

        if (!sameStructure) {
            const auto status = solver.passModel(model);

            if (status != HighsStatus::kOk) {
                return std::unexpected(std::format("HiGHS model pass was not exact: status={} columns={} rows={} nonzeros={}", static_cast<int>(status), model.num_col_, model.num_row_, model.a_matrix_.index_.size()));
            }
        } else {
            auto status = HighsStatus::kOk;

            if (m_implementation->model->col_cost_ != model.col_cost_) {
                status = solver.changeColsCost(0, model.num_col_ - 1, model.col_cost_.data());
            }

            if (status == HighsStatus::kOk && (m_implementation->model->col_lower_ != model.col_lower_ || m_implementation->model->col_upper_ != model.col_upper_)) {
                status = solver.changeColsBounds(0, model.num_col_ - 1, model.col_lower_.data(), model.col_upper_.data());
            }

            if (status == HighsStatus::kOk && (m_implementation->model->row_lower_ != model.row_lower_ || m_implementation->model->row_upper_ != model.row_upper_)) {
                status = solver.changeRowsBounds(0, model.num_row_ - 1, model.row_lower_.data(), model.row_upper_.data());
            }

            for (HighsInt row = 0; status == HighsStatus::kOk && row < model.num_row_; ++row) {
                for (auto entry = model.a_matrix_.start_[row]; entry < model.a_matrix_.start_[row + 1]; ++entry) {
                    if (m_implementation->model->a_matrix_.value_[entry] == model.a_matrix_.value_[entry]) {
                        continue;
                    }

                    status = solver.changeCoeff(row, model.a_matrix_.index_[entry], model.a_matrix_.value_[entry]);

                    if (status != HighsStatus::kOk) {
                        break;
                    }
                }
            }

            if (status != HighsStatus::kOk) {
                return std::unexpected("could not update a structure-stable HiGHS model");
            }
        }

        m_implementation->model = model;

        if (settings.reuseBasis && m_implementation->basis) {
            result.diagnostics.basisReuseAttempted = true;

            if (m_implementation->basis->valid && m_implementation->basis->col_status.size() == static_cast<std::size_t>(model.num_col_) && m_implementation->basis->row_status.size() == static_cast<std::size_t>(model.num_row_)) {
                if (solver.setBasis(*m_implementation->basis, "PathTempo SCP reuse") != HighsStatus::kOk) {
                    return std::unexpected(std::format("could not apply a dimension-checked HiGHS basis: columns={} rows={}", model.num_col_, model.num_row_));
                }

                result.diagnostics.basisReuseApplied = true;
            } else {
                result.diagnostics.basisDimensionMismatch = true;
                m_implementation->basis.reset();
            }
        }

        // HiGHS otherwise applies its time limit to cumulative runtime across repeated calls.
        solver.zeroAllClocks();
        const auto solveStatus = solver.run();
        const auto &solveInfo = solver.getInfo();

        if (solveInfo.simplex_iteration_count > 0) {
            result.diagnostics.simplexIterations = static_cast<std::size_t>(solveInfo.simplex_iteration_count);
        }

        const auto modelStatus = solver.getModelStatus();

        if (solveStatus != HighsStatus::kError && modelStatus == HighsModelStatus::kOptimal && solveStatus == HighsStatus::kOk) {
            result.status = LinearSolveStatus::Optimal;
        } else if (solveStatus != HighsStatus::kError && modelStatus == HighsModelStatus::kTimeLimit) {
            result.status = LinearSolveStatus::TimeLimit;
        } else if (solveStatus != HighsStatus::kError && modelStatus == HighsModelStatus::kIterationLimit) {
            result.status = LinearSolveStatus::IterationLimit;
        } else {
            return std::unexpected(std::format("HiGHS solve failed: {}", solver.modelStatusToString(modelStatus)));
        }

        if (result.status != LinearSolveStatus::Optimal) {
            return result;
        }

        const auto &solution = solver.getSolution();

        if (!solution.value_valid || solution.col_value.size() != program.columnCount()) {
            return std::unexpected("HiGHS solution has no primal values");
        }

        result.values = solution.col_value;

        if (settings.reuseBasis) {
            const auto &basis = solver.getBasis();

            if (basis.valid && basis.col_status.size() == static_cast<std::size_t>(model.num_col_) && basis.row_status.size() == static_cast<std::size_t>(model.num_row_)) {
                m_implementation->basis = basis;
            } else {
                m_implementation->basis.reset();
            }
        }

        return result;
    }
}
