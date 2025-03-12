#include <iostream>
#include <cmath>
#include <vector>
#include <limits> 

int n = 2; //размерность вектора x
double eps1 = 0.01; //точность решения задачи
double eps2 = 0.01;//точность решения задачи безусловной минимизации
double r = 100.0; // начальное значение параметра штрафа
double C = 10.0; // число для уменьшения параметра штрафа
double alpha = 0.01; //шаг
int max_iter = 1000; //максимальное кол-во итераций

// Целевая функция
double f(const std::vector<double>& x) {
    double x1 = x[0];
    double x2 = x[1];
    return pow(x1 - 2, 2) + pow(x2 - 3, 2);
}

double g(const std::vector<double>& x, int i) {  //ограничения-неравенства
    if (i == 0) return x[0] + x[1] - 9;
    if (i == 1) return x[0] + 2 * x[1] - 12;
    return 0;
}

double phi(const std::vector<double>& x, double r) {
    double penalty = 0.0;
    for (int i = 0; i < 2; ++i) {
        if (g(x, i) >= 0) {
            return std::numeric_limits<double>::infinity();
        }
        penalty += log(-g(x, i));  //штрафная функция
    }
    return f(x) - r * penalty;
}

std::vector<double> gradient_phi(const std::vector<double>& x, double r, double h = 1e-5) {   //градиент штрафной функции через конечные разности
    int n = x.size();
    std::vector<double> grad(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        std::vector<double> x_plus_h = x;
        x_plus_h[i] += h;
        std::vector<double> x_minus_h = x;
        x_minus_h[i] -= h;

        double phi_plus = phi(x_plus_h, r);
        double phi_minus = phi(x_minus_h, r);

        if (std::isnan(phi_plus) || std::isinf(phi_plus) || std::isnan(phi_minus) || std::isinf(phi_minus)) {
            return std::vector<double>(n, 0.0);
        }

        grad[i] = (phi_plus - phi_minus) / (2 * h);
    }
    return grad;
}


// Метод градиентного спуска
std::vector<double> gradient_descent(const std::vector<double>& x0, double r, double alpha, int max_iter) {
    std::vector<double> x = x0;
    for (int iter = 0; iter < max_iter; ++iter) {
        auto grad = gradient_phi(x, r);

        bool nan_grad = false;
        for (double val : grad) {
            if (std::isnan(val) || std::isinf(val)) {
                std::cout << "Gradient became NaN or infinite. Stopping descent.\n";
                nan_grad = true;
                break;
            }
        }
        if (nan_grad) {
            break;
        }

        for (size_t i = 0; i < x.size(); ++i) {
            x[i] -= alpha * grad[i];
        }
        bool constraint_violated = false;
        for (int i = 0; i < 2; ++i) {
            if (g(x, i) >= 0) {
                std::cout << "Constraint " << i << " violated.  Stopping descent.\n";
                constraint_violated = true;
                break;
            }
        }
        if (constraint_violated) break;
    }
    return x;
}

int main() {
    setlocale(LC_ALL, "RUS");

    std::vector<double> x0 = { 0.0, 0.0 }; // начальное приближение
    std::vector<double> x = x0;

    for (int k = 0; k < 10; ++k) {
        x = gradient_descent(x, r, alpha, max_iter);
        r /= C;
        std::cout << "Номер итерациия " << k << ": x = (" << x[0] << ", " << x[1] << ")" << std::endl;
    }

    std::cout << "Точка экстремума x = (" << x[0] << ", " << x[1] << ")\n";
    return 0;
}
