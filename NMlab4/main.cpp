#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

// Функція
double equation(double x) {
    return x * x + 5 * sin(x) - 1;
}

// Генерування вузлів Чебишова
vector<double> chebyshevNodes(int n, double a, double b) {
    vector<double> nodes;
    for (int i = 0; i < n; ++i) {
        double node = ((a+b)/2)+(((b-a)/2)*cos(((2.0 * i + 1.0) * M_PI) / (2.0 * n)));
        nodes.push_back(node);
    }
    return nodes;
}

// Функція для обчислення полінома Ньютона
double newtonInterpolation(const vector<double>& x, const vector<double>& y, double xi) {
    double result = 0.0;
    vector<double> coefficients(y.size(), 0.0);

    for(int k = 0; k < x.size(); k++){
        for(int j = 0; j <= k; j++){
            double temp = 1.0;
            for(int i = 0; i <= k; i++){
                if(i != j){
                    temp *= (x[j] - x[i]);
                }
            }
            coefficients[k] +=  y[j] / temp;
        }
    }

    for (int i = 0; i < coefficients.size(); ++i) {
        double term = 1;
        for (int j = 0; j < i; ++j) {
            term = term * (xi - x[j]);
        }
        result += term * coefficients[i];
    }

    return result;
}

// Функція для обчислення полінома Лагранжа
double lagrangeInterpolation(const vector<double>& x, const vector<double>& y, double xi) {
    double result = 0.0;
    for (int i = 0; i < x.size(); ++i) {
        double term = 1;
        for (int j = 0; j < x.size(); ++j) {
            if (j != i) {
                term = term * (xi - x[j]) / (x[i] - x[j]);
            }
        }
        result += term * y[i];
    }
    return result;
}

// Метод дихотомії для знаходження кореня полінома Лагранжа
double bisectionMethodLang(const vector<double>& x, const vector<double>& y, double a, double b, double epsilon) {
    double xi = 0.0;
    double fa = lagrangeInterpolation(x, y, a);
    double fb = lagrangeInterpolation(x, y, b);

    if (fa * fb > 0) {
        cout << "It is not possible to apply the bisection method on this interval." << endl;
        return xi;
    }

    while (fabs(b - a) >= epsilon) {
        xi = (a + b) / 2;
        double fxi = lagrangeInterpolation(x, y, xi);

        if (fabs(fxi) < epsilon) {
            break;
        } else {
            if (fa * fxi < 0) {
                b = xi;
            } else {
                a = xi;
            }
        }
    }

    return xi;
}

// Метод дихотомії для знаходження кореня полінома Ньютона
double bisectionMethodNew(const vector<double>& x, const vector<double>& y, double a, double b, double epsilon) {
    double xi = 0.0;
    double fa = newtonInterpolation(x, y, a);
    double fb = newtonInterpolation(x, y, b);

    if (fa * fb > 0) {
        cout << "It is not possible to apply the bisection method on this interval." << endl;
        return xi;
    }

    while (fabs(b - a) >= epsilon) {
        xi = (a + b) / 2;
        double fxi = newtonInterpolation(x, y, xi);

        if (fabs(fxi) < epsilon) {
            break;
        } else {
            if (fa * fxi < 0) {
                b = xi;
            } else {
                a = xi;
            }
        }
    }

    return xi;
}

int main() {
    // Визначення проміжку та кількості вузлів
    double a = -3.0;
    double b = -1.0;
    int numNodes = 10;
    double epsilon = 1e-6;

    vector<double> x_values = chebyshevNodes(numNodes, a, b);
    vector<double> y_values;

    // Обчислення значень функції для відповідних x
    for (double x : x_values) {
        y_values.push_back(equation(x));
    }

    cout << "xk: {";
    for(int i = 0; i < x_values.size() - 1; i++) {
        cout << x_values[i] << ", ";
    }
    cout << x_values[x_values.size()-1] << "}\n";

    cout << "yk: {";
    for(int i = 0; i < y_values.size() - 1; i++) {
        cout << y_values[i] << ", ";
    }
    cout << y_values[y_values.size()-1] << "}\n";

    // Задання точки, в якій шукаємо значення функції
    double xi = bisectionMethodLang(x_values, y_values, a, b, epsilon);

    double func = equation(xi);

    cout << "Lagrange x: " << xi << endl;

    xi = bisectionMethodNew(x_values, y_values, a, b, epsilon);

    cout << "Newton x: " << xi << endl;

    cout << "f(x): " << func << endl;

    return 0;
}
