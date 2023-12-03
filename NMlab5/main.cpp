#include <iostream>
#include <cmath>

using namespace std;

// Функція, яку інтегруємо
double f(double x) {
    return 1.0 / (2 + x);
}

// Метод Сімпсона для обчислення інтегралу
double simpson(double a, double b, int n) {
    double h = (b - a) / n;
    double integral = f(a) + f(b);

    for (int i = 1; i < n; i += 2) {
        integral += 4 * f(a + i * h);
    }

    for (int i = 2; i < n - 1; i += 2) {
        integral += 2 * f(a + i * h);
    }

    return (h / 3) * integral;
}

// Правило Рунге для оцінки точності
double runge(double I1, double I2, int p) {
    return abs(I2 - I1) / (pow(2, p) - 1);
}

int main() {
    double a = 1.0;  // Нижня межа інтегрування
    double b = 5.0;  // Верхня межа інтегрування
    double tolerance = 0.005;  // Точність
    int n = 2;  // Початкова кількість підінтервалів
    double I1, I2, R;

    // Обчислення інтегралу та оцінка точності за допомогою правила Рунге
    do {
        cout << "Iteration:\n";
        I1 = simpson(a, b, n);
        I2 = simpson(a, b, 2 * n);
        R = runge(I1, I2, 4);
        cout << "h1: " << (b - a) / n << endl;
        cout << "I1: " << I1 << endl;
        cout << "h2: " << (b - a) / (n*2) << endl;
        cout << "I2: " << I2 << endl;
        cout << "R: " << R << endl;

        n *= 2;  // Подвоєння кількості підінтервалів для найбільшої точності

    } while (R > tolerance);

    // Виведення результатів
    cout << "\nResult: " << I2 << endl;
    cout << "n: " << n << endl;
    cout << "h: " << (b - a) / n << endl;
    cout << "R: " << R << endl;

    return 0;
}