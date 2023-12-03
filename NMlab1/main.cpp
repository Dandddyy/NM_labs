#include <iostream>
#include <cmath>

using namespace std;

// f(x)
double f(double x) {
    return pow(x, 3) + 4 * x - 6;
}

// g(x)
double g(double x) {
    return (6 - pow(x, 3)) / 4;
}

// f`(x)
double df(double x) {
    return 3 * pow(x, 2) + 4;
}

// f``(x)
double ddf(double x) {
    return 6 * x;
}

// g`(x)
double dg(double x) {
    return 0-(3 * pow(x, 2))/4;
}

//double SIpriori(double x0, double epsilon, double q) {
//    return ceil(log(fabs(g(x0) - x0)/((1 - q) * epsilon))/log(1/q)) + 1;
//}

// Апріорна оцінка для методу Релаксацій
double REpriori(double z0, double epsilon, double q) {
    return ceil(log(fabs(z0)/epsilon)/log(1/q)) + 1;
}

// Апріорна оцінка для методу Ньютона
double NEpriori(double z0, double epsilon, double q) {
    return ceil(log2(log(fabs(z0)/epsilon)/log(1/q)) + 1) + 1;
}

//void simpleIterations(double x0, double epsilon) {
//
//    double q;
//
//    double x = x0;
//
//    int iter = 0;
//
//    cout << "Simple iterations Method:" << endl;
//
//    while (true) {
//        double x_new = g(x);
//        iter++;
//        if (fabs(x_new - x) < epsilon) {
//            cout << "Iteration: " << iter << "   x: " << x_new << "   |x(k) - x(k - 1)|: " << fabs(x_new - x) << endl;
//            q = fabs(dg(1.15));
//            x = x_new;
//            break;
//        }
//        cout << "Iteration: " << iter << "   x: " << x_new << "   |x(k) - x(k - 1)|: " << fabs(x_new - x) << endl;
//        x = x_new;
//    }
//
//    cout << "Root: " << x << endl;
//    cout << "Posterior Estimate: " << iter << endl << endl;
//    cout << "Priori Estimate: " << SIpriori(x0, epsilon, q) << endl << endl;
//}

void relaxation(double x0, double epsilon, double relaxParam) {

    double q;

    double x = x0;

    int iter = 0;

    cout << "Relaxation Method:" << endl;

    while (true) {
        double x_new = x - relaxParam * f(x);
        iter++;
        if (fabs(x_new - x) < epsilon) {
            cout << "Iteration: " << iter << "   x: " << x_new << "   |x(k) - x(k - 1)|: " << fabs(x_new - x) << endl;
            //q = fabs(1 - relaxParam * df(1.15));
            q = (fabs(df(1.15)) - fabs(df(1.1)))/(fabs(df(1.15)) + fabs(df(1.1)));
            cout << q << endl;
            x = x_new;
            break;
        }
        cout << "Iteration: " << iter << "   x: " << x_new << "   |x(k) - x(k - 1)|: " << fabs(x_new - x) << endl;
        x = x_new;
    }

    cout << "Root: " << x << endl;
    cout << "Posterior Estimate: " << iter << endl;
    cout << "Priori Estimate: " << REpriori(x0 - 1.15, epsilon, q) << endl << endl;
}

void newton(double x0, double epsilon) {

    double q;

    double x = x0;

    int iter = 0;

    cout << "Newton Method:" << endl;

    while (true) {
        double x_new = x - f(x) / df(x);
        iter++;

        if (fabs(x_new - x) < epsilon) {
            cout << "Iteration: " << iter << "   x: " << x_new << "   |x(k) - x(k - 1)|: " << fabs(x_new - x) << endl;
            q = (fabs(ddf(1.15))*fabs(x0 - 1.15))/(2*fabs(ddf(1.1)));
            x = x_new;
            break;
        }
        cout << "Iteration: " << iter << "   x: " << x_new << "   |x(k) - x(k - 1)|: " << fabs(x_new - x) << endl;
        x = x_new;
    }

    cout << "Root: " << x << endl;
    cout << "Posterior Estimate: " << iter << endl;
    cout << "Priori Estimate: " << NEpriori(x0 - 1.15, epsilon, q) << endl << endl;
}

int main() {

    const double epsilon = 1e-4;
    double x0 = 1.14;  // Початкове наближення

    // Метод релаксації
    double relaxParam = 0.13;  // Параметр релаксації
    relaxation(x0, epsilon, relaxParam);

    // Метод простих ітерацій
    //simpleIterations(x0, epsilon);

    // Метод Ньютона
    newton(x0, epsilon);

    return 0;
}
