#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// скалярне множення
double scalarMult(vector<double>& a, vector<double>& b){
    double result = 0.0;
    for(int i = 0; i < a.size(); i++){
        result += a[i] * b[i];
    }
    return result;
}

// Пошук максимального значення методом скалярних добутків
double computeMaxEigenvalue(int size, vector<vector<double>>& data) {
    vector<double> x(size, 1.0);
    vector<double> e(size, 0.0);
    double maxEigenvalue = 0.0;
    double prevEigenvalue = 0.0;
    double epsilon = 0.2;
    int iterations = 0;
    double norm;

    do {
        prevEigenvalue = maxEigenvalue;
        maxEigenvalue = 0.0;

        norm = 0.0;
        for(int i = 0; i < size; i++){
            norm = max(norm, abs(x[i]));
        }

        for(int i = 0; i < size; i++){
            e[i] = x[i] / norm;
        }

        cout << "x" << iterations << " = (";
        for(int i = 0; i < size; i++) {
            if(i != size - 1) {
                cout << x[i] << ", ";
            }
            else{
                cout << x[i];
            }
        }
        cout << ")\ne" << iterations << " = (";
        for(int i = 0; i < size; i++) {
            if(i != size - 1) {
                cout << e[i] << ", ";
            }
            else{
                cout << e[i];
            }
        }
        cout << ")\n";

        for (int i = 0; i < size; i++) {
            x[i] = 0.0;

            for (int j = 0; j < size; j++) {
                x[i] += data[i][j] * e[j];
            }
        }

        maxEigenvalue = scalarMult(x, e)/ scalarMult(e, e);

        cout << "Eigenvalue" << iterations + 1 << " = " << maxEigenvalue << endl;

        iterations++;
    } while (abs(maxEigenvalue - prevEigenvalue) > epsilon && iterations < 100);
    return maxEigenvalue;
}

double calculateDeterminant(vector<vector<double>>& matrix, int n) {
    if(n == 1){
        return matrix[0][0];
    }

    double det = 0;
    vector<vector<double>> submatrix(n, vector<double>(n));
    if (n == 2)
        return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
    else {
        for (int x = 0; x < n; x++) {
            int subi = 0;
            for (int i = 1; i < n; i++) {
                int subj = 0;
                for (int j = 0; j < n; j++) {
                    if (j == x)
                        continue;
                    submatrix[subi][subj] = matrix[i][j];
                    subj++;
                }
                subi++;
            }
            det = det + (pow(-1, x) * matrix[0][x] * calculateDeterminant( submatrix, n - 1 ));
        }
    }
    return det;
}

// Переврка збіжності визначник
bool checkDet(vector<vector<double>>& matrix){
    int n = matrix.size();

    for(int i = 1; i <= n; i++){
        if(calculateDeterminant(matrix, i) <= 0){
            cout << "det: " << calculateDeterminant(matrix, i) << endl;
            return false;
        }
    }

    return true;
}

// Переврка збіжності транспанована матриця
bool checkT(vector<vector<double>>& matrix){
    int n = matrix.size();
    vector<vector<double>> matrixT(n, vector<double>(n));

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            matrixT[i][j] = matrix[j][i];
        }
    }
    if(matrixT != matrix){
        cout << "A != A^T" << endl;
        return false;
    }
    return true;
}

// Функція для виводу матриці
void printMatrix(vector<vector<double>>& mat) {
    for (int i = 0; i < mat.size(); i++) {
        for (int j = 0; j < mat[i].size(); j++) {
            cout << mat[i][j] << "\t";
        }
        cout << endl;
    }
}

// Пошук мінімального значення
double computeMinEigenvalue(int size, vector<vector<double>>& data) {
    double minEigenvalue = 0.0;
    if(checkT(data) && checkDet(data)) {
        vector<double> sum(size, 0.0);
        double norm = 0.0;
        for(int i = 0; i < size; i++){
            for(int j = 0; j < size; j++){
                sum[i] += abs(data[i][j]);
            }
            norm = max(norm, sum[i]);
        }

        vector<vector<double>> B(size, vector<double>(size, 0.0));

        for(int i = 0; i < size; i++){
            B[i][i] = norm;
            for(int j = 0; j < size; j++){
                B[i][j] -= data[i][j];
            }
        }
        cout << "B = {\n";
        printMatrix(B);
        cout << "}\n";

        double maxEigenvalue = computeMaxEigenvalue(size, B);
        minEigenvalue = norm - maxEigenvalue;

        cout << "||A|| = " << norm << endl;
    }

    return minEigenvalue;
}

int main() {

    vector<vector<double>> matrixData = {
            {2, 1, 0},
            {1, 2, 1},
            {0, 1, 2}
    };

    cout << "Maximum eigenvalue:" << endl;
    double maxEigenvalue = computeMaxEigenvalue(matrixData.size(), matrixData);
    cout << "\nMinimum eigenvalue:" << endl;
    double minEigenvalue = computeMinEigenvalue(matrixData.size(), matrixData);

    cout << "\nMaximum eigenvalue: " << maxEigenvalue << endl;
    cout << "Minimum eigenvalue: " << minEigenvalue << endl;

    return 0;
}