#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Функція для знаходження максимального елемента в стовпці
int findMaxRow(vector<vector<double> >& mat, int col, int startRow) {
    int maxRow = startRow;
    double maxVal = fabs(mat[startRow][col]);

    for (int i = startRow + 1; i < mat.size(); i++) {
        double curVal = fabs(mat[i][col]);
        if (curVal > maxVal) {
            maxVal = curVal;
            maxRow = i;
        }
    }

    return maxRow;
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

// Функція для знахлдження числа обумовленості
double conditionNumber(vector<vector<double>> data, vector<vector<double>>& inverse) {
    // Підрахунок норми Фробеніуса матриці
    double norm = 0.0;
    for (const auto& row : data) {
        for (double element : row) {
            norm += element * element;
        }
    }
    norm = sqrt(norm);

    // Підрахунок норми Фробеніуса оберненої матриці
    double inverseNorm = 0.0;
    for (const auto& row : inverse) {
        for (double element : row) {
            inverseNorm += element * element;
        }
    }
    inverseNorm = sqrt(inverseNorm);

    cout << "Frobenius norm matrix A: " << norm << endl;
    cout << "Frobenius norm inverse matrix A: " << inverseNorm << endl << endl;

    // Число обумовленості
    double conditionNumber = norm * inverseNorm;

    return conditionNumber;
}

vector<vector<double>> mult(vector<vector<double>>& A, vector<vector<double>>& M, int n){
    vector<vector<double>> result(n, vector<double>(n + 1));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += M[i][k] * A[k][j];
            }
        }
    }

    return result;
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
bool checkDet(vector<vector<double>> matrix){
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
bool checkT(vector<vector<double>> matrix){
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

void gaussMethod(vector<vector<double>>& A, int& swaps, vector<vector<double>>& inverse, vector<double>& X, double n, double& det) {
    vector<vector<double>> augmentedMatrix(n, vector<double>(2 * n));
    vector<vector<double>> M(n, vector<double>(n));
    vector<vector<double>> P(n, vector<double>(n));
    det = 1.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmentedMatrix[i][j] = A[i][j];
        }
        augmentedMatrix[i][i + n] = 1.0;
    }

    for (int i = 0; i < n; i++) {

        for (int k = 0; k < n; k++) {
            P[k][k] = 1.0;
            M[k][k] = 1.0;
            for(int j = 0; j < n; j++){
                if(j != k){
                    P[k][j] = 0.0;
                    M[k][j] = 0.0;
                }
            }
        }

        int maxRow = findMaxRow(A, i, i);

        if (maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(augmentedMatrix[i], augmentedMatrix[maxRow]);
            swap(P[i], P[maxRow]);
            swaps++;
            det *= -(1);
        }

        cout << "P" << i + 1 << ":" << endl;
        printMatrix(P);
        cout << endl;
        cout << "P" << i + 1 << "A" << i << ":" << endl;
        printMatrix(A);
        cout << endl;

        det *= A[i][i];

        double pivotElement = augmentedMatrix[i][i];
        for (int j = i; j < 2 * n; j++) {
            augmentedMatrix[i][j] /= pivotElement;
        }

        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = augmentedMatrix[j][i];
                for (int k = i; k <= 2 * n; k++) {
                    augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
                }
            }
        }

        M[i][i] = 1/A[i][i];
        for(int j = 0; j < n; j++){
            if(j != i) {
                M[j][i] = -(A[j][i] / A[i][i]);
            }
        }

        cout << "M" << i + 1 << ":" << endl;
        printMatrix(M);
        cout << endl;

        A = mult(A, M, n);

    }

    cout << "A" << n << ":" << endl;
    printMatrix(A);
    cout << endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse[i][j] = augmentedMatrix[i][j + n];
        }
    }

    cout << "Inverse:\n";
    printMatrix(inverse);
    cout << endl;

    if (fabs(det) == 0) {
        cout << "The system has many solutions or no solutions" << endl;
        return;
    }

    // Обертаємо матрицю та знаходимо розв'язок
    for (int i = n - 1; i >= 0; i--) {
        X[i] = A[i][n] / A[i][i];
        for (int j = i - 1; j >= 0; j--) {
            A[j][n] -= A[j][i] * X[i];
        }
    }
}

vector<double> seidelMethod(vector<vector<double>> A, vector<double> b, int maxIterations, double epsilon, int n) {
    vector<double> previousSolution(n, 0.0);
    vector<double> currentSolution(n, 0.0);

    for (int iteration = 0; iteration < maxIterations; iteration++) {
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;

            for (int j = 0; j < i; j++) {
                if (j != i) {
                    sum1 += A[i][j] * currentSolution[j];
                }
            }
            for (int j = i + 1; j < n; j++) {
                if (j != i) {
                    sum2 += A[i][j] * currentSolution[j];
                }
            }
            currentSolution[i] = (b[i] - sum1 - sum2)/A[i][i];

        }

        double error = 0.0;

        for (int i = 0; i < n; i++) {
            error = max(error, abs(previousSolution[i] - currentSolution[i]));
        }

        for(int i = 0; i < n; i++){
            cout << "x" << i + 1 << ":" << currentSolution[i] << " ";
        }
        cout << "error:" << error << endl;

        if (error < epsilon) {
            return currentSolution;
        }

        previousSolution = currentSolution;
    }
    return currentSolution;
}

int main() {
    // Вхідні дані - матриця та вектор вільних членів
    //const vector<vector<double> > A = {
    //        {10, 1, -1, 1, 3},
    //        {1, 2, 0, 2, 2},
    //        {-1, 0, 8, 3, 4},
    //        {1, 2, 3, 10, 8},
    //        {3, 2, 4, 8, 10}
    //};
//
    //const vector<double> B = {4, 10, 4, 12, 8};

    const vector<vector<double> > A = {
            {2, 1, -1, 1, 3},
            {1, 2, 3, 2, 7},
            {1, 0, 1, 3, 4},
            {3, 2, -3, 4, 8},
            {5, 2, 2, 1, 10}
    };

    const vector<double> B = {4, 10, 4, 12, 8};

    vector<vector<double>> matrix = A;

    int n = matrix.size();

    // Створюємо початкову розширену матрицю [A | B]
    for (int i = 0; i < n; i++) {
        matrix[i].push_back(B[i]);
    }

    // Розв'язуємо систему рівнянь методом Гауса з вибором головного елемента
    cout << "Gauss:" << endl;
    vector<vector<double>> inverse(n, vector<double>(n));
    vector<double> X(n);
    double det;
    int swaps = 0;
    gaussMethod(matrix, swaps, inverse, X, n, det);

    // Обчислюємо число обумовленості
    double cond = conditionNumber(A, inverse);

    // Виводимо розв'язок
    cout << "Solution:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x[" << i + 1 << "] = " << X[i] << endl;
    }

    // Розв'язуємо систему рівнянь методом Зейделя
    cout << "\nSeidel:" << endl;
    if(checkDet(A) && checkT(A)) {
        vector<double> x = seidelMethod(A, B, 100, 0.0001, n);

        cout << "\nSolution:" << endl;
        for (int i = 0; i < n; i++) {
            cout << "x[" << i + 1 << "] = " << x[i] << endl;
        }
    }
    else{
        cout << "Condition of convergence is not fulfilled" << endl;
    }

    // Виводимо визначник та число обумовленості
    cout << "\nDeterminant: " << det << endl;
    cout << "Condition Number Frobenius norm: " << cond << endl << endl;

    return 0;
}
