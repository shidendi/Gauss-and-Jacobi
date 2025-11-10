#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
using namespace std;

// вывод
void printMatrix(const vector<vector<double>>& a) {
    cout << "\nТекущая матрица:\n";
    cout << fixed << setprecision(6);
    for (const auto& row : a) {
        for (int j = 0; j < (int)row.size() - 1; ++j)
            cout << setw(12) << row[j] << " ";
        cout << " | " << setw(12) << row.back() << "\n";
    }
    cout << "\n";
}

// диагональное преобладание для Якобииии
bool isDiagonallyDominant(vector<vector<double>>& a, double EPS) {
    int n = a.size(); // количество уравнений (строк)
    int m = a[0].size() - 1; // количество неизвестных (без последнего столбца)
    for (int i = 0; i < n; ++i) {
        double diag = fabs(a[i][i]); // модуль диагонального элемента
        double sum = 0; // сумма модулей внедиагональных элементов
        for (int j = 0; j < m; ++j)
            if (i != j) 
                sum += fabs(a[i][j]);
        // если диагональный элемент меньше суммы остальных тогда преобладания нет
        if (diag < sum - EPS)
            return false;
    }
    return true;
}

// try перестановок для преобладания
void makeDiagonalDominant(vector<vector<double>>& a) {
    int n = a.size();
    int m = a[0].size() - 1;
    for (int i = 0; i < n; ++i) {
        int bestRow = i; // предполагаем, что текущая строка — лучшая
        double maxDiag = fabs(a[i][i]);  // текущий модуль диагонального элемента
        // поиск строки ниже, у которой в текущем столбце элемент по модулю больше
        for (int k = i + 1; k < n; ++k) {
            if (fabs(a[k][i]) > maxDiag) {
                maxDiag = fabs(a[k][i]);
                bestRow = k; // строка с бОльшим элементом
            }
        }
        if (bestRow != i) 
            swap(a[i], a[bestRow]); // свап
    }
}

// Гаусссс
void solveGauss(vector<vector<double>> a, double EPS) {
    cout << "\n=== РЕШЕНИЕ МЕТОДОМ ГАУССА ===\n";
    cout << "\nИсходная матрица:"; printMatrix(a);
    int n = a.size();
    int m = a[0].size() - 1;
    int row = 0; // текущая обрабатываемая строка
    vector<int> where(m, -1); // вектор для отслеживания позиций ведущих элементов

    // прямой ход метода Гаусса
    for (int col = 0; col < m && row < n; ++col) {
        cout << "==== Обработка столбца " << (col + 1) << " ====\n";
        int sel = row; // строка с максимальным элементом
        // поиск строки, где текущий элемент в столбце максимален по модулю
        for (int i = row; i < n; ++i)
            if (fabs(a[i][col]) > fabs(a[sel][col])) 
                sel = i;
        // если элемент слишком мал скип столбца
        if (fabs(a[sel][col]) < EPS) {
            cout << "Столбец " << col + 1 << " почти нулевой, пропускаем.\n\n";
            continue;
        }
        // свап если нашли лучший ведущий элемент
        if (sel != row) {
            cout << "Меняем строки " << (row + 1) << " и " << (sel + 1)
                << " (частичный выбор главного элемента)\n";
            swap(a[sel], a[row]);
            printMatrix(a);
        }
        // для этого столбца ведущая строка = row
        where[col] = row;

        // нормализация текущей строки (деление на ведущий элемент)
        double div = a[row][col];
        cout << "Нормализуем строку " << (row + 1) << " (делим на " << div << ")\n";
        for (int j = col; j <= m; ++j) a[row][j] /= div;
        printMatrix(a);
        
        // обнуление остальных в столбце
        cout << "Обнуляем остальные элементы в столбце " << (col + 1) << "\n";
        for (int i = 0; i < n; ++i) {
            if (i == row) 
                continue;
            double factor = a[i][col];
            if (fabs(factor) < EPS) 
                continue;
            for (int j = col; j <= m; ++j) 
                a[i][j] -= factor * a[row][j];
        }
        printMatrix(a);
        ++row;
    }

    // чек на несовместность
    for (int i = 0; i < n; ++i) {
        bool allzero = true;
        for (int j = 0; j < m; ++j)
            if (fabs(a[i][j]) > EPS) { 
                allzero = false;
                break; 
            }
        // если все коэффициенты нули, но правая часть не ноль = система несовместна
        if (allzero && fabs(a[i][m]) > EPS) {
            cout << "Система несовместна (решений нет).\n";
            return;
        }
    }

    // поиск решеения
    vector<double> ans(m, 0.0);
    for (int i = 0; i < m; ++i) 
        if (where[i] != -1)
            ans[i] = a[where[i]][m]; // значение переменной из строки с ведущим элементом

    // кол-во свободных переменных
    int free_vars = 0;
    for (int i = 0; i < m; ++i) 
        if (where[i] == -1) 
            ++free_vars;

    cout << "\n========================\n";
    cout << "Итоговая матрица после преобразований:"; printMatrix(a);
    cout << setprecision(10);
    if (free_vars > 0) {
        cout << "Бесконечно много решений. Число свободных переменных: " << free_vars << "\n";
        cout << "Одно частичное решение (свободные переменные = 0):\n";
        for (int i = 0; i < m; ++i)
            cout << "x" << (i + 1) << " = " << ans[i] << (i + 1 == m ? "\n" : " ");
        cout << "\n";
    }
    else {
        cout << "Единственное решение:\n";
        for (int i = 0; i < m; ++i) 
            cout << "x" << (i + 1) << " = " << ans[i] << (i + 1 == m ? "\n" : " ");
    }
}

// Якобииии
void solveJacobi(vector<vector<double>> a, double EPS, int maxIter = 1000) {
    cout << "\n=== РЕШЕНИЕ МЕТОДОМ ЯКОБИ ===\n";

    // диаг преобладание
    if (!isDiagonallyDominant(a, EPS)) {
        cout << "Матрица не имеет диагонального преобладания. Пробуем переставить строки...\n";
        makeDiagonalDominant(a);
        if (!isDiagonallyDominant(a, EPS)) {
            cout << "После перестановки диагональное преобладание не достигнуто. Метод Якоби может не сойтись.\n";
        }
        else {
            cout << "После перестановки строки матрица имеет диагональное преобладание.\n";
        }
    }
    else {
        cout << "Матрица имеет диагональное преобладание.\n";
    }

    int n = a.size();
    int m = a[0].size() - 1;
    vector<double> x(m, 0.0); // начальное приближение (нулевое)
    vector<double> x_new(m, 0.0); // новое приближение на каждой итерации

    cout << fixed << setprecision(8);
    for (int iter = 1; iter <= maxIter; ++iter) {
        // новое приближение x_new
        for (int i = 0; i < n; ++i) {
            double sum = a[i][m]; // правая часть
            for (int j = 0; j < m; ++j) 
                if (j != i) 
                    sum -= a[i][j] * x[j]; // перенос всё кроме диагонального
            x_new[i] = sum / a[i][i]; // деление на диагональный элемент
        }

        // максимальное изменение между старым и новым вектором
        double diff = 0.0;
        for (int i = 0; i < m; ++i) 
            diff = max(diff, fabs(x_new[i] - x[i]));

        cout << "Итерация " << setw(3) << iter << ": ";
        for (int i = 0; i < m; ++i) 
            cout << setw(12) << x_new[i] << " ";
        cout << " | макс. разница = " << scientific << diff << fixed << "\n";

        // проверка критерия сходимости
        if (diff < EPS) {
            cout << "\nМетод Якоби сошёлся за " << iter << " итераций.\n";
            cout << "Решение:\n";
            for (int i = 0; i < m; ++i) 
                cout << setw(12) << x_new[i] << " ";
            cout << "\n";
            return;
        }
        x = x_new;
    }
    cout << "Метод Якоби не сошёлся за " << maxIter << " итераций.\n";
}

int main() {
    setlocale(LC_ALL, "Russian");
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    cout << "Введите число уравнений n и число неизвестных m (через пробел): \n";
    if (!(cin >> n >> m)) 
        return 0;
    if (n <= 0 || m <= 0) {
        cerr << "n и m должны быть положительными целыми числами\n";
        return 0;
    }

    vector<vector<double>> a(n, vector<double>(m + 1));
    cout << "Введите расширенную матрицу (n строк, m+1 столбец — правые части):\n";
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= m; ++j)
            cin >> a[i][j];

    double EPS;
    cout << "Введите EPS (малое число для сравнения с нулём), или 0 чтобы взять значение по умолчанию 1e-9:\n";
    cin >> EPS;
    if (EPS <= 0) 
        EPS = 1e-9;

    cout << "\nВыберите метод решения:\n";
    cout << "1 — Метод Гаусса\n";
    cout << "2 — Метод Якоби\n";
    int choice;
    cin >> choice;
    if (choice == 1) 
        solveGauss(a, EPS);
    else if 
        (choice == 2) 
        solveJacobi(a, EPS);
    else cout << "Неверный выбор метода!\n";

    return 0;
}
