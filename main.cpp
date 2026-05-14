#include <iostream>
#include <vector>
#include <chrono>

#include <random>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

using namespace chrono;

// типы для удобства
using Matr = vector<vector<double>>;
using Vec = vector<double>;

void print_mat(const Matr& a) {
    for (size_t i = 0; i < a.size(); i++) {
        for (size_t j = 0; j < a[i].size(); j++) {
            cout << setw(10) << setprecision(4) << a[i][j] << " ";
        }
        cout << endl;
    }
}

void print_vec(const Vec& v) {
    for (double val : v) {
        cout << setw(12) << setprecision(6) << val << " ";
    }
    cout << endl;
}

// генерация случ матрицы
Matr rand_mat(int n, unsigned seed = 42) {
    mt19937 rng(seed);
    uniform_real_distribution<double> dist(-1.0, 1.0);
    
    Matr a(n, Vec(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i][j] = dist(rng);
        }
    }
    return a;
}

// генерация правой части
Vec rand_rhs(int n, unsigned seed = 123) {
    mt19937 rng(seed);
    uniform_real_distribution<double> dist(-1.0, 1.0);
    
    Vec b(n);
    for (int i = 0; i < n; ++i) {
        b[i] = dist(rng);
    }
    return b;
}

// матрица Гильберта
Matr hilbert_mat(int n) {
    Matr h(n, Vec(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            h[i][j] = 1.0 / (i + j + 1);
        }
    }
    return h;
}

// считаем невязку
double norm_resid(const Matr& a, const Vec& x, const Vec& b) {
    int n = a.size();
    double max_err = 0.0;
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) {
            s += a[i][j] * x[j];
        }
        double err = fabs(s - b[i]);
        if (err > max_err) max_err = err;
    }
    return max_err;
}

// копируем матрицу
Matr copy_mat(const Matr& src) {
    int n = src.size();
    Matr dst(n, Vec(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dst[i][j] = src[i][j];
        }
    }
    return dst;
}

//  метод гаусса без выбора 

bool gauss1(Matr a, Vec b, Vec& x) {
    int n = a.size();
    x.assign(n, 0.0);
    
    // прямой ход
    for (int k = 0; k < n; ++k) {
        if (fabs(a[k][k]) < 1e-12) {
            return false;
        }
        
        for (int i = k + 1; i < n; ++i) {
            double m = a[i][k] / a[k][k];
            for (int j = k; j < n; ++j) {
                a[i][j] -= m * a[k][j];
            }
            b[i] -= m * b[k];
        }
    }
    
    // обратный ход
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
    return true;
}

// метод гаусса с выбором 

bool gauss2(Matr a, Vec b, Vec& x) {
    int n = a.size();
    x.assign(n, 0.0);
    
    for (int k = 0; k < n; ++k) {
        // ищем главный элемен
        int best = k;
        double best_val = fabs(a[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (fabs(a[i][k]) > best_val) {
                best_val = fabs(a[i][k]);
                best = i;
            }
        }
        
        if (best != k) {
            swap(a[k], a[best]);
            swap(b[k], b[best]);
        }
        
        if (fabs(a[k][k]) < 1e-12) {
            return false;
        }
        
        for (int i = k + 1; i < n; ++i) {
            double m = a[i][k] / a[k][k];
            for (int j = k; j < n; ++j) {
                a[i][j] -= m * a[k][j];
            }
            b[i] -= m * b[k];
        }
    }
    
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
    return true;
}

// lu разложение 

bool lu_decomp(const Matr& a, Matr& l, Matr& u) {
    int n = a.size();
    l.assign(n, Vec(n, 0.0));
    u.assign(n, Vec(n, 0.0));
    
    for (int i = 0; i < n; ++i) {
        l[i][i] = 1.0;
    }
    
    for (int i = 0; i < n; ++i) {
        // считаем U[i][j]
        for (int j = i; j < n; ++j) {
            double s = 0.0;
            for (int k = 0; k < i; ++k) {
                s += l[i][k] * u[k][j];
            }
            u[i][j] = a[i][j] - s;
        }
        
        // считаем L[j][i]
        for (int j = i + 1; j < n; ++j) {
            double s = 0.0;
            for (int k = 0; k < i; ++k) {
                s += l[j][k] * u[k][i];
            }
            if (fabs(u[i][i]) < 1e-12) {
                return false;
            }
            l[j][i] = (a[j][i] - s) / u[i][i];
        }
    }
    return true;
}

Vec forward_solve(const Matr& l, const Vec& b) {
    int n = l.size();
    Vec y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < i; ++j) {
            s += l[i][j] * y[j];
        }
        y[i] = b[i] - s;
    }
    return y;
}

Vec back_solve(const Matr& u, const Vec& y) {
    int n = u.size();
    Vec x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        if (fabs(u[i][i]) < 1e-12) {
            return Vec();
        }
        double s = 0.0;
        for (int j = i + 1; j < n; ++j) {
            s += u[i][j] * x[j];
        }
        x[i] = (y[i] - s) / u[i][i];
    }
    return x;
}

bool solve_lu(const Matr& a, const Vec& b, Vec& x) {
    Matr L, U;
    if (!lu_decomp(a, L, U)) {
        return false;
    }
    Vec y = forward_solve(L, b);
    x = back_solve(U, y);
    return true;
}

bool solve_lu_multi(const Matr& a, const vector<Vec>& all_b, vector<Vec>& all_x) {
    Matr L, U;
    if (!lu_decomp(a, L, U)) {
        return false;
    }
    
    all_x.clear();
    for (const auto& bb : all_b) {
        Vec y = forward_solve(L, bb);
        Vec xx = back_solve(U, y);
        all_x.push_back(xx);
    }
    return true;
}

// замер времени 

template<typename F>
double time_me(F&& f) {
    auto start = high_resolution_clock::now();
    f();
    auto end = high_resolution_clock::now();
    return duration_cast<microseconds>(end - start).count() / 1000.0;
}

// эксперимент 1 

void exp1() {
    cout << "\n" << string(80, '=') << endl;
    cout << "эксперимент 1: время одной системы" << endl;
    cout << string(80, '=') << endl;
    
    vector<int> sizes = {100, 200, 500, 1000};
    
    cout << left << setw(10) << "n"
         << setw(20) << "Гаусс(без)"
         << setw(20) << "Гаусс(с)"
         << setw(20) << "LU разл"
         << setw(20) << "LU подст"
         << setw(20) << "LU всего" << endl;
    cout << string(110, '-') << endl;
    
    for (int n : sizes) {
        cout << "генерация " << n << "..." << endl;
        
        Matr A = rand_mat(n, 42);
        Vec b = rand_rhs(n, 123);
        
        double t1 = 0, t2 = 0, t_lu_fact = 0, t_lu_solve = 0, t_lu_tot = 0;
        
        // гаусс без выбора
        {
            Matr Ac = copy_mat(A);
            Vec bc = b;
            Vec x;
            t1 = time_me([&]() { gauss1(Ac, bc, x); });
        }
        
        // гаусс с выбором
        {
            Matr Ac = copy_mat(A);
            Vec bc = b;
            Vec x;
            t2 = time_me([&]() { gauss2(Ac, bc, x); });
        }
        
        // LU
        {
            Matr L, U;
            t_lu_fact = time_me([&]() { lu_decomp(A, L, U); });
            
            t_lu_solve = time_me([&]() {
                Matr L2, U2;
                lu_decomp(A, L2, U2);
                Vec y = forward_solve(L2, b);
                Vec x = back_solve(U2, y);
            });
            
            t_lu_tot = t_lu_fact + t_lu_solve;
        }
        
        cout << left << setw(10) << n
             << setw(20) << fixed << setprecision(4) << t1
             << setw(20) << t2
             << setw(20) << t_lu_fact
             << setw(20) << t_lu_solve
             << setw(20) << t_lu_tot << endl;
    }
}

// эксперимент 2 

void exp2() {
    cout << "\n" << string(80, '=') << endl;
    cout << "ЭКСПЕРИМЕНТ 2: много правых частей (n=500)" << endl;
    cout << string(80, '=') << endl;
    
    const int n = 500;
    vector<int> ks = {1, 10, 100};
    
    cout << "генерация матрицы..." << endl;
    Matr A = rand_mat(n, 42);
    
    cout << left << setw(15) << "k"
         << setw(25) << "Гаусс (всего мс)"
         << setw(25) << "LU (всего мс)"
         << setw(20) << "ускорение" << endl;
    cout << string(85, '-') << endl;
    
    for (int k : ks) {
        vector<Vec> all_b;
        for (int i = 0; i < k; ++i) {
            all_b.push_back(rand_rhs(n, 123 + i));
        }
        
        double tg = time_me([&]() {
            for (int i = 0; i < k; ++i) {
                Matr Ac = copy_mat(A);
                Vec bc = all_b[i];
                Vec x;
                gauss2(Ac, bc, x);
            }
        });
        
        double tlu = time_me([&]() {
            vector<Vec> xx;
            solve_lu_multi(A, all_b, xx);
        });
        
        double sp = tg / tlu;
        
        cout << left << setw(15) << k
             << setw(25) << fixed << setprecision(4) << tg
             << setw(25) << tlu
             << setw(20) << setprecision(2) << sp << "x" << endl;
    }
}

// эксперимент 3 

void exp3() {
    cout << "\n" << string(80, '=') << endl;
    cout << "ЭКСПЕРИМЕНТ 3: точность (матрица Гильберта)" << endl;
    cout << string(80, '=') << endl;
    
    vector<int> sizes = {5, 10, 15};
    
    cout << left << setw(10) << "n"
         << setw(25) << "Гаусс(без)"
         << setw(25) << "Гаусс(с)"
         << setw(25) << "LU" << endl;
    cout << string(85, '-') << endl;
    
    for (int n : sizes) {
        Matr H = hilbert_mat(n);
        Vec exact(n, 1.0);
        
        Vec b(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                b[i] += H[i][j] * exact[j];
            }
        }
        
        double r1 = -1, r2 = -1, r3 = -1;
        
        {
            Matr Hc = copy_mat(H);
            Vec bc = b;
            Vec x;
            if (gauss1(Hc, bc, x)) r1 = norm_resid(H, x, b);
        }
        
        {
            Matr Hc = copy_mat(H);
            Vec bc = b;
            Vec x;
            if (gauss2(Hc, bc, x)) r2 = norm_resid(H, x, b);
        }
        
        {
            Vec x;
            if (solve_lu(H, b, x)) r3 = norm_resid(H, x, b);
        }
        
        cout << left << setw(10) << n;
        
        if (r1 < 0) cout << setw(25) << "провал";
        else cout << setw(25) << scientific << setprecision(4) << r1;
        
        if (r2 < 0) cout << setw(25) << "провал";
        else cout << setw(25) << scientific << setprecision(4) << r2;
        
        if (r3 < 0) cout << setw(25) << "провал";
        else cout << setw(25) << scientific << setprecision(4) << r3;
        
        cout << endl;
    }
}

// ========================== эксперимент 4 ==========================

void exp4() {
    cout << "\n" << string(80, '=') << endl;
    cout << "эксперимент 4: масштабирование" << endl;
    cout << string(80, '=') << endl;
    
    vector<int> sizes = {50, 100, 200, 400, 600, 800, 1000};
    
    cout << left << setw(10) << "n"
         << setw(20) << "время (мс)"
         << setw(25) << "теор. O(n^3)" << endl;
    cout << string(55, '-') << endl;
    
    double base_t = 0;
    int base_n = 50;
    
    for (int n : sizes) {
        Matr A = rand_mat(n, 42);
        Vec b = rand_rhs(n, 123);
        
        double tt = time_me([&]() {
            Matr Ac = copy_mat(A);
            Vec bc = b;
            Vec x;
            gauss2(Ac, bc, x);
        });
        
        if (n == base_n) base_t = tt;
        
        double theo = base_t * pow(double(n) / base_n, 3.0);
        
        cout << left << setw(10) << n
             << setw(20) << fixed << setprecision(4) << tt
             << setw(25) << setprecision(1) << theo << endl;
    }
}

// ========================== main ==========================

int main() {
    cout << "ЛАБОРАТОРНАЯ РАБОТА: Сравнение методов решения СЛАУ" << endl;
    exp1();
    exp2();
    exp3();
    exp4();
    cout << "\n" << string(80, '=') << endl;
    cout << "готово" << endl;
    cout << string(80, '=') << endl;
    
    return 0;
}
