#include <iostream>
#include <Windows.h>
#include <cmath>
using namespace std;

double epsilon1 = 0.01;
double epsilon2 = 0.01;
double fi = 10;  //значение штрафной функции в текущей точке
double r = 10;  //штрафной коэффициент
double C = 10;  //число для увеличения параметра штрафа
int k = 0;

struct Point
{
    double x1;
    double x2;
};

double function(Point X) {
    return 2 * pow(X.x1 - 1, 2) + 3 * pow(X.x2 - 3, 2);
}

double H(Point X) {  //равенство
    return X.x1 + X.x2 - 6;
}

double G(Point X) { //неравенство
    return X.x1 * X.x1 + X.x2 * X.x2 - 10;
}

double F(Point X, double r) {  //штрафная функция
    return (H(X) * H(X) / (2 * r) + r / G(X)); // Исправлено знак в штрафной функции
}

double P(Point X, double r) {  //вспомогательная функция
    return function(X) + F(X, r);
}

Point gradientF(Point X) {
    return{ 4 * (X.x1 - 1), 6 * (X.x2 - 3) };
}

Point gradientH(Point X) {
    return{ 1, 1 };
}

Point gradientG(Point X) {
    return{ 2 * X.x1, 2 * X.x2 };
}

Point Shtrafgradient(Point X, double r) // градиент штрафной функции
{
    Point gFi;
    gFi.x1 = H(X) * gradientH(X).x1 / r - 2 * r * gradientG(X).x1 / (G(X) * G(X)); // Исправлено знак
    gFi.x2 = H(X) * gradientH(X).x2 / r - 2 * r * gradientG(X).x2 / (G(X) * G(X)); // Исправлено знак
    return gFi;
}

Point gradientP(Point X, double r)
{
    return{ gradientF(X).x1 + Shtrafgradient(X, r).x1, gradientF(X).x2 + Shtrafgradient(X, r).x2 };
}

double Norma(Point X)
{
    return sqrt(X.x1 * X.x1 + X.x2 * X.x2);
}

Point GradientSpusk(Point p0, double r)
{
    Point X, gf = gradientP(p0, r);
    int k_local = 0; // Исправлено имя переменной для избежания конфликтов
    double t = 1.0;
    do
    {
        k_local++;

        X.x1 = p0.x1 - t * gf.x1;
        X.x2 = p0.x2 - t * gf.x2;

        if ((P(X, r) - P(p0, r)) >= 0) t = t / 2;
        p0 = X;
        gf = gradientP(p0, r);

    } while (Norma(gf) > epsilon2);

    return X;
}

int main()
{
    setlocale(LC_ALL, "rus");
    Point p0, X;

    p0.x1 = 1, p0.x2 = 1; // Начальные значения

    do
    {
        k++;
        X = GradientSpusk(p0, r);
        fi = F(X, r);
        r = r / C;
        p0 = X; // Обновление начальной точки для следующей итерации
    } while (abs(fi) > epsilon1);

    cout << "Точка условного экстремума =  " << round(X.x1) << ", " << round(X.x2) << "\n";
    cout << "Количество итераций K = " << k << endl;
    //cout << "function = " << round(function(X)) << endl;
    //cout << "Проверка:\n" << "h = " << H(X) << endl;
    //cout << "g = " << G(X) << endl;

    return 0;
}
