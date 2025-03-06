#include <iostream>
#include <Windows.h>
#include <cmath>
using namespace std;

double epsilon1 = 0.01;  //точности решений
double epsilon2 = 0.01;
double fi = 10;  //значение штрафной функции в текущей точке
double r = 10;  //штрафной коэффициент
double C = 10;  //число для увеличения параметра штрафа
int k = 0;  //номер итерации

struct Point
{
    double x1;
    double x2;
};

double function(Point x) {  //исходная функция
    return 2 * pow(x.x1 - 1, 2) + 3 * pow(x.x2 - 3, 2);
}

double H(Point x) {  //равенство
    return x.x1 + x.x2 - 6;
}

double G(Point x) { //неравенство
    return x.x1 * x.x1 + x.x2 * x.x2 - 10;
}

double F(Point x, double r) {  //штрафная функция
    return (H(x) * H(x) / (2 * r) + r / G(x));
}

double P(Point x, double r) {  //вспомогательная функция
    return function(x) + F(x, r);
}

Point gradientFunction(Point x) {
    return{ 4 * (x.x1 - 1), 6 * (x.x2 - 3) };
}

Point gradientH(Point x) {
    return{ 1, 1};
}

Point gradientG(Point x) {
    return{ 2 * x.x1, 2 * x.x2 };
}

Point Shtrafgradient(Point x, double r) // градиент штрафной функции
{
    Point gFi;
    gFi.x1 = H(x) * gradientH(x).x1 / r - 2 * r * gradientG(x).x1 / (G(x) * G(x));
    gFi.x2 = H(x) * gradientH(x).x2 / r - 2 * r * gradientG(x).x2 / (G(x) * G(x));
    return gFi;
}

Point gradientP(Point x, double r)
{
    return{ gradientFunction(x).x1 + Shtrafgradient(x, r).x1, gradientFunction(x).x2 + Shtrafgradient(x, r).x2 };
}

double Norma(Point x)
{
    return sqrt(x.x1 * x.x1 + x.x2 * x.x2);
}

Point GradientSpusk(Point p0, double r)
{
    Point x, gf = gradientP(p0, r);
    int k_local = 0;
    double t = 1.0;
    do
    {
        k_local++;

        x.x1 = p0.x1 - t * gf.x1;
        x.x2 = p0.x2 - t * gf.x2;

        if ((P(x, r) - P(p0, r)) >= 0) t = t / 2;
        p0 = x;
        gf = gradientP(p0, r);

    } while (Norma(gf) > epsilon2);

    return x;
}
void IterationProcess(Point& p0, Point& x) {
    do {
        k++;
        x = GradientSpusk(p0, r);
        fi = F(x, r);
        r = r / C;
        p0 = x; // Обновление начальной точки для следующей итерации
    } while (abs(fi) > epsilon1);
}

int main()
{
    setlocale(LC_ALL, "rus");
    Point p0, x;

    p0.x1 = 1, p0.x2 = 1; // Начальное приближение

    IterationProcess(p0, x);

    cout << "Точка условного экстремума =  " << round(x.x1) << ", " << round(x.x2) << "\n";
    cout << "Количество итераций K = " << k << endl;
    return 0;
}
