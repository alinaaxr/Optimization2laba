#include <iostream>
#include <Windows.h>
#include <cmath>
using namespace std;

double epsilon1 = 0.01;  //точности решений
double epsilon2 = 0.01;
double fi = 10;  //значение штрафной функции в текущей точке
double r = 10;  //штрафной коэффициент
double C = 12;  //число для уменьшения параметра штрафа
int k = 0;  //номер итерации

struct Point
{
    double x1;
    double x2;
    double x3;
};

double function(Point x) {  //исходная функция
    return x.x1 * x.x2 * x.x3;
}

double H(Point x) {  //равенство
    return x.x1 * x.x1 + x.x2 * x.x2 + x.x3 * x.x3 - 1;
}

double G(Point x) { //неравенство
    return x.x1 + x.x2 + x.x3;
}

double F(Point x, double r) {  //штрафная функция с проверкой
    double hPart = H(x) * H(x) / (2 * r);

    double gPart;
    if (abs(G(x)) < epsilon1) {
        gPart = 0; // Если G(x) близко к нулю, штраф не добавляем
    }
    else if (G(x) > 0) {
        gPart = r / G(x); // Используем r/G(x) для положительных G(x)
    }
    else {
        gPart = 0; // Для отрицательных G(x) штраф не добавляем
    }

    return hPart + gPart;
}

double P(Point x, double r) {  //вспомогательная функция
    return function(x) + F(x, r);
}

Point gradientFunction(Point x) {
    return{ x.x2 * x.x3, x.x1 * x.x3, x.x1 * x.x2 };
}

Point gradientH(Point x) {
    return{ 2 * x.x1, 2 * x.x2, 2 * x.x3 };
}

Point gradientG(Point x) {
    return{ 1, 1, 1 };
}

Point Shtrafgradient(Point x, double r) // градиент штрафной функции
{
    Point gFi;

    double hGrad = H(x);
    double gGrad;

    if (abs(G(x)) < epsilon1) {
        gGrad = 0;
    }
    else if (G(x) > 0) {
        gGrad = -r / (G(x) * G(x)); // Производная r/G(x)
    }
    else {
        gGrad = 0;
    }

    gFi.x1 = hGrad * gradientH(x).x1 / r + gGrad * gradientG(x).x1;
    gFi.x2 = hGrad * gradientH(x).x2 / r + gGrad * gradientG(x).x2;
    gFi.x3 = hGrad * gradientH(x).x3 / r + gGrad * gradientG(x).x3;

    return gFi;
}

Point gradientP(Point x, double r)
{
    Point gf = gradientFunction(x);
    Point sgf = Shtrafgradient(x, r);
    return{ gf.x1 + sgf.x1, gf.x2 + sgf.x2, gf.x3 + sgf.x3 };
}

double Norma(Point x)
{
    return sqrt(x.x1 * x.x1 + x.x2 * x.x2 + x.x3 * x.x3);
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
        x.x3 = p0.x3 - t * gf.x3;

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
    } while (abs(H(x)) > epsilon1 && abs(max(0.0, G(x))) > epsilon1);
}

int main()
{
    setlocale(LC_ALL, "rus");
    Point p0, x;

    p0.x1 = -0.5, p0.x2 = -0.5, p0.x3 = -0.5; // Начальное приближение

    IterationProcess(p0, x);

    cout << "Точка условного экстремума =  " << round(x.x1 * 100) / 100 << ", " << round(x.x2 * 100) / 100.0 << ", " << round(x.x3 * 100) / 100.0 << "\n";
    cout << "Значение функции в точке экстремума = " << function(x) << endl;
    cout << "Количество итераций K = " << k << endl;
    return 0;
}
