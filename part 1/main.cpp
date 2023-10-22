#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>

const double epsilon = 1e-4;

//34
double f(double x, double alpha)
{
    double tmp = std::atan(std::exp(-x*x));
    return std::sin(1/(1+tmp*tmp)) + x - alpha;
}

double phi(double x, double alpha)
{
    return x - f(x, alpha);
}

//delta = m/(M-m)*epsilon
double secant(double (*f)(double, double), double x0, double alpha, double a, double b, double delta)
{
    std::ofstream file("output.txt", std::ios::app);
    double x = x0;
    do
    {
        x0 = x;
        x = x0 - f(x0, alpha)/(f(x0, alpha) - f(a, alpha))*(x0-a);
        file << x << " " << f(x, alpha) << "\n";
    } while (std::abs(x-x0) >= delta);
    file << "\n";
    file.close();
    return x;
}

//delta = epsilon
double bisection(double (*f)(double, double), double x0, double alpha, double a, double b, double delta)
{
    std::ofstream file("output.txt", std::ios::app);
    double x = (a+b)/2;
    bool sa = f(a, alpha) > 0, s;
    int N = std::floor(std::log((b-a)/delta)/std::log(2))+1;
    for (int i = 0; i < N; i++)
    {
        s = (f(x, alpha) > 0);
        if (s ^ sa)
        {
            b = x;
        }
        else
        {
            a = x;
            sa = s;
        }
        x = (a+b)/2;
        file << x << " " << f(x, alpha) << "\n";
    }
    file << "\n";
    file.close();
    return x;
}

//delta = (1-q)/q*epsilon
double simple_iter(double (*phi)(double, double), double x0, double alpha, double delta)
{
    std::ofstream file("output.txt", std::ios::app);
    double x = x0;
    do
    {
        x0 = x;
        x = phi(x0, alpha);
        file << x << " " << f(x, alpha) << "\n";
    } while (std::abs(x-x0) >= delta);
    file << "\n";
    file.close();
    return x;
}

//delta = m/(M-m)*epsilon
double newton(double (*f)(double, double), double x0, double alpha, double delta)
{
    std::ofstream file("output.txt", std::ios::app);
    const double h = std::cbrt(std::numeric_limits<double>::epsilon());
    double x = x0;
    do
    {
        x0 = x;
        x = x0 - f(x0, alpha)/(f(x0+h, alpha) - f(x0-h, alpha))*2*h;
        file << x << " " << f(x, alpha) << "\n";
    } while (std::abs(x-x0) >= delta);
    file << "\n";
    file.close();
    return x;
}

//delta = (1-q)/q*epsilon
double aitken(double (*phi)(double, double), double x0, double alpha, double delta)
{
    std::ofstream file("output.txt", std::ios::app);
    double x = x0, x1, x2;
    do
    {
        x0 = x;
        x1 = phi(x0, alpha);
        x2 = phi(x1, alpha);
        x = (x2*x0-x1*x1)/(x2-2*x1+x0);
        file << x << " " << f(x, alpha) << "\n";
    } while (std::abs(x-x0) >= delta);
    file << "\n";
    file.close();
    return x;
}

int main()
{
    const double a = 0, b = 5, m = 1, M = 1.3, q = 0.3;
    double x0 = (a+b)/2;
    std::ofstream file("output.txt");
    file.close();
    for (int alpha = 1; alpha < 6; alpha++)
    {
        secant(f, x0, 1, a, b, m/(M-m)*epsilon);
        bisection(f, x0, 1, a, b, epsilon);
        simple_iter(phi, x0, 1, (1-q)/q*epsilon);
        newton(f, x0, 1, m/(M-m)*epsilon);
        aitken(phi, x0, 1, (1-q)/q*epsilon);
    }
    return 0;
}