#include <fstream>
#include <iostream>
#include <cmath>
#include <array>

const double epsilon = 1e-4;

std::array<double, 2> f(std::array<double, 2> x)
{
    std::array<double, 2> result = {x[0]*x[0]+x[1]*x[1]-2*x[0], x[0]-std::exp(-x[1])};
    return result;
}

std::array<double, 2> phi(std::array<double, 2> x)
{
    std::array<double, 2> tmp = f(x);
    std::array<double, 2> result = {x[0]-tmp[1], x[1]-tmp[0]};
    return result;
}

double phi0(std::array<double, 2> x)
{
    return x[0] - (x[0] - std::exp(-x[1]))/2;
}

double phi1(std::array<double, 2> x)
{
    return x[1]-(x[0]*x[0]+x[1]*x[1]-2*x[0])/2;
}

double norm(std::array<double, 2> x)
{return std::max(std::abs(x[0]), std::abs(x[1]));}

double norm(std::array<double, 2> x, std::array<double, 2> x0)
{return std::max(std::abs(x[0]-x0[0]), std::abs(x[1]-x0[1]));}

//delta = (1-q)/q*epsilon
std::array<double, 2> simple_iter(std::array<double, 2> (*phi)(std::array<double, 2>), std::array<double, 2> x0, double delta)
{
    std::ofstream file("output.txt", std::ios::app);
    std::array<double, 2> x = x0;
    do
    {
        x0 = x;
        x = phi(x0);
        file << x[0] << " " << x[1] << " " << norm(f(x)) << "\n";
    } while (norm(x, x0) >= delta);
    file << "\n";
    file.close();
    return x;
}

//delta = epsilon/M
std::array<double, 2> newton(std::array<double, 2> (*f)(std::array<double, 2>), std::array<double, 2> x0, double delta)
{
    std::ofstream file("output.txt", std::ios::app);
    const double h = std::cbrt(std::numeric_limits<double>::epsilon());
    std::array<double, 2> x = x0, xp1, xm1, rp1, rm1, r;
    std::array<std::array<double, 2>, 2> jac, invjac;
    double det;
    do
    {
        x0 = x;
        xp1 = x0; xm1 = x0;
        xp1[0] += h; xm1[0] -= h; rp1 = f(xp1); rm1 = f(xm1);
        jac[0][0] = (rp1[0]-rm1[0])/h/2; jac[1][0] = (rp1[1]-rm1[1])/h/2;
        xp1 = x0; xm1 = x0;
        xp1[1] += h; xm1[1] -= h; rp1 = f(xp1); rm1 = f(xm1);
        jac[0][1] = (rp1[0]-rm1[0])/h/2; jac[1][1] = (rp1[1]-rm1[1])/h/2;
        det = jac[0][0]*jac[1][1]-jac[1][0]*jac[0][1];
        invjac[0][0] = jac[1][1]/det; invjac[1][1] = jac[0][0]/det; invjac[1][0] = -jac[1][0]/det; invjac[0][1] = -jac[0][1]/det;
        r = f(x);
        x[0] = x0[0] - r[0]*invjac[0][0] - r[1]*invjac[0][1];
        x[1] = x0[1] - r[0]*invjac[1][0] - r[1]*invjac[1][1];
        file << x[0] << " " << x[1] << " " << norm(f(x)) << "\n";
    } while (norm(x, x0) >= delta);
    file << "\n";
    file.close();
    return x;
}

//delta = (1-q)/q*epsilon
std::array<double, 2> seidel(double (*phi0)(std::array<double, 2>), double (*phi1)(std::array<double, 2>), std::array<double, 2> x0, double delta)
{
    std::ofstream file("output.txt", std::ios::app);
    std::array<double, 2> x = x0;
    do
    {
        x0 = x;
        x[0] = phi0(x);
        x[1] = phi1(x);
        file << x[0] << " " << x[1] << " " << norm(f(x)) << "\n";
    } while (norm(x, x0) >= delta);
    file << "\n";
    file.close();
    return x;
}

int main()
{
    const std::array<double, 2> a = {0.4, 0.8}, b = {0.5, 0.85};
    const double M = 1.3, q = 0.3;
    std::array<double, 2> x0 = {(a[0]+b[0])/2, (a[1]+b[1])/2};
    std::ofstream file("output.txt");
    file.close();
    simple_iter(phi, x0, (1-q)/q*epsilon);
    newton(f, x0, epsilon/M);
    seidel(phi0, phi1, x0, (1-q)/q*epsilon);
    return 0;
}