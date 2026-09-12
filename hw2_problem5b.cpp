#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return x*x*x + x - 4;
}

int main() {
    double a = 1.0;
    double b = 4.0;
    double tol = 1e-3;

    int n_bound = 0;
    double error_bound = b - a;

    while (error_bound > tol) {
        error_bound /= 2.0;
        n_bound++;
    }

    cout << "Theoretical upper bound: " << n_bound << " iterations\n";
    cout << "Error bound after " << n_bound << " iterations: "
         << error_bound << "\n\n";

    int iterations = 0;

    while (true) {
        double d = (a + b) / 2.0;

        double current_bound = (b - a) / 2.0;

        iterations++;

        if (current_bound <= tol) {
            break;
        }

        if (f(a) * f(d) < 0) {
            b = d;
        } else {
            a = d;
        }
    }

    double root = (a + b) / 2.0;

    cout << "Bisection approximation: " << root << "\n";
    cout << "Iterations used: " << iterations << "\n";

    return 0;
}