#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return 2 * x - 1 - sin(x);
}

int main() {
    double a = 0;
    double b = 1;
    double tol = 1e-10;

    int iterations = 0;
    
    while ((b-a) / 2 > tol) {
        double c = (a+b) / 2;

        if (f(c) == 0) {
            break;
        }
        else if  (f(a) * f(c) < 0){
            b = c;
        }
        else {
            a = c;
        }

        iterations++;
    }

    double root = (a+b) / 2;

    cout << "Root: " << root << endl;
    cout << "Iterations: " << iterations << endl;

    return 0;
}