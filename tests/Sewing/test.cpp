    #include <iostream>
    #include <Eigen/Core>
    #include <autodiff/reverse.hpp>

    using namespace std;
    using namespace autodiff;

    var f1(var x)
    {
        return x*x;
    }

    var f2(var x)
    {
        return pow(x, 2);
    }

    int main(){
        var x = 0.0;                           // zero is special case but when power >= 1 derivative should exist
        var u1 = f1(x);                        // f1 = x * x
        var u2 = f2(x);                        // f2 = pow(x, 2)

        Derivatives dud1 = derivatives(u1);
        Derivatives dud2 = derivatives(u2);

        var dudx1 = dud1(x);
        var dudx2 = dud2(x);

        cout << "du/dx1 = " << dudx1 << endl;  // here du/dx = 0, that's OK
        cout << "du/dx2 = " << dudx2 << endl;  // but here du/dx = -nan, that's possible bug
    }