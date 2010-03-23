#include <iostream>
#include <stdexcept>

#include "matrix.h"
#include "python_solvers.h"

#define EPS 1e-12

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_solver1()
{
    CooMatrix A(4);
    A.add(0, 0, -1);
    A.add(1, 1, -1);
    A.add(2, 2, -1);
    A.add(3, 3, -1);
    A.add(0, 1, 2);
    A.add(1, 0, 2);
    A.add(1, 2, 2);
    A.add(2, 1, 2);
    A.add(2, 3, 2);
    A.add(3, 2, 2);

    double res[4] = {1., 1., 1., 1.};

    solve_linear_system_dense_lu(&A, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);
}

void test_solver2()
{
    CooMatrix A(4);
    A.add(0, 0, -1);
    A.add(1, 1, -1);
    A.add(2, 2, -1);
    A.add(3, 3, -1);
    A.add(0, 1, 2);
    A.add(1, 0, 2);
    A.add(1, 2, 2);
    A.add(2, 1, 2);
    A.add(2, 3, 2);
    A.add(3, 2, 2);

    double res[4] = {1., 1., 1., 1.};

    solve_linear_system_dense_lu(&A, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);

    DenseMatrix B(&A);
    for (int i=0; i < 4; i++) res[i] = 1.;

    solve_linear_system_dense_lu(&B, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);
}

void test_solver3()
{
    CooMatrix A(4);
    A.add(0, 0, -1);
    A.add(1, 1, -1);
    A.add(2, 2, -1);
    A.add(3, 3, -1);
    A.add(0, 1, 2);
    A.add(1, 0, 2);
    A.add(1, 2, 2);
    A.add(2, 1, 2);
    A.add(2, 3, 2);
    A.add(3, 2, 2);

    double res[4] = {1., 1., 1., 1.};

    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);

    for (int i=0; i < 4; i++) res[i] = 1.;
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);

    for (int i=0; i < 4; i++) res[i] = 1.;
    solve_linear_system_scipy_cg(&A, res);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);

    for (int i=0; i < 4; i++) res[i] = 1.;
    _assert(solve_linear_system_cg(&A, res, EPS, 2) == 1);
    _assert(fabs(res[0] - 0.2) < EPS);
    _assert(fabs(res[1] - 0.6) < EPS);
    _assert(fabs(res[2] - 0.6) < EPS);
    _assert(fabs(res[3] - 0.2) < EPS);
}

void test_solver4()
{
    CooMatrix A(5);
    A.add(0, 0, 2);
    A.add(0, 1, 3);
    A.add(1, 0, 3);
    A.add(1, 2, 4);
    A.add(1, 4, 6);
    A.add(2, 1, -1);
    A.add(2, 2, -3);
    A.add(2, 3, 2);
    A.add(3, 2, 1);
    A.add(4, 1, 4);
    A.add(4, 2, 2);
    A.add(4, 4, 1);

    double res[5] = {8., 45., -3., 3., 19.};
    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0] - 1.) < EPS);
    _assert(fabs(res[1] - 2.) < EPS);
    _assert(fabs(res[2] - 3.) < EPS);
    _assert(fabs(res[3] - 4.) < EPS);
    _assert(fabs(res[4] - 5.) < EPS);

    res[0] = 8.;
    res[1] = 45.;
    res[2] = -3.;
    res[3] = 3.;
    res[4] = 19.;
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0] - 1.) < EPS);
    _assert(fabs(res[1] - 2.) < EPS);
    _assert(fabs(res[2] - 3.) < EPS);
    _assert(fabs(res[3] - 4.) < EPS);
    _assert(fabs(res[4] - 5.) < EPS);

    res[0] = 8.;
    res[1] = 45.;
    res[2] = -3.;
    res[3] = 3.;
    res[4] = 19.;
    solve_linear_system_scipy_gmres(&A, res);
    _assert(fabs(res[0] - 1.) < EPS);
    _assert(fabs(res[1] - 2.) < EPS);
    _assert(fabs(res[2] - 3.) < EPS);
    _assert(fabs(res[3] - 4.) < EPS);
    _assert(fabs(res[4] - 5.) < EPS);
}

void test_solver5()
{
    CooMatrix A(5);
    A.add(0, 0, 2);
    A.add(0, 1, 3);
    A.add(1, 0, 3);
    A.add(1, 2, 4);
    A.add(1, 4, 6);
    A.add(2, 1, -1);
    A.add(2, 2, -3);
    A.add(2, 3, 2);
    A.add(3, 2, 1);
    A.add(4, 1, 4);
    A.add(4, 2, 2);
    A.add(4, 4, 1);

    cplx res[5];

    res[0] = cplx(8.);
    res[1] = cplx(45.);
    res[2] = cplx(-3.);
    res[3] = cplx(3.);
    res[4] = cplx(19.);
    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0].real() - 1.) < EPS);
    _assert(fabs(res[1].real() - 2.) < EPS);
    _assert(fabs(res[2].real() - 3.) < EPS);
    _assert(fabs(res[3].real() - 4.) < EPS);
    _assert(fabs(res[4].real() - 5.) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - 0.) < EPS);
    _assert(fabs(res[2].imag() - 0.) < EPS);
    _assert(fabs(res[3].imag() - 0.) < EPS);
    _assert(fabs(res[4].imag() - 0.) < EPS);

    res[0] = cplx(8.);
    res[1] = cplx(45.);
    res[2] = cplx(-3.);
    res[3] = cplx(3.);
    res[4] = cplx(19.);
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0].real() - 1.) < EPS);
    _assert(fabs(res[1].real() - 2.) < EPS);
    _assert(fabs(res[2].real() - 3.) < EPS);
    _assert(fabs(res[3].real() - 4.) < EPS);
    _assert(fabs(res[4].real() - 5.) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - 0.) < EPS);
    _assert(fabs(res[2].imag() - 0.) < EPS);
    _assert(fabs(res[3].imag() - 0.) < EPS);
    _assert(fabs(res[4].imag() - 0.) < EPS);
}

void test_solver6()
{
    CooMatrix A(2, true);
    A.add(0, 0, cplx(1, 1));
    A.add(0, 1, cplx(2, 2));
    A.add(1, 0, cplx(3, 3));
    A.add(1, 1, cplx(4, 4));

    cplx res[2];

    //----------------------

    res[0] = cplx(1);
    res[1] = cplx(2);
    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0].real() - 0) < EPS);
    _assert(fabs(res[1].real() - 0.25) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - (-0.25)) < EPS);

    res[0] = cplx(1);
    res[1] = cplx(2);
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0].real() - 0) < EPS);
    _assert(fabs(res[1].real() - 0.25) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - (-0.25)) < EPS);

    //----------------------

    res[0] = cplx(1, 1);
    res[1] = cplx(2, 2);
    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0].real() - 0) < EPS);
    _assert(fabs(res[1].real() - 0.5) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - 0.) < EPS);

    res[0] = cplx(1, 1);
    res[1] = cplx(2, 2);
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0].real() - 0) < EPS);
    _assert(fabs(res[1].real() - 0.5) < EPS);
    _assert(fabs(res[0].imag() - 0.) < EPS);
    _assert(fabs(res[1].imag() - 0.) < EPS);

    //----------------------

    res[0] = cplx(2, 1);
    res[1] = cplx(2, 2);
    solve_linear_system_numpy(&A, res);
    _assert(fabs(res[0].real() - (-1)) < EPS);
    _assert(fabs(res[1].real() - 1.25) < EPS);
    _assert(fabs(res[0].imag() - 1.) < EPS);
    _assert(fabs(res[1].imag() - (-0.75)) < EPS);

    res[0] = cplx(2, 1);
    res[1] = cplx(2, 2);
    solve_linear_system_scipy_umfpack(&A, res);
    _assert(fabs(res[0].real() - (-1)) < EPS);
    _assert(fabs(res[1].real() - 1.25) < EPS);
    _assert(fabs(res[0].imag() - 1.) < EPS);
    _assert(fabs(res[1].imag() - (-0.75)) < EPS);
}

int main(int argc, char* argv[])
{
    try {
        test_solver1();
        test_solver2();
        test_solver3();
        test_solver4();
        test_solver5();
        test_solver6();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
