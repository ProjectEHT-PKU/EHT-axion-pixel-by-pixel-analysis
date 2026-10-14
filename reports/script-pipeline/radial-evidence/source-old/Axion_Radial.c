//
//  Axion_Radial.c
//  Radial_function
//
//

#include "Axion_Radial.h"

double get_max(double R[200])
{
    double max = 0.0;

    for (int n = 0; n < 200; n++)
    {
        if (R[n] > max)
        {
            max = R[n];
        }
    }
    return max;
}

double complex get_omega(double mu, double a)
{
    FILE *fp1 = NULL, *fp2 = NULL;
    double omega_mu_ima[500][2];
    double omega_mu_rea[500][2];
    if (FLAG_L_M == 11)
    {
        if (a == 0.99)
        {
            fp1 = fopen("source/a=0.99_l1m1_ima.txt", "r");
        }
        else if (a == 0.95)
        {
            fp1 = fopen("source/a=0.95_l1m1_ima.txt", "r");
        }
        else if (a == 0.9375)
        {
            fp1 = fopen("source/a=0.9375_l1m1_ima.txt", "r");
        }
        else if (a == 0.9)
        {
            fp1 = fopen("source/a=0.9_l1m1_ima.txt", "r");
        }
        else if (a == 0.8)
        {
            fp1 = fopen("source/a=0.8_l1m1_ima.txt", "r");
        }
        else if (a == 0.5)
        {
            fp1 = fopen("source/a=0.5_l1m1_ima.txt", "r");
        }
        else
        {
            printf("I'm wrong!!");
        }
    }
    if (FLAG_L_M == 22)
    {
        if (a == 0.99)
        {
            fp1 = fopen("source/a=0.99_l2m2_ima.txt", "r");
        }
        else if (a == 0.95)
        {
            fp1 = fopen("source/a=0.95_l2m2_ima.txt", "r");
        }
        else
        {
            printf("I'm wrong!!");
        }
    }
    int rows;
    fscanf(fp1, "%d", &rows);
    for (int n = 0; n < rows; n++)
    {
        fscanf(fp1, "%lf,", &omega_mu_ima[n][0]);
        fscanf(fp1, "%lf\n", &omega_mu_ima[n][1]);
    }
    fclose(fp1);
    if (omega_mu_ima[rows - 1][0] < mu)
    {
        printf("Sorry, no data for mu of this big!\n");
        return 0;
    }
    int i = 0;
    while (omega_mu_ima[i][0] < mu)
    {
        ++i;
    }
    --i;

    double omega_imag = ((omega_mu_ima[i + 1][1] - omega_mu_ima[i][1]) * (mu - omega_mu_ima[i][0]) / (omega_mu_ima[i + 1][0] - omega_mu_ima[i][0]) + omega_mu_ima[i][1]) * mu;

    if (FLAG_L_M == 11)
    {
        if (a == 0.99)
        {
            fp2 = fopen("source/a=0.99_l1m1_real.txt", "r");
        }
        else if (a == 0.95)
        {
            fp2 = fopen("source/a=0.95_l1m1_real.txt", "r");
        }
        else if (a == 0.9375)
        {
            fp2 = fopen("source/a=0.9375_l1m1_real.txt", "r");
        }
        else if (a == 0.9)
        {
            fp2 = fopen("source/a=0.9_l1m1_real.txt", "r");
        }
        else if (a == 0.8)
        {
            fp2 = fopen("source/a=0.8_l1m1_real.txt", "r");
        }
        else if (a == 0.5)
        {
            fp2 = fopen("source/a=0.5_l1m1_real.txt", "r");
        }
        else
        {
            printf("I'm wrong!!");
        }
    }
    if (FLAG_L_M == 22)
    {
        if (a == 0.99)
        {
            fp2 = fopen("source/a=0.99_l1m1_real.txt", "r");
        }
        else if (a == 0.95)
        {
            fp2 = fopen("source/a=0.95_l1m1_real.txt", "r");
        }
        else
        {
            printf("I'm wrong!!");
        }
    }
    int rows2;
    fscanf(fp2, "%d", &rows2);
    for (int n = 0; n < rows2; n++)
    {
        fscanf(fp2, "%lf,", &omega_mu_rea[n][0]);
        fscanf(fp2, "%lf\n", &omega_mu_rea[n][1]);
    }
    fclose(fp2);
    if (omega_mu_rea[rows2 - 1][0] < mu)
    {
        printf("Sorry, no data for mu of this big!\n");
        return 0;
    }
    int j = 0;
    while (omega_mu_rea[j][0] < mu)
    {
        ++j;
    }
    --j;

    double omega_real = ((omega_mu_rea[j + 1][1] - omega_mu_rea[j][1]) * (mu - omega_mu_rea[j][0]) / (omega_mu_rea[j + 1][0] - omega_mu_rea[j][0]) + omega_mu_rea[j][1]) * mu;

    double complex omega_complex = omega_real + I * omega_imag;
    return omega_complex;
}

struct Inter_variable get_Intermediate_variables(double M, double a, double complex omega, double mu)
{
    double complex b = csqrt(1 - cpow(a, 2.0));
    struct Inter_variable A;
    A.r_plus = M + csqrt(cpow(M, 2.0) - cpow(a, 2.0));
    A.r_minus = M - csqrt(cpow(M, 2.0) - cpow(a, 2.0));
    A.omega_c = (AZIMUTHAL_M * a) / (2.0 * M * A.r_plus);
    A.sigma = (2.0 * A.r_plus * (omega - A.omega_c)) / (A.r_plus - A.r_minus);
    A.q = -csqrt(cpow(mu, 2.0) - cpow(omega, 2.0));
    A.chi = (cpow(mu, 2.0) - 2.0 * cpow(omega, 2.0)) / A.q;
    A.c0 = 1.0 - 2.0 * I * omega - 2.0 * I * (omega - a * AZIMUTHAL_M / 2.0) / b;
    A.c1 = -4.0 + 4.0 * I * (omega - I * A.q * (1.0 + b)) + 4.0 * I * (omega - a * AZIMUTHAL_M / 2.0) / b - 2.0 * (cpow(omega, 2.0) + cpow(A.q, 2.0)) / A.q;
    A.c2 = 3.0 - 2.0 * I * omega - 2.0 * (cpow(A.q, 2.0) - cpow(omega, 2.0)) / A.q - 2.0 * I * (omega - a * AZIMUTHAL_M / 2.0) / b;
    A.c3 = 2.0 * I * cpow((omega - I * A.q), 3.0) / A.q + 2.0 * cpow((omega - I * A.q), 2.0) * b + cpow(A.q * a, 2.0) + 2.0 * I * A.q * a * AZIMUTHAL_M - LAMBDA_L_M - 1.0 - cpow((omega - I * A.q), 2.0) / A.q + 2.0 * A.q * b + 2.0 * I * (cpow((omega - I * A.q), 2.0) / A.q + 1.0) * (omega - a * AZIMUTHAL_M / 2.0) / b;
    A.c4 = cpow((omega - I * A.q), 4.0) / cpow(A.q, 2.0) + 2 * I * omega * cpow((omega - I * A.q), 2.0) / A.q - 2.0 * I * (cpow((omega - I * A.q), 2.0) / A.q) * (omega - a * AZIMUTHAL_M / 2.0) / b;
    A.power1 = -1.0 * I * A.sigma;
    A.power2 = I * A.sigma + A.chi - 1.0;
    return A;
}

void change_coefficient(int n, double complex *a, double complex *a_previous, struct Inter_variable A) // change the circulation coefficient to the next round
{
    double complex alpha = cpow(n, 2.0) + (A.c0 + 1.0) * n + A.c0;
    double complex beta = cpow(n, 2.0) * (-2.0) + (A.c1 + 2.0) * n + A.c3;
    double complex gamma = cpow(n, 2.0) + (A.c2 - 3.0) * n + A.c4;
    double complex a_n = *a;            // a_n
    double complex a_n_1 = *a_previous; // a_(n-1)
    *a = -(beta * a_n + gamma * a_n_1) / alpha;
    *a_previous = a_n;
}

double get_radius(double X[4])
{
    double r = X[1];
    return r;
}

double complex Circulation_R(double r, struct Inter_variable A) // the circulation part of R(r)
{
    double complex a, a_previous; // a is the coefficients of the circulation a_n.And a_previous is a_(n-1)
    a = CIRCULATION_A_0;
    a_previous = 0.0;
    double complex result = 0.0;
    for (int n = 0; n < CIRCLE_TIMES; n++)
    {
        result += a * cpow((r - A.r_plus) / (r - A.r_minus), n); // caculate the nth term of the circulation

        change_coefficient(n, &a, &a_previous, A); // make the coefficient a_n changed into a_(n+1)
    }
    return result;
}

double Function_R(double X[4], double a, double mu)
{
    double M = 1.0;
    double omega = get_omega(mu, a);
    struct Inter_variable Inter_var = get_Intermediate_variables(M, a, omega, mu);
    double r = get_radius(X);
    double complex result, result1;
    double complex circluation = Circulation_R(r, Inter_var);
    result = cpow((r - Inter_var.r_plus), Inter_var.power1) * cpow(r - Inter_var.r_minus, Inter_var.power2) * cpow(NUMBER_OF_E, (Inter_var.q * r)) * circluation;
    result1 = cabs(result);

    return result1;
}

void print_cable(double a, double Mu_min, double Mu_max)
{
    double X[4];
    double Result[40][1200][3];
    double complex omega_complex, R_r;
    double max_temp, mu;
    double temp[1200];
    for (int i = 0; i < 40; i++)
    {
        printf("i=%d\n", i);
        mu = (Mu_max - Mu_min) / 39.0 * i + Mu_min;
        Result[i][0][0] = mu;
        Result[i][0][1] = 1.0;
        Result[i][0][2] = 0.0;
        for (int j = 1; j < 1200; j++)
        {
            X[1] = (600 - 2) / 598 * (j - 1) + 2;
            omega_complex = get_omega(mu, a);
            R_r = Function_R(X, a, mu);
            Result[i][j][0] = mu;
            Result[i][j][1] = X[1];
            Result[i][j][2] = R_r;
            temp[j] = Result[i][j][2];
        }
        /*
        max_temp = get_max(temp);
        for(int j=0;j<1200;j++)
        {
            Result[i][j][2] = temp[j]/max_temp;
            temp[j]=0;
        }*/
    }
    FILE *result;
    if (FLAG_L_M == 11)
    {
        if (a == 0.99)
        {
            result = fopen("result/50.txt", "w");
        }
        else if (a == 0.95)
        {
            result = fopen("result/result_a=0.95_l1m1.txt", "w");
        }
        else if (a == 0.9375)
        {
            result = fopen("result/result_a=0.9375_l1m1.txt", "w");
        }
        else if (a == 0.9)
        {
            result = fopen("result/result_a=0.9_l1m1.txt", "w");
        }
        else if (a == 0.8)
        {
            result = fopen("result/result_a=0.8_l1m1.txt", "w");
        }
        else if (a == 0.5)
        {
            result = fopen("result/result_a=0.5_l1m1.txt", "w");
        }
        else
        {
            printf("I'm wrong when printing!!");
        }
    }
    if (FLAG_L_M == 22)
    {
        if (a == 0.99)
        {
            result = fopen("result/result_a=0.99_l2m2.txt", "w");
        }
        else if (a == 0.95)
        {
            result = fopen("result/result_a=0.95_l2m2.txt", "w");
        }
        else
        {
            printf("I'm wrong in printing!!");
        }
    }
    for (int i = 0; i < 40; i++)
    {
        for (int j = 0; j < 1200; j++)
        {
            fprintf(result, "%f,%f,%f\n", Result[i][j][0], Result[i][j][1], Result[i][j][2]);
        }
    }
    fclose(result);
    printf("Mission print complete!");
    return;
}

void print_all_cable()
{
    double mu = 0.4;
    if (FLAG_L_M == 11)
    {
        double a = 0.99;
        double Mu_min = 0.1;
        double Mu_max = 0.446;
        print_cable(a, Mu_min, Mu_max);
        a = 0.95;
        Mu_max = 0.367;
        print_cable(a, Mu_min, Mu_max);
        a = 0.9375;
        Mu_max = 0.316;
        print_cable(a, Mu_min, Mu_max);
        a = 0.9;
        Mu_max = 0.316;
        print_cable(a, Mu_min, Mu_max);
        a = 0.8;
        Mu_max = 0.249;
        print_cable(a, Mu_min, Mu_max);
        a = 0.5;
        Mu_max = 0.132;
        print_cable(a, Mu_min, Mu_max);
    }
    if (FLAG_L_M == 22)
    {
        double a = 0.99;
        double Mu_min = 0.320;
        double Mu_max = 0.928;
        print_cable(a, Mu_min, Mu_max);
        a = 0.95;
        Mu_max = 0.752;
        print_cable(a, Mu_min, Mu_max);
    }
}
