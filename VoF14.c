// Volume of Fluids Algorithm with Piecewise Linear Interface Construction (PLIC)
// Drake Jones
// San Diego State University
// Thesis Title: Machine Accurate Mass Conserving Interface Capturing Methods For Two Phase Flows
// Fall 2023

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "dimc.h"
#include "cvoftools.h"
#include "cmesh.h"
#include "cuservoftools.h"

double *VoF_tools(double C, double *m_i, double dx0, int m, int n)
{
    // using analytical tools developed by Lopez and Hernandez to solve plane equation
    // (n_x)(x) + (n_y)(y) + alpha = 0 and calculate interface segments of cell: (x1,y1), (x2,y2)
    // output rows: X = [x1; x2; y1; y2]

    int ntp, ntv, icontn, icontp, ipv[nv], cc, ym, ij, cnt, num;
    double v, vertp[nv * 2], vt, a, n_x, n_y, mag;

    double *A;
    A = malloc(sizeof(double *) * 6);

    mag = sqrt(m_i[0] * m_i[0] + m_i[1] * m_i[1]); // magnitude of "normalized" normal vector

    n_x = m_i[0] / mag; // conventional normalizing of interface normal
    n_y = m_i[1] / mag;

    v = C * dx0 * dx0; // area (v0lume) of reference fluid in cell i,j

    squaremesh_(ipv, &ntp, &ntv, vertp);
    enforv2dsz_(&a, &dx0, &dx0, &v, vertp, &n_x, &n_y);                // solves for alpha
    inte2d_(&a, &icontn, &icontp, ipv, &ntp, &ntv, vertp, &n_x, &n_y); // outputs local cell coordinates

    // output matrix: A = [x1; x2; y1; y2; n_x; n_y]
    A[0] = vertp[4];
    A[1] = vertp[5];
    A[2] = vertp[124];
    A[3] = vertp[125];
    A[4] = n_x;
    A[5] = n_y;

    return A;
}

double error_el(double **C_g, double *A, double dx0, int i, int j, int ecase)
{

    // A input values: A = [x1 x2 y1 y2 n_x n_y] from PLIC cell
    int ip, jp, n, nx_l;
    double m, b, err_el, Cs, n_x, n_y, y1, y2, dx_l, Q, x_l, y_l, E_l;

    nx_l = 16;         // number of local grid cells
    dx_l = 1.0 / nx_l; // local grid spacing
    err_el = 0.0;      // initialzing error

    double dif[nx_l + 1];

    n_x = A[4];
    n_y = A[5]; // interface normal components

    m = -n_x / n_y;                  // slope for y = mx + b
    b = (j + A[2]) - m * (i + A[0]); // y intercept for y = mx + b

    for (int k = -1; k < 2; k++)
    { // looping over 3x3 stencil centered about cell (i,j)
        for (int l = -1; l < 2; l++)
        {

            ip = i + l; // i,j for 3x3 stencil arround cell i,j
            jp = j + k;

            if (ip != i || jp != j)
            {
                y1 = m * (ip + 1) + b; // y value at right side of cell (ip,jp)
                y2 = m * (ip) + b;     // y value at left side of cell (ip,jp)

                if (y1 >= jp + 1.0 && y2 >= jp + 1.0)
                { // if y = mx + b is entirely above cell (ip,jp)
                    if (n_y > 0.0)
                    {
                        Cs = 0.0;
                    }
                    else
                    {
                        Cs = 1.0;
                    }
                }
                else if (y1 <= jp && y2 <= jp)
                { // if y = mx + b is entirely below cell (ip,jp)
                    if (n_y > 0.0)
                    {
                        Cs = 1.0;
                    }
                    else
                    {
                        Cs = 0.0;
                    }
                }
                else
                { // y = mx + b passes through cell (ip,jp)

                    if (n_y > 0.0)
                    { // if n_y pointing up, integrate vs upper bond of cell

                        Q = jp + 1.0;

                        for (int n = 0; n < nx_l + 1; n++)
                        {

                            x_l = ip + dx_l * n; // local discretization
                            y_l = m * x_l + b;

                            if (y_l > Q)
                            {
                                y_l = Q; // boundary condition if PLIC is above cell
                            }
                            else if (y_l < jp)
                            {
                                y_l = jp; // boundary condition if PLIC is below cell
                            }

                            dif[n] = Q - y_l;
                        }
                    }
                    else
                    { // if n_y pointing down, integrate vs lower bound of cell

                        Q = jp;

                        for (int n = 0; n < nx_l + 1; n++)
                        {

                            x_l = ip + dx_l * n; // local discretization
                            y_l = m * x_l + b;

                            if (y_l > jp + 1.0)
                            {
                                y_l = jp + 1.0; // boundary condition if PLIC is above cell
                            }
                            else if (y_l < jp)
                            {
                                y_l = jp; // boundary condition if PLIC is below cell
                            }

                            dif[n] = y_l - Q;
                        }
                    }

                    // integrating to estimate volume fraction over cell
                    Cs = (0.5) * dx_l * dif[0];

                    for (int n = 1; n < nx_l; n++)
                    {
                        Cs = Cs + dx_l * dif[n];
                    }

                    Cs = Cs + (0.5) * dx_l * dif[nx_l];
                }

                //  least squares error
                if (ecase == 1)
                {
                    E_l = pow(Cs - C_g[jp + 1][ip + 1], 2.0); // using C_g instead of C to handle boundary conditions
                }
                else
                {
                    E_l = fabs(Cs - C_g[jp + 1][ip + 1]); // using C_g instead of C to handle boundary conditions
                }

                err_el = err_el + E_l;
            }
        }
    }

    if (ecase == 1)
    {
        err_el = pow(err_el, 0.5);
    }

    return err_el;
}

double **Youngs_method(double **C, double **S, int m, int n, double dx0)
{
    // calculation of "normalized" normal vector components using Youngs method
    // "Normalized": |m_x| + |m_y| = 1
    // output columns: my_ip = [my_xp my_yp]

    int ij, ycase, bc;

    double *my_ip, my_x, my_y, my_xp, my_yp, dx0i, i16, i6, dx3i;
    double **myv_x, **myv_y, *A;

    ycase = 1; // 1 = Original Youngs Method	2 = Inward Scheme	3 = Outward scheme

    bc = 2; // 1 = no penetration       2 = periodic

    dx0i = 1.0 / dx0; // reciperol of values to prevent repeated division
    dx3i = 1.0 / (3.0 * dx0);
    i16 = 1.0 / 16.0;
    i6 = 1.0 / 6.0;

    my_ip = malloc(sizeof(double *) * 2);
    myv_x = malloc(sizeof(double *) * (n + 1));
    myv_y = malloc(sizeof(double *) * (n + 1));

    for (int j = 0; j < n + 1; j++)
    {
        myv_x[j] = malloc(sizeof(double *) * (m + 1));
        myv_y[j] = malloc(sizeof(double *) * (m + 1));
    }

    if (ycase == 1)
    { // Original

        double C_h1, C_h2, C_h3, C_v1, C_v2, C_v3;

        for (int j = 0; j < m + 1; j++)
        {
            for (int i = 0; i < n + 1; i++)
            {
                if (bc == 1) // No penetration
                {
                    if (i > 0 && i < n && j > 0 && j < m)
                    {

                        C_h1 = 0.5 * (C[j - 1][i - 1] + C[j][i - 1]); // Linear interpolating C to horizontal cell faces
                        C_h2 = 0.5 * (C[j - 1][i] + C[j][i]);

                        C_v1 = 0.5 * (C[j - 1][i - 1] + C[j - 1][i]); // Linear interpolating C to verticle faces
                        C_v2 = 0.5 * (C[j][i - 1] + C[j][i]);

                        myv_x[j][i] = -dx0i * (C_h1 - C_h2); // -dC/dx
                        myv_y[j][i] = -dx0i * (C_v1 - C_v2); // -dC/dy
                    }
                    else if (i == 0 && j > 0 && j < m)
                    {

                        C_h1 = 0.5 * (C[j - 1][0] + C[j][0]);
                        C_h2 = 0.5 * (C[j - 1][1] + C[j][1]);
                        C_h3 = 0.5 * (C[j - 1][2] + C[j][2]);

                        myv_x[j][i] = dx0i * (-2.0 * C_h1 + 3.0 * C_h2 - C_h3); // -dC/dx
                        myv_y[j][i] = 0.0;
                    }
                    else if (i == n && j > 0 && j < m)
                    {

                        C_h1 = 0.5 * (C[j - 1][n - 1] + C[j][n - 1]);
                        C_h2 = 0.5 * (C[j - 1][n - 2] + C[j][n - 2]);
                        C_h3 = 0.5 * (C[j - 1][n - 3] + C[j][n - 3]);

                        myv_x[j][i] = -dx0i * (-2.0 * C_h1 + 3.0 * C_h2 - C_h3); //-dC/dx
                        myv_y[j][i] = 0.0;
                    }
                    else if (i > 0 && i < n && j == 0)
                    {

                        C_v1 = 0.5 * (C[0][i - 1] + C[0][i]);
                        C_v2 = 0.5 * (C[1][i - 1] + C[1][i]);
                        C_v3 = 0.5 * (C[2][i - 1] + C[2][i]);

                        myv_x[j][i] = 0.0;
                        myv_y[j][i] = dx0i * (-2.0 * C_v1 + 3.0 * C_v2 - C_v3); // -dC/dy
                    }
                    else if (i > 0 && i < n && j == m)
                    {
                        C_v1 = 0.5 * (C[m - 1][i - 1] + C[m - 1][i]);
                        C_v2 = 0.5 * (C[m - 2][i - 1] + C[m - 2][i]);
                        C_v3 = 0.5 * (C[m - 3][i - 1] + C[m - 3][i]);

                        myv_x[j][i] = 0.0;
                        myv_y[j][i] = -dx0i * (-2.0 * C_v1 + 3.0 * C_v2 - C_v3); // -dC/dy
                    }
                    else if (i == 0 && j == 0)      // corner normal values
                    {
                        myv_x[j][i] = 0.0;
                        myv_y[j][i] = 0.0;
                    }
                    else if (i == n && j == 0)
                    {
                        myv_x[j][i] = 0.0;
                        myv_y[j][i] = 0.0;
                    }
                    else if (i == 0 && j == m)
                    {
                        myv_x[j][i] = 0.0;
                        myv_y[j][i] = 0.0;
                    }
                    else if (i == n && j == m)
                    {
                        myv_x[j][i] = 0.0;
                        myv_y[j][i] = 0.0;
                    }
                }
                else if (bc == 2)
                { // periodic

                    if (i > 0 && i < n && j > 0 && j < m)
                    {

                        C_h1 = 0.5 * (C[j - 1][i - 1] + C[j][i - 1]); // Linear interpolating C to horizontal cell faces
                        C_h2 = 0.5 * (C[j - 1][i] + C[j][i]);

                        C_v1 = 0.5 * (C[j - 1][i - 1] + C[j - 1][i]); // Linear interpolating C to verticle faces
                        C_v2 = 0.5 * (C[j][i - 1] + C[j][i]);

                        myv_x[j][i] = -dx0i * (C_h1 - C_h2); // -dC/dx
                        myv_y[j][i] = -dx0i * (C_v1 - C_v2); // -dC/dy
                    }
                    else if (i == 0 && j > 0 && j < m)
                    {

                        C_h1 = 0.5 * (C[j - 1][n - 1] + C[j][n - 1]);
                        C_h2 = 0.5 * (C[j - 1][0] + C[j][0]);

                        C_v1 = 0.5 * (C[j - 1][n - 1] + C[j - 1][0]);
                        C_v2 = 0.5 * (C[j][n - 1] + C[j][0]);

                        myv_x[j][i] = -dx0i * (C_h1 - C_h2); // -dC/dx
                        myv_y[j][i] = -dx0i * (C_v1 - C_v2); // -dC/dy
                    }
                    else if (i == n && j > 0 && j < m)
                    {
                        myv_x[j][i] = myv_x[j][0]; // periodic boundary conditons
                        myv_y[j][i] = myv_y[j][0];
                    }
                    else if (i > 0 && i < n && j == 0)
                    {

                        C_h1 = 0.5 * (C[0][i - 1] + C[m - 1][i - 1]);
                        C_h2 = 0.5 * (C[0][i] + C[m - 1][i]);

                        C_v1 = 0.5 * (C[m - 1][i - 1] + C[m - 1][i]);
                        C_v2 = 0.5 * (C[0][i - 1] + C[0][i]);

                        myv_x[j][i] = -dx0i * (C_h1 - C_h2); // -dC/dx
                        myv_y[j][i] = -dx0i * (C_v1 - C_v2); // -dC/dy
                    }
                    else if (i > 0 && i < n && j == m)
                    {
                        myv_x[j][i] = myv_x[0][i]; // periodic boundary conditions
                        myv_y[j][i] = myv_y[0][i];
                    }
                    else if (i == 0 && j == 0)
                    {
                        myv_x[j][i] = 0.0;
                        myv_y[j][i] = 0.0;
                    }
                    else if (i == n && j == 0)
                    {
                        myv_x[j][i] = 0.0;
                        myv_y[j][i] = 0.0;
                    }
                    else if (i == 0 && j == m)
                    {
                        myv_x[j][i] = 0.0;
                        myv_y[j][i] = 0.0;
                    }
                    else if (i == n && j == m)
                    {
                        myv_x[j][i] = 0.0;
                        myv_y[j][i] = 0.0;
                    }
                }
            }
        }
    }
    else if (ycase == 2)
    { // Inward Scheme
        printf("ycase = %i\n", ycase);
        double C1, C2, C3, C4, C_h1, C_h2, C_h3, C_h4, C_v1, C_v2, C_v3, C_v4;

        for (int j = 0; j < m + 1; j++)
        {
            for (int i = 0; i < n + 1; i++)
            {
                if (i == 0 && j == 0)
                { // bottom left corner.

                    C1 = 0.5 * (2.0 * C[m - 1][0] - C[m - 1][1] + 2.0 * C[0][n - 1] - C[1][n - 1]);
                    C2 = C[m - 1][0];
                    C3 = C[0][n - 1];
                    C4 = C[0][0];
                }
                else if (i == 0 && j == m)
                { // top left corner

                    C3 = 0.5 * (2.0 * C[m - 1][n - 1] - C[m - 2][n - 1] + 2.0 * C[0][0] - C[0][1]);
                    C1 = C[m - 1][n - 1];
                    C2 = C[m - 1][0];
                    C4 = C[0][0];
                }
                else if (i == n && j == m)
                { // top right corner

                    C4 = 0.5 * (2.0 * C[0][n - 1] - C[0][n - 2] + 2.0 * C[m - 1][0] - C[m - 2][0]);
                    C1 = C[m - 1][n - 1];
                    C2 = C[m - 1][0];
                    C3 = C[0][n - 1];
                }
                else if (i == n && j == 0)
                { // bottom right corner

                    C2 = 0.5 * (2.0 * C[m - 1][n - 1] - C[m - 1][n - 2] + 2.0 * C[0][0] - C[1][0]);
                    C1 = C[m - 1][n - 1];
                    C3 = C[0][n - 1];
                    C4 = C[0][0];
                }
                else if (i == 0)
                { // periodic bondary conditions

                    C1 = C[j - 1][n - 1];
                    C2 = C[j - 1][0];
                    C3 = C[j][n - 1];
                    C4 = C[j][0];
                }
                else if (j == 0)
                {

                    C1 = C[m - 1][i - 1];
                    C2 = C[m - 1][i];
                    C3 = C[0][i - 1];
                    C4 = C[0][i];
                }
                else if (i == m)
                {

                    C1 = C[j - 1][n - 1];
                    C2 = C[j - 1][0];
                    C3 = C[j][n - 1];
                    C4 = C[j][0];
                }
                else if (j == m)
                {

                    C1 = C[m - 1][i - 1];
                    C2 = C[m - 1][i];
                    C3 = C[0][i - 1];
                    C4 = C[0][i];
                }
                else
                { // cells needed for bilinear interpolation

                    C1 = C[j - 1][i - 1];
                    C2 = C[j - 1][i];
                    C3 = C[j][i - 1];
                    C4 = C[j][i];
                }

                C_v1 = 0.5 * (C1 + C2); // bilinear interpolation to verticle cell boundary
                C_v2 = 0.125 * (3.0 * C1 + 3.0 * C2 + C3 + C4);
                C_v3 = 0.125 * (C1 + C2 + 3.0 * C3 + 3.0 * C4);
                C_v4 = 0.5 * (C3 + C4);

                C_h1 = 0.5 * (C1 + C3); // bilinear interpolation to horizontal cell boundary
                C_h2 = 0.125 * (3.0 * C1 + C2 + 3.0 * C3 + C4);
                C_h3 = 0.125 * (C1 + 3.0 * C2 + C3 + 3.0 * C4);
                C_h4 = 0.5 * (C2 + C4);

                myv_x[j][i] = dx3i * (C_h1 - 8.0 * C_h2 + 8.0 * C_h3 - C_h4); // -dC/dx
                myv_y[j][i] = dx3i * (C_v1 - 8.0 * C_v2 + 8.0 * C_v3 - C_v4); // -dC/dy
            }
        }
    }
    else if (ycase == 3)
    { // Outward scheme

        double C1, C2, C3, C4, C_h1, C_h2, C_h3, C_h4, C_v1, C_v2, C_v3, C_v4;
        printf("ycase = %i\n", ycase);

        for (int j = 0; j < m + 1; j++)
        {
            for (int i = 0; i < n + 1; i++)
            {

                if (i > 1 && i < n - 1 && j > 1 && j < m - 1)
                {

                    C_h1 = 0.25 * (C[j - 1][i - 2] + C[j - 1][i - 1] + C[j][i - 2] + C[j][i - 1]);      // bilinear interpolation to horizontal face
                    C_h2 = 0.5 * (C[j - 1][i - 1] + C[j][i - 1]);
                    C_h3 = 0.5 * (C[j][i] + C[j - 1][i]);
                    C_h4 = 0.25 * (C[j - 1][i] + C[j - 1][i + 1] + C[j][i] + C[j][i + 1]);

                    C_v1 = 0.25 * (C[j - 2][i - 1] + C[j - 2][i] + C[j - 1][i - 1] + C[j - 1][i]);      // bilinear interpolation to verticel cell face
                    C_v2 = 0.5 * (C[j - 1][i - 1] + C[j - 1][i]);
                    C_v3 = 0.5 * (C[j][i - 1] + C[j][i]);
                    C_v4 = 0.25 * (C[j][i - 1] + C[j][i] + C[j + 1][i - 1] + C[j + 1][i]);
                }
                else if (i == 0 && j == 0)
                {

                    C_h1 = 0.25 * (C[m - 1][n - 2] + C[m - 1][n - 1] + C[0][n - 1] + C[0][n - 2]);
                    C_h2 = 0.5 * (C[0][n - 1] + C[m - 1][n - 1]);
                    C_h3 = 0.5 * (C[0][0] + C[m - 1][0]);
                    C_h4 = 0.25 * (C[m - 1][0] + C[m - 1][1] + C[0][0] + C[0][1]);

                    C_v1 = 0.25 * (C[m - 2][n - 1] + C[m - 2][0] + C[m - 1][n - 1] + C[m - 1][0]);
                    C_v2 = 0.5 * (C[m - 1][n - 1] + C[0][n - 1]);
                    C_v3 = 0.5 * (C[0][n - 1] + C[0][0]);
                    C_v4 = 0.25 * (C[0][n - 1] + C[0][0] + C[1][n - 1] + C[1][0]);
                }
                else if (i == 1 && j == 0)
                {

                    C_h1 = 0.25 * (C[m - 1][n - 1] + C[m - 1][0] + C[0][n - 1] + C[0][0]);
                    C_h2 = 0.5 * (C[m - 1][0] + C[0][0]);
                    C_h3 = 0.5 * (C[m - 1][1] + C[0][1]);
                    C_h4 = 0.25 * (C[m - 1][1] + C[m - 1][2] + C[0][1] + C[0][2]);

                    C_v1 = 0.25 * (C[m - 2][0] + C[m - 2][1] + C[m - 1][0] + C[m - 1][1]);
                    C_v2 = 0.5 * (C[m - 1][0] + C[m - 1][1]);
                    C_v3 = 0.5 * (C[0][0] + C[0][1]);
                    C_v4 = 0.25 * (C[0][0] + C[0][1] + C[1][0] + C[1][1]);
                }
                else if (i == 0 && j == 1)
                {

                    C_h1 = 0.25 * (C[0][n - 2] + C[0][n - 1] + C[1][n - 2] + C[1][n - 1]);
                    C_h2 = 0.5 * (C[0][n - 1] + C[1][n - 1]);
                    C_h3 = 0.5 * (C[0][0] + C[1][0]);
                    C_h4 = 0.25 * (C[0][0] + C[1][0] + C[0][1] + C[1][1]);

                    C_v1 = 0.25 * (C[m - 1][n - 1] + C[m - 1][0] + C[0][n - 1] + C[0][0]);
                    C_v2 = 0.5 * (C[0][m - 1] + C[0][0]);
                    C_v3 = 0.5 * (C[1][n - 1] + C[1][0]);
                    C_v4 = 0.25 * (C[1][n - 1] + C[1][0] + C[2][n - 1] + C[2][0]);
                }
                else if (i == 1 && j == 1)
                {

                    C_h1 = 0.25 * (C[0][n - 1] + C[0][0] + C[1][n - 1] + C[1][0]);
                    C_h2 = 0.5 * (C[0][0] + C[1][0]);
                    C_h3 = 0.5 * (C[0][1] + C[1][1]);
                    C_h4 = 0.25 * (C[0][1] + C[0][2] + C[1][1] + C[1][2]);

                    C_v1 = 0.25 * (C[m - 1][0] + C[m - 1][1] + C[0][0] + C[0][1]);
                    C_v2 = 0.5 * (C[0][0] + C[0][1]);
                    C_v3 = 0.5 * (C[1][0] + C[1][1]);
                    C_v4 = 0.25 * (C[1][0] + C[1][1] + C[2][0] + C[2][1]);
                }
                else if (i == n - 1 && j == 1)
                {

                    C_h1 = 0.25 * (C[0][n - 3] + C[0][n - 2] + C[1][n - 3] + C[1][n - 2]);
                    C_h2 = 0.5 * (C[0][n - 2] + C[1][n - 2]);
                    C_h3 = 0.5 * (C[0][n - 1] + C[1][n - 1]);
                    C_h4 = 0.25 * (C[0][n - 1] + C[0][0] + C[1][n - 1] + C[1][0]);

                    C_v3 = 0.25 * (C[m - 1][n - 2] + C[m - 1][n - 1] + C[0][n - 2] + C[0][n - 1]);
                    C_v2 = 0.5 * (C[0][n - 2] + C[0][n - 1]);
                    C_v3 = 0.5 * (C[1][n - 2] + C[1][n - 1]);
                    C_v4 = 0.25 * (C[1][n - 2] + C[1][n - 1] + C[2][n - 2] + C[2][n - 1]);
                }
                else if (i == n - 1 && j == 0)
                {

                    C_h1 = 0.25 * (C[m - 1][n - 3] + C[m - 1][n - 2] + C[0][n - 3] + C[0][n - 2]);
                    C_h2 = 0.5 * (C[m - 1][n - 2] + C[0][n - 2]);
                    C_h3 = 0.5 * (C[m - 1][n - 1] + C[0][n - 1]);
                    C_h4 = 0.25 * (C[m - 1][n - 1] + C[m - 1][0] + C[0][n - 1] + C[0][0]);

                    C_v1 = 0.25 * (C[m - 2][n - 2] + C[m - 2][n - 1] + C[m - 1][n - 2] + C[m - 1][n - 1]);
                    C_v2 = 0.5 * (C[m - 1][n - 2] + C[m - 1][n - 1]);
                    C_v3 = 0.5 * (C[0][n - 2] + C[0][n - 1]);
                    C_v4 = 0.25 * (C[0][n - 2] + C[0][n - 1] + C[1][n - 2] + C[1][n - 1]);
                }
                else if (i == n && j == 0)
                {

                    C_h1 = 0.25 * (C[m - 1][n - 2] + C[m - 1][n - 1] + C[0][n - 2] + C[0][n - 1]);
                    C_h2 = 0.5 * (C[m - 1][n - 1] + C[0][n - 1]);
                    C_h3 = 0.5 * (C[m - 1][0] + C[0][0]);
                    C_h4 = 0.25 * (C[m - 1][0] + C[m - 1][1] + C[0][0] + C[0][1]);

                    C_v1 = 0.25 * (C[m - 2][n - 1] + C[m - 2][0] + C[m - 1][n - 1] + C[m - 1][0]);
                    C_v2 = 0.5 * (C[m - 1][n - 1] + C[m - 1][0]);
                    C_v3 = 0.5 * (C[0][n - 1] + C[0][0]);
                    C_v4 = 0.25 * (C[0][n - 1] + C[0][0] + C[1][n - 1] + C[1][0]);
                }
                else if (i == n && j == 1)
                {

                    C_h1 = 0.25 * (C[0][n - 2] + C[0][n - 1] + C[1][n - 2] + C[1][n - 1]);
                    C_h2 = 0.5 * (C[0][n - 1] + C[1][n - 1]);
                    C_h3 = 0.5 * (C[0][0] + C[1][0]);
                    C_h4 = 0.25 * (C[0][0] + C[0][1] + C[1][0] + C[1][1]);

                    C_v1 = 0.25 * (C[m - 1][n - 1] + C[m - 1][0] + C[0][n - 1] + C[0][0]);
                    C_v2 = 0.5 * (C[0][n - 1] + C[0][0]);
                    C_v3 = 0.5 * (C[1][n - 1] + C[1][0]);
                    C_v4 = 0.25 * (C[1][n - 1] + C[1][0] + C[2][n - 1] + C[2][0]);
                }
                else if (i == 0 && j == m - 1)
                {

                    C_h1 = 0.25 * (C[m - 2][n - 2] + C[m - 2][n - 1] + C[m - 1][n - 2] + C[m - 1][n - 2]);
                    C_h2 = 0.5 * (C[m - 2][n - 1] + C[m - 2][0]);
                    C_h3 = 0.5 * (C[m - 2][0] + C[m - 1][0]);
                    C_h4 = 0.25 * (C[m - 2][0] + C[m - 2][1] + C[m - 1][0] + C[m - 1][1]);

                    C_v1 = 0.25 * (C[m - 3][n - 1] + C[m - 3][0] + C[m - 2][n - 1] + C[m - 2][0]);
                    C_v2 = 0.5 * (C[m - 2][n - 1] + C[m - 2][0]);
                    C_v3 = 0.5 * (C[m - 1][n - 1] + C[m - 1][0]);
                    C_v4 = 0.25 * (C[m - 1][n - 1] + C[m - 1][0] + C[0][n - 1] + C[0][0]);
                }
                else if (i == 1 && j == m - 1)
                {

                    C_h1 = 0.25 * (C[m - 2][n - 1] + C[m - 2][0] + C[m - 1][n - 1] + C[m - 1][0]);
                    C_h2 = 0.5 * (C[m - 2][0] + C[m - 1][0]);
                    C_h3 = 0.5 * (C[m - 2][1] + C[m - 1][1]);
                    C_h4 = 0.25 * (C[m - 2][1] + C[m - 2][2] + C[m - 1][1] + C[m - 1][2]);

                    C_v1 = 0.25 * (C[m - 3][0] + C[m - 3][1] + C[m - 2][0] + C[m - 2][1]);
                    C_v2 = 0.5 * (C[m - 2][0] + C[m - 2][1]);
                    C_v3 = 0.5 * (C[m - 1][0] + C[m - 1][1]);
                    C_v4 = 0.25 * (C[m - 1][0] + C[m - 1][1] + C[0][0] + C[0][1]);
                }
                else if (i == 0 && j == m)
                {

                    C_h1 = 0.25 * (C[m - 1][n - 2] + C[m - 1][n - 1] + C[0][n - 2] + C[0][n - 1]);
                    C_h2 = 0.5 * (C[m - 1][n - 1] + C[0][n - 1]);
                    C_h3 = 0.5 * (C[m - 1][0] + C[0][0]);
                    C_h4 = 0.25 * (C[m - 1][0] + C[m - 1][1] + C[0][0] + C[0][1]);

                    C_v1 = 0.25 * (C[m - 2][n - 1] + C[m - 2][0] + C[m - 1][n - 1] + C[m - 1][0]);
                    C_v2 = 0.5 * (C[m - 1][n - 1] + C[m - 1][0]);
                    C_v3 = 0.5 * (C[0][n - 1] + C[0][0]);
                    C_v4 = 0.25 * (C[0][n - 1] + C[0][0] + C[1][n - 1] + C[1][0]);
                }
                else if (i == 1 && j == m)
                {

                    C_h1 = 0.25 * (C[m - 1][n - 1] + C[m - 1][0] + C[0][n - 1] + C[0][0]);
                    C_h2 = 0.5 * (C[m - 1][0] + C[0][0]);
                    C_h3 = 0.5 * (C[m - 1][1] + C[0][1]);
                    C_h4 = 0.25 * (C[m - 1][1] + C[m - 1][2] + C[0][1] + C[0][2]);

                    C_v1 = 0.25 * (C[m - 2][0] + C[m - 2][1] + C[m - 1][0] + C[m - 1][1]);
                    C_v2 = 0.5 * (C[m - 1][0] + C[m - 1][1]);
                    C_v3 = 0.5 * (C[0][0] + C[0][1]);
                    C_v4 = 0.25 * (C[0][0] + C[0][1] + C[1][0] + C[1][1]);
                }
                else if (i == n && j == m - 1)
                {

                    C_h1 = 0.25 * (C[m - 2][n - 2] + C[m - 2][n - 1] + C[m - 1][n - 2] + C[m - 1][n - 1]);
                    C_h2 = 0.5 * (C[m - 2][n - 1] + C[m - 1][n - 1]);
                    C_h3 = 0.5 * (C[m - 2][0] + C[m - 1][0]);
                    C_h4 = 0.25 * (C[m - 2][0] + C[m - 2][1] + C[m - 1][0] + C[m - 1][1]);

                    C_v1 = 0.25 * (C[m - 3][n - 1] + C[m - 3][0] + C[m - 2][n - 1] + C[m - 2][0]);
                    C_v2 = 0.5 * (C[m - 2][n - 1] + C[m - 2][0]);
                    C_v3 = 0.5 * (C[m - 1][n - 1] + C[m - 1][0]);
                    C_v4 = 0.25 * (C[m - 1][n - 1] + C[m - 1][0] + C[0][n - 1] + C[0][0]);
                }
                else if (i == n - 1 && j == m - 1)
                {

                    C_h1 = 0.25 * (C[m - 2][n - 3] + C[m - 2][n - 2] + C[m - 1][n - 3] + C[m - 2][n - 2]);
                    C_h2 = 0.5 * (C[m - 2][n - 2] + C[m - 1][n - 2]);
                    C_h3 = 0.5 * (C[m - 2][n - 1] + C[m - 1][n - 1]);
                    C_h4 = 0.25 * (C[m - 2][n - 1] + C[m - 2][0] + C[m - 1][n - 1] + C[m - 1][0]);

                    C_v1 = 0.25 * (C[m - 3][n - 2] + C[m - 3][n - 1] + C[m - 2][n - 2] + C[m - 2][n - 1]);
                    C_v2 = 0.5 * (C[m - 2][n - 2] + C[m - 2][n - 1]);
                    C_v3 = 0.5 * (C[m - 1][n - 2] + C[m - 1][n - 1]);
                    C_v4 = 0.25 * (C[m - 1][n - 2] + C[m - 1][n - 1] + C[0][n - 2] + C[0][n - 1]);
                }
                else if (i == n - 1 && j == m)
                {

                    C_h1 = 0.25 * (C[m - 1][n - 3] + C[m - 1][n - 2] + C[0][n - 3] + C[0][n - 2]);
                    C_h2 = 0.5 * (C[m - 1][n - 2] + C[0][n - 2]);
                    C_h3 = 0.5 * (C[m - 1][n - 1] + C[0][n - 1]);
                    C_h4 = 0.5 * (C[m - 1][n - 1] + C[m - 1][0] + C[0][n - 1] + C[0][0]);

                    C_v1 = 0.25 * (C[m - 2][n - 2] + C[m - 2][n - 1] + C[m - 1][n - 2] + C[m - 1][n - 1]);
                    C_v2 = 0.5 * (C[m - 1][n - 2] + C[m - 1][n - 1]);
                    C_v3 = 0.5 * (C[0][n - 2] + C[0][n - 1]);
                    C_v4 = 0.25 * (C[0][n - 2] + C[0][n - 1] + C[1][n - 2] + C[1][n - 1]);
                }
                else if (i == n && j == m)
                {

                    C_h1 = 0.25 * (C[m - 1][n - 2] + C[m - 1][n - 1] + C[0][n - 2] + C[0][n - 1]);
                    C_h2 = 0.5 * (C[m - 1][n - 1] + C[0][n - 1]);
                    C_h3 = 0.5 * (C[m - 1][0] + C[0][0]);
                    C_h4 = 0.25 * (C[m - 1][0] + C[m - 1][1] + C[0][0] + C[0][1]);

                    C_v1 = 0.25 * (C[m - 2][n - 1] + C[m - 2][0] + C[m - 1][n - 1] + C[m - 1][0]);
                    C_v2 = 0.5 * (C[m - 1][n - 1] + C[m - 1][0]);
                    C_v3 = 0.25 * (C[0][n - 1] + C[0][0]);
                    C_v4 = 0.5 * (C[0][n - 1] + C[0][0] + C[1][n - 1] + C[1][0]);
                }
                else if (i == 0)
                {

                    C_h1 = 0.25 * (C[j - 1][n - 2] + C[j - 1][n - 1] + C[j][n - 2] + C[j][n - 1]);
                    C_h2 = 0.5 * (C[j - 1][n - 1] + C[j][n - 1]);
                    C_h3 = 0.5 * (C[j - 1][0] + C[j][0]);
                    C_h4 = 0.25 * (C[j - 1][0] + C[j - 1][1] + C[j][0] + C[j][1]);

                    C_v1 = 0.25 * (C[j - 2][n - 1] + C[j - 2][0] + C[j - 1][n - 1] + C[j - 1][0]);
                    C_v2 = 0.5 * (C[j - 1][n - 1] + C[j - 1][0]);
                    C_v3 = 0.5 * (C[j][n - 1] + C[j][0]);
                    C_v4 = 0.25 * (C[j][n - 1] + C[j][0] + C[j + 1][n - 1] + C[j + 1][0]);
                }
                else if (i == 1)
                {

                    C_h1 = 0.25 * (C[j - 1][n - 1] + C[j - 1][0] + C[j][n - 1] + C[j][0]);
                    C_h2 = 0.5 * (C[j - 1][0] + C[j][0]);
                    C_h3 = 0.5 * (C[j - 1][1] + C[j][1]);
                    C_h4 = 0.25 * (C[j - 1][i] + C[j - 1][i + 1] + C[j][i] + C[j][i + 1]);

                    C_v1 = 0.25 * (C[j - 2][i - 1] + C[j - 1][i] + C[j - 1][i - 1] + C[j - 1][i]);
                    C_v2 = 0.5 * (C[j - 1][i - 1] + C[j - 1][i]);
                    C_v3 = 0.5 * (C[j][i - 1] + C[j][i]);
                    C_v4 = 0.25 * (C[j][i - 1] + C[j][i] + C[j + 1][i - 1] + C[j + 1][i]);
                }
                else if (i == m - 1)
                {

                    C_h1 = 0.25 * (C[j - 1][i - 2] + C[j - 1][i - 1] + C[j][i - 2] + C[j][i - 1]);
                    C_h2 = 0.5 * (C[j - 1][i - 1] + C[j][i - 1]);
                    C_h3 = 0.5 * (C[j - 1][i] + C[j][i]);
                    C_h4 = 0.25 * (C[j - 1][i] + C[j - 1][0] + C[j][0]);

                    C_v1 = 0.25 * (C[j - 2][i - 1] + C[j - 1][i] + C[j - 1][i - 1] + C[j - 1][i]);
                    C_v2 = 0.5 * (C[j - 1][i - 1] + C[j - 1][i]);
                    C_v3 = 0.5 * (C[j][i - 1] + C[j][i]);
                    C_v4 = 0.25 * (C[j][i - 1] + C[j][i] + C[j + 1][i - 1] + C[j + 1][i]);
                }
                else if (i == m)
                {

                    C_h1 = 0.25 * (C[j - 1][n - 2] + C[j - 1][n - 1] + C[j][n - 2] + C[j][n - 1]);
                    C_h2 = 0.5 * (C[j - 1][i - 1] + C[j][i - 1]);
                    C_h3 = 0.5 * (C[j - 1][0] + C[j][0]);
                    C_h4 = 0.5 * (C[j - 1][0] + C[j - 1][1] + C[j][0] + C[j][1]);

                    C_v1 = 0.25 * (C[j - 2][n - 1] + C[j - 2][0] + C[j - 1][n - 1] + C[j - 1][0]);
                    C_v2 = 0.5 * (C[j - 1][n - 1] + C[j - 1][0]);
                    C_v3 = 0.5 * (C[j][n - 1] + C[j][0]);
                    C_v4 = 0.25 * (C[j][n - 1] + C[j][0] + C[j + 1][n - 1] + C[j + 1][0]);
                }
                else if (j == 0)
                {

                    C_h1 = 0.25 * (C[m - 1][i - 2] + C[m - 1][i - 1] + C[0][i - 2] + C[0][i - 1]);
                    C_h2 = 0.5 * (C[m - 1][i - 1] + C[0][i - 1]);
                    C_h3 = 0.5 * (C[m - 1][i] + C[0][i]);
                    C_h4 = 0.25 * (C[m - 1][i] + C[m - 1][i + 1] + C[0][i] + C[0][i + 1]);

                    C_v1 = 0.25 * (C[m - 2][i - 1] + C[m - 2][i] + C[m - 1][i - 1] + C[m - 1][i]);
                    C_v2 = 0.5 * (C[m - 1][i - 1] + C[m - 1][i]);
                    C_v3 = 0.5 * (C[0][i - 1] + C[0][i]);
                    C_v4 = 0.25 * (C[0][i - 1] + C[0][i] + C[1][i - 1] + C[1][i]);
                }
                else if (j == 1)
                {

                    C_h1 = 0.25 * (C[j - 1][i - 2] + C[j - 1][i - 1] + C[j][i - 2] + C[j][i - 1]);
                    C_h2 = 0.5 * (C[j - 1][i - 1] + C[j][i - 1]);
                    C_h3 = 0.5 * (C[j - 1][i] + C[j][i]);
                    C_h4 = 0.25 * (C[j - 1][i] + C[j - 1][i + 1] + C[j][i] + C[j][i + 1]);

                    C_v1 = 0.25 * (C[m - 1][i - 1] + C[m - 1][i] + C[j - 1][i - 1] + C[j - 1][i]);
                    C_v2 = 0.5 * (C[j - 1][i - 1] + C[j - 1][i]);
                    C_v3 = 0.5 * (C[j][i - 1] + C[j][i]);
                    C_v4 = 0.5 * (C[j][i - 1] + C[j][i] + C[j + 1][i - 1] + C[j + 1][i]);
                }
                else if (j == n - 1)
                {

                    C_h1 = 0.25 * (C[j - 1][i - 2] + C[j - 1][i - 1] + C[j][i - 2] + C[j][i - 1]);
                    C_h2 = 0.5 * (C[j - 1][i - 1] + C[j][i - 1]);
                    C_h3 = 0.5 * (C[j - 1][i] + C[j][i]);
                    C_h4 = 0.25 * (C[j - 1][i] + C[j - 1][i + 1] + C[j][i] + C[j][i + 1]);

                    C_v1 = 0.25 * (C[j - 2][i - 1] + C[j - 2][i] + C[j - 1][i - 1] + C[j - 1][i]);
                    C_v2 = 0.5 * (C[j - 1][i - 1] + C[j - 1][i]);
                    C_v3 = 0.5 * (C[j][i - 1] + C[j][i]);
                    C_v4 = 0.25 * (C[j][i - 1] + C[j][i] + C[0][i - 1] + C[0][i]);
                }
                else if (j == n)
                {

                    C_h1 = 0.25 * (C[j - 1][i - 2] + C[j - 1][i - 1] + C[0][i - 2] + C[0][i - 2]);
                    C_h2 = 0.5 * (C[j - 1][i - 1] + C[0][i - 1]);
                    C_h3 = 0.5 * (C[j - 1][i] + C[0][i]);
                    C_h4 = 0.25 * (C[j - 1][i] + C[j - 1][i + 1] + C[0][i] + C[0][i + 1]);

                    C_v1 = 0.25 * (C[j - 2][i - 1] + C[j - 2][i] + C[j - 1][i - 1] + C[j - 1][i]);
                    C_v2 = 0.5 * (C[j - 1][i - 1] + C[j - 1][i]);
                    C_v3 = 0.5 * (C[0][i - 1] + C[0][i]);
                    C_v4 = 0.25 * (C[0][i - 1] + C[0][i] + C[1][i - 1] + C[1][i]);
                }

                myv_x[j][i] = dx0i * i6 * (C_h1 - 8.0 * C_h2 + 8.0 * C_h3 - C_h4);      // myv_x = dC/dx
                myv_y[j][i] = dx0i * i6 * (C_v1 - 8.0 * C_v2 + 8.0 * C_v3 - C_v4);      // myv_y = dC/dy
            }
        }
    }

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            if (C[j][i] > 1e-10 && C[j][i] < 1.0 - 1e-10)
            {

                my_x = (myv_x[j][i] + myv_x[j + 1][i] + myv_x[j][i + 1] + myv_x[j + 1][i + 1]) * 0.25; // average slope at cell centers
                my_y = (myv_y[j][i] + myv_y[j + 1][i] + myv_y[j][i + 1] + myv_y[j + 1][i + 1]) * 0.25;

                my_xp = my_x / (fabs(my_x) + fabs(my_y)); // "normalizing" my_x and my_y
                my_yp = my_y / (fabs(my_x) + fabs(my_y)); // so ||my_x|| + ||my_y|| = 1

                if (fabs(my_yp) < 1e-10)
                { // ensuring slope never approaches infinity
                    if (my_yp > 0.0)
                    {
                        my_yp = 1e-10;
                    }
                    else
                    {
                        my_yp = -1e-10;
                    }
                }

                ij = i + j * m; // linearizing matrix

                my_ip[0] = my_xp;
                my_ip[1] = my_yp;

                A = VoF_tools(C[j][i], my_ip, dx0, m, n); // VoF tools

                S[0][ij] = A[0]; // x1
                S[1][ij] = A[1]; // x2
                S[2][ij] = A[2]; // y1
                S[3][ij] = A[3]; // y2
                S[4][ij] = A[4]; // n_x
                S[5][ij] = A[5]; // n_y
                S[6][ij] = 2.0;  // indicating youngs method

                free(A);
            }
        }
    }

    for (int i = 0; i < m + 1; i++)
    {
        free(myv_x[i]);
        free(myv_y[i]);
    }
    free(my_ip);
    free(myv_x);
    free(myv_y);

    return S;
}

double **CC_method(double **C, double **S, int m, int n, double dx0)
{
    // calculation of "normalized" normal vector components using centered columns method
    // "Normalized": |m_x| + |m_y| = 1
    // output rows mci_jp = [mcy_xp; mcy_yp; mcx_xp; mcx_yp]

    int ij, ccase, bc;
    double dx0i, dx0i2, C_g[m + 2][n + 2], hgt_x[m][n], hgt_y[m][n], mcy_y, mcy_x, sgn_mcy_y, sgn_mcx_x;
    double mcx_x, mcx_y, mcx_xp, mcx_yp, mcy_yp, mcy_xp, *m_i, *A1, **mci_jp, dir;

    m_i = malloc(sizeof(double *) * 2);

    mci_jp = malloc(sizeof(double *) * 4);
    for (int j = 0; j < 4; j++)
    {
        mci_jp[j] = malloc(sizeof(double *) * m * n);
    }

    ccase = 1; // 1 = 3x3 stencil		2 = extendted stencil

    bc = 2; // 1 = no penetration       2 = periodic

    dx0i = 1.0 / dx0; // inverting dx to prevent repeated division (for computational speed)
    dx0i2 = 1.0 / (2.0 * dx0);

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            C_g[0][i + 1] = C[m - 1][i]; // adding periodic ghost values to function
            C_g[m + 1][i + 1] = C[0][i];
            C_g[j + 1][0] = C[j][n - 1];
            C_g[j + 1][n + 1] = C[j][0];
            C_g[j + 1][i + 1] = C[j][i];
        }
    }

    C_g[0][0] = 0.5 * (C_g[0][1] + C_g[1][0]); // Corner Ghost values
    C_g[0][n + 1] = 0.5 * (C_g[1][n + 1] + C_g[0][n]);
    C_g[m + 1][n + 1] = 0.5 * (C_g[m][n + 1] + C_g[m + 1][n]);
    C_g[m + 1][0] = 0.5 * (C_g[m][0] + C_g[m + 1][1]);

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            hgt_x[j][i] = 0.0; // initialzing arrays
            hgt_y[j][i] = 0.0;
        }
    }

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            for (int k = -1; k < 2; k++)
            {
                hgt_x[j][i] = hgt_x[j][i] + dx0 * C_g[j + 1][i + 1 + k]; // calculating height functions
                hgt_y[j][i] = hgt_y[j][i] + dx0 * C_g[j + 1 + k][i + 1];
            }
        }
    }

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            if (C[j][i] > 1e-10 && C[j][i] < 1.0 - 1e-10)
            {

                if (i == 0)
                { // periodic boundary conditions

                    mcy_x = dx0i2 * (hgt_y[j][i + 1] - hgt_y[j][n - 1]);
                    mcx_x = dx0i2 * (C[j][n - 1] - C[j][i + 1]);
                }
                else if (i == n - 1)
                {

                    mcy_x = dx0i2 * (hgt_y[j][0] - hgt_y[j][i - 1]);
                    mcx_x = dx0i2 * (C[j][i - 1] - C[j][0]);
                }
                else
                {
                    mcy_x = dx0i2 * (hgt_y[j][i + 1] - hgt_y[j][i - 1]); // mcy_x = d(hgt_y)/dx
                    mcx_x = dx0i2 * (C[j][i - 1] - C[j][i + 1]);         // mcx_x = -dC/dx
                }

                if (j == 0)
                { // periodic boundary conditions

                    mcy_y = dx0i2 * (C[m - 1][i] - C[j + 1][i]);
                    mcx_y = dx0i2 * (hgt_x[j + 1][i] - hgt_x[m - 1][i]);
                }
                else if (j == m - 1)
                {

                    mcy_y = dx0i2 * (C[j - 1][i] - C[0][i]);
                    mcx_y = dx0i2 * (hgt_x[0][i] - hgt_x[j - 1][i]);
                }
                else
                {
                    mcy_y = dx0i2 * (C[j - 1][i] - C[j + 1][i]);         // mcy_y = -dC/dy
                    mcx_y = dx0i2 * (hgt_x[j + 1][i] - hgt_x[j - 1][i]); // mcx_y = d(hgt_x)/dy
                }

                if (mcy_y < 0.0)
                {
                    sgn_mcy_y = 1.0; // calculating the sign of mcy_y
                }
                else
                {
                    sgn_mcy_y = -1.0;
                }

                mcy_yp = sgn_mcy_y / (1.0 + fabs(mcy_x)); // "normalizing" mcy_y and mcy_x so
                mcy_xp = mcy_x / (1.0 + fabs(mcy_x));     // ||mcy_y|| + ||mcy_x|| = 1

                if (fabs(mcy_yp) < 1e-10)
                { // ensuring slope never approaches infinity
                    if (mcy_yp > 0.0)
                    {
                        mcy_yp = 1e-10;
                    }
                    else
                    {
                        mcy_yp = -1e-10;
                    }
                }

                if (mcx_x < 0.0)
                {
                    sgn_mcx_x = 1.0; // calculating sgn of mcx_x
                }
                else
                {
                    sgn_mcx_x = -1.0;
                }

                mcx_xp = sgn_mcx_x / (1.0 + fabs(mcx_y)); // "normalizing" mcx_x and mcx_y so
                mcx_yp = mcx_y / (1.0 + fabs(mcx_y));     // ||mcx_x|| + ||mcx_y|| = 1

                if (fabs(mcx_yp) < 1e-10)
                { // ensuring slope never approaches infinity
                    if (mcx_yp > 0.0)
                    {
                        mcx_yp = 1e-10;
                    }
                    else
                    {
                        mcx_yp = -1e-10;
                    }
                }

                ij = i + j * m; // linearizing matrix

                mci_jp[0][ij] = mcy_xp;         // interface normal compoents for both
                mci_jp[1][ij] = mcy_yp;         // interface represenations
                mci_jp[2][ij] = mcx_xp;
                mci_jp[3][ij] = mcx_yp;
            }
        }
    }
    int cnt = 0;
    int cnt2 = 0;

    if (ccase == 1)
    { // normal Centered columns method
        for (int j = 0; j < m; j++)
        {
            for (int i = 0; i < n; i++)
            {
                if (C[j][i] > 1e-10 && C[j][i] < 1.0 - 1e-10)
                {

                    ij = i + j * m;

                    if (fabs(mci_jp[1][ij]) > fabs(mci_jp[2][ij]))
                    { // selecting between interface norms
                        m_i[0] = mci_jp[0][ij];
                        m_i[1] = mci_jp[1][ij];
                        dir = 1.0; // y representaiton
                    }
                    else
                    {
                        m_i[0] = mci_jp[2][ij];
                        m_i[1] = mci_jp[3][ij];
                        dir = 2.0; // x representation
                    }

                    A1 = VoF_tools(C[j][i], m_i, dx0, m, n); // VoF tools

                    S[0][ij] = A1[0]; // x1
                    S[1][ij] = A1[1]; // x2
                    S[2][ij] = A1[2]; // y1
                    S[3][ij] = A1[3]; // y2
                    S[4][ij] = A1[4]; // n_x
                    S[5][ij] = A1[5]; // n_y
                    S[6][ij] = dir;   // interface represenations

                    free(A1);
                }
            }
        }
    }
    else if (ccase == 2)
    { // adding cells to height function if needed

        int w;
        double *A2, m0, b0, y3, hm1, hp1;

        for (int j = 0; j < m; j++)
        {
            for (int i = 0; i < n; i++)
            {
                if (C[j][i] > 1e-10 && C[j][i] < 1.0 - 1e-10)
                {

                    cnt2++;
                    ij = i + j * m;

                    if (i > 0 && i < n - 1 && j > 0 && j < m - 1)
                    {
                        m_i[0] = mci_jp[0][ij]; // mcy_xp
                        m_i[1] = mci_jp[1][ij]; // mcy_yp

                        A1 = VoF_tools(C[j][i], m_i, dx0, m, n); // VoF tools

                        m0 = -A1[4] / A1[5];                     // parameters for y = mx + b form
                        b0 = (1.0 + A1[2]) - m0 * (1.0 + A1[0]); // b0 also is y value at left side of 3x3 stencil

                        y3 = m0 * 3.0 + b0; // y value at right side of 3x3 stencil

                        if (b0 <= 0.0 || b0 >= 3.0 || y3 <= 0.0 || y3 >= 3.0)
                        { // intersecting adjecent sides of 3x3 stencil

                            cnt++;

                            w = 1; // indicating Vof tools might need to be run again for this interface representation

                            hm1 = dx0 * (C_g[j - 1][i] + C_g[j][i] + C_g[j + 1][i] + C_g[j + 2][i] + C_g[j + 3][i]);                     // height function using 5x3 stencil for cell (i-1,j)
                            hp1 = dx0 * (C_g[j - 1][i + 2] + C_g[j][i + 2] + C_g[j + 1][i + 2] + C_g[j + 2][i + 2] + C_g[j + 3][i + 2]); // height function using 5x3 stencil for cell (i+1,j)

                            mcy_x = dx0i2 * (hp1 - hm1);                 // mcy_x = d(hgt_y)/dx
                            mcy_y = dx0i2 * (C[j - 1][i] - C[j + 1][i]); // mcy_y = -dC/dy

                            if (mcy_y < 0.0)
                            {
                                sgn_mcy_y = 1.0; // calculating the sign of mcy_y
                            }
                            else
                            {
                                sgn_mcy_y = -1.0;
                            }

                            mcy_yp = sgn_mcy_y / (1.0 + fabs(mcy_x)); // "normalizing" mcy_y and mcy_x so
                            mcy_xp = mcy_x / (1.0 + fabs(mcy_x));     // ||mcy_y|| + ||mcy_x|| = 1

                            if (fabs(mcy_yp) < 1e-10)
                            { // ensuring slope never approaches infinity
                                if (mcy_yp > 0.0)
                                {
                                    mcy_yp = 1e-10;
                                }
                                else
                                {
                                    mcy_yp = -1e-10;
                                }
                            }

                            mci_jp[0][ij] = mcy_xp;
                            mci_jp[1][ij] = mcy_yp;
                        }
                        else
                        {
                            w = 0;
                        }

                        m_i[0] = mci_jp[2][ij]; // mcx_xp
                        m_i[1] = mci_jp[3][ij]; // mcx_yp

                        A2 = VoF_tools(C[j][i], m_i, dx0, m, n); // VoF tools

                        m0 = -A2[4] / A2[5];                     // parameters for y = mx + b form
                        b0 = (1.0 + A2[2]) - m0 * (1.0 + A2[0]); // b0 also is y value at left side of 3x3 stencil

                        y3 = m0 * 3.0 + b0; // y value at right side of 3x3 stencil

                        if (b0 <= 0.0 || b0 >= 3.0 || y3 <= 0.0 || y3 >= 3.0)
                        { // intersecting adjecent sides of 3x3 stencil

                            cnt++;

                            w = 1; // indicating VOF tools might need to be rerun for this interface representation

                            hm1 = dx0 * (C_g[j][i - 1] + C_g[j][i] + C_g[j][i + 1] + C_g[j][i + 2] + C_g[j][i + 3]);
                            hp1 = dx0 * (C_g[j + 2][i - 1] + C_g[j + 2][i] + C_g[j + 2][i + 1] + C_g[j + 2][i + 2] + C_g[j + 2][i + 3]);

                            mcx_y = dx0i2 * (hp1 - hm1);                 // mcx_y = d(hgt_x)/dx
                            mcx_x = dx0i2 * (C[j][i - 1] - C[j][i + 1]); // mcx_x = -dC/dx

                            if (mcx_x < 0.0)
                            {
                                sgn_mcx_x = 1.0; // calculating sgn of mcx_x
                            }
                            else
                            {
                                sgn_mcx_x = -1.0;
                            }

                            mcx_xp = sgn_mcx_x / (1.0 + fabs(mcx_y)); // "normalizing" mcx_x and mcx_y so
                            mcx_yp = mcx_y / (1.0 + fabs(mcx_y));     // ||mcx_x|| + ||mcx_y|| = 1

                            if (fabs(mcx_yp) < 1e-10)
                            { // ensuring slope never approaches infinity
                                if (mcx_yp > 0.0)
                                {
                                    mcx_yp = 1e-10;
                                }
                                else
                                {
                                    mcx_yp = -1e-10;
                                }
                            }

                            mci_jp[2][ij] = mcx_xp;
                            mci_jp[3][ij] = mcx_yp;
                        }
                        else
                        {
                            w = 0;
                        }

                        if (w == 0 || i == 1 || i == n - 2 || j == 1 || j == m - 2)
                        { // no extra cells added; Vof tools already run for these configuration
                            if (fabs(mci_jp[1][ij]) > fabs(mci_jp[2][ij]))
                            { // selecting between interface norms

                                S[0][ij] = A1[0]; // x1
                                S[1][ij] = A1[1]; // x2
                                S[2][ij] = A1[2]; // y1
                                S[3][ij] = A1[3]; // y2
                                S[4][ij] = A1[4]; // n_x
                                S[5][ij] = A1[5]; // n_y
                                S[6][ij] = 1.0;   // y representation

                                free(A1);
                            }
                            else
                            {

                                S[0][ij] = A2[0]; // x1
                                S[1][ij] = A2[1]; // x2
                                S[2][ij] = A2[2]; // y1
                                S[3][ij] = A2[3]; // y2
                                S[4][ij] = A2[4]; // n_x
                                S[5][ij] = A2[5]; // n_y
                                S[6][ij] = 2.0;   // x representation

                                free(A2);
                            }
                        }
                        else if (w == 1 && i > 1 && i < n - 2 && j > 1 && j < m - 2)
                        {
                            if (fabs(mci_jp[1][ij]) > fabs(mci_jp[2][ij]))
                            { // selecting between interface norms
                                m_i[0] = mci_jp[0][ij];
                                m_i[1] = mci_jp[1][ij];
                                dir = 1.0;
                            }
                            else
                            {
                                m_i[0] = mci_jp[2][ij];
                                m_i[1] = mci_jp[3][ij];
                                dir = 2.0;
                            }

                            A1 = VoF_tools(C[j][i], m_i, dx0, m, n); // VoF tools

                            S[0][ij] = A1[0]; // x1
                            S[1][ij] = A1[1]; // x2
                            S[2][ij] = A1[2]; // y1
                            S[3][ij] = A1[3]; // y2
                            S[4][ij] = A1[4]; // n_x
                            S[5][ij] = A1[5]; // n_y
                            S[6][ij] = dir;

                            free(A1);
                        }
                    }
                    else
                    {
                        if (bc == 1)
                        { // no penetration
                            if (i == 0 && j > 0 && j < m - 1)
                            { // no penetration boundary condition
                                m_i[0] = 1.0;
                                m_i[1] = 1e-10;
                                dir = 0.0;
                            }
                            else if (i == n - 1 && j > 0 && j < m - 1)
                            {
                                m_i[0] = -1.0;
                                m_i[1] = -1e-10;
                                dir = 0.0;
                            }
                            else if (i > 0 && i < n - 1 && j == 0)
                            {
                                m_i[0] = 0.0;
                                m_i[1] = 1.0;
                                dir = 0.0;
                            }
                            else if (i > 0 && i < n - 1 && j == m - 1)
                            {
                                m_i[0] = 0.0;
                                m_i[1] = -1.0;
                                dir = 0.0;
                            }
                            else if (i == 0 && j == 0)
                            {
                                m_i[0] = 1.0;
                                m_i[1] = 1.0;
                            }
                            else if (i == n - 1 && j == 0)
                            {
                                m_i[0] = -1.0;
                                m_i[0] = 1.0;
                            }
                            else if (i == 0 && j == m - 1)
                            {
                                m_i[0] = 1.0;
                                m_i[1] = -1.0;
                            }
                            else if (i == n - 1 && j == n - 1)
                            {
                                m_i[0] = -1.0;
                                m_i[1] = -1.0;
                            }
                        }
                        else
                        { // periodic
                            if (fabs(mci_jp[1][ij]) > fabs(mci_jp[2][ij]))
                            { // selecting between interface norms
                                m_i[0] = mci_jp[0][ij];
                                m_i[1] = mci_jp[1][ij];
                                dir = 1.0;
                            }
                            else
                            {
                                m_i[0] = mci_jp[2][ij];
                                m_i[1] = mci_jp[3][ij];
                                dir = 2.0;
                            }
                        }

                        A1 = VoF_tools(C[j][i], m_i, dx0, m, n); // VoF tools

                        S[0][ij] = A1[0]; // x1
                        S[1][ij] = A1[1]; // x2
                        S[2][ij] = A1[2]; // y1
                        S[3][ij] = A1[3]; // y2
                        S[4][ij] = A1[4]; // n_x
                        S[5][ij] = A1[5]; // n_y
                        S[6][ij] = dir;

                        free(A1);
                    }
                }
            }
        }
    }

    for (int i = 0; i < 4; i++)
    {
        free(mci_jp[i]);
    }

    free(mci_jp);
    free(m_i);

    return S;
}

double **myc_method(double **C, double **S, int m, int n, double dx0)
{
    // calculation of "normalized" normal vector components using Mixed Youngs Centered Method
    // "Normalized": |m_x| + |m_y| = 1
    // output rows m_i = [m_x; m_y]

    int ij, cc, ym, mcase, ecase;
    double **m_y, **m_c, err_el[2], *m_i, *A, **C_g;

    mcase = 1; // 1 = original MYC		2 = error check MYC

    ecase = 2; // case of error calculation

    m_i = malloc(sizeof(double) * 2);
    C_g = malloc(sizeof(double *) * (n + 2));
    A = malloc(sizeof(double) * 6);

    for (int j = 0; j < n + 2; j++)
    {
        C_g[j] = malloc(sizeof(double *) * (m + 2));
    }

    m_y = Youngs_method(C, S, m, n, dx0); // row output: my_ip = [my_xp; my_yp]

    m_c = CC_method(C, S, m, n, dx0); // row outputs: mci_jp = [mcy_xp; mcy_yp; mcx_xp; mcx_yp]

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            C_g[0][i + 1] = C[0][i]; // adding homogenous ghost values to color function
            C_g[m + 1][i + 1] = C[m - 1][i];
            C_g[j + 1][0] = C[j][0];
            C_g[j + 1][n + 1] = C[j][n - 1];
            C_g[j + 1][i + 1] = C[j][i];
        }
    }

    C_g[0][0] = 0.5 * (C_g[0][1] + C_g[1][0]); // Corner Ghost values
    C_g[0][n + 1] = 0.5 * (C_g[1][n + 1] + C_g[0][n]);
    C_g[m + 1][n + 1] = 0.5 * (C_g[m][n + 1] + C_g[m + 1][n]);
    C_g[m + 1][0] = 0.5 * (C_g[m][0] + C_g[m + 1][1]);

    cc = 0;
    ym = 0; // initialzing counters

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {

            if (C[j][i] > 1e-10 && C[j][i] < 1.0 - 1e-10)
            {

                ij = i + j * m;

                if (i == 0 || i == n - 1 || j == 0 || j == m - 1) // if on cell boundary, use youngs method. Guarantees proper boundary conditions
                {

                    S[0][ij] = m_y[0][ij]; // x1
                    S[1][ij] = m_y[1][ij]; // x2
                    S[2][ij] = m_y[2][ij]; // y1
                    S[3][ij] = m_y[3][ij]; // y2
                    S[4][ij] = m_y[4][ij]; // n_x
                    S[5][ij] = m_y[5][ij]; // n_y
                    S[6][ij] = 2.0;        // indicating Youngs' method

                    ym++;
                }
                else
                {

                    if (mcase == 1)
                    {
                        if (m_c[6][ij] == 1.0)
                        { // y representation from CC

                            if (fabs(m_c[5][ij]) < fabs(m_y[5][ij]))
                            {

                                S[0][ij] = m_c[0][ij]; // x1
                                S[1][ij] = m_c[1][ij]; // x2
                                S[2][ij] = m_c[2][ij]; // y1
                                S[3][ij] = m_c[3][ij]; // y2
                                S[4][ij] = m_c[4][ij]; // n_x
                                S[5][ij] = m_c[5][ij]; // n_y
                                S[6][ij] = 1.0;        // indicating centered columns

                                cc++;
                            }
                            else
                            {

                                S[0][ij] = m_y[0][ij]; // x1
                                S[1][ij] = m_y[1][ij]; // x2
                                S[2][ij] = m_y[2][ij]; // y1
                                S[3][ij] = m_y[3][ij]; // y2
                                S[4][ij] = m_y[4][ij]; // n_x
                                S[5][ij] = m_y[5][ij]; // n_y
                                S[6][ij] = 2.0;        // indicating Youngs' method

                                ym++;
                            }
                        }
                        else
                        { // x representation from cc
                            if (fabs(m_c[4][ij]) < fabs(m_y[4][ij]))
                            {

                                S[0][ij] = m_c[0][ij]; // x1
                                S[1][ij] = m_c[1][ij]; // x2
                                S[2][ij] = m_c[2][ij]; // y1
                                S[3][ij] = m_c[3][ij]; // y2
                                S[4][ij] = m_c[4][ij]; // n_x
                                S[5][ij] = m_c[5][ij]; // n_y
                                S[6][ij] = 1.0;        // indicating centered columns method

                                cc++;
                            }
                            else
                            {

                                S[0][ij] = m_y[0][ij]; // x1
                                S[1][ij] = m_y[1][ij]; // x2
                                S[2][ij] = m_y[2][ij]; // y1
                                S[3][ij] = m_y[3][ij]; // y2
                                S[4][ij] = m_y[4][ij]; // n_x
                                S[5][ij] = m_y[5][ij]; // n_y
                                S[6][ij] = 2.0;        // indicating Youngs' Method

                                ym++;
                            }
                        }
                    }
                    else
                    {

                        A[0] = m_c[0][ij]; // x1 for current cell using centered columns method
                        A[1] = m_c[1][ij]; // x2 "                                             "
                        A[2] = m_c[2][ij]; // y1 "												"
                        A[3] = m_c[3][ij]; // y2 "												"
                        A[4] = m_c[4][ij]; // n_x "											"
                        A[5] = m_c[5][ij]; // n_y "											"

                        err_el[0] = error_el(C_g, A, dx0, i, j, ecase); // error "					"

                        A[0] = m_y[0][ij]; // x1 for current cell using Youngs' Method
                        A[1] = m_y[1][ij]; // x2 "										"
                        A[2] = m_y[2][ij]; // y1 "										"
                        A[3] = m_y[3][ij]; // y2 "										"
                        A[4] = m_y[4][ij]; // n_x "									"
                        A[5] = m_y[5][ij]; // n_y "									"

                        err_el[1] = error_el(C_g, A, dx0, i, j, ecase);

                        if (err_el[0] < err_el[1])
                        {

                            S[0][ij] = m_c[0][ij]; // x1
                            S[1][ij] = m_c[1][ij]; // x2
                            S[2][ij] = m_c[2][ij]; // y1
                            S[3][ij] = m_c[3][ij]; // y2
                            S[4][ij] = m_c[4][ij]; // n_x
                            S[5][ij] = m_c[5][ij]; // n_y
                            S[6][ij] = 1.0;        // indicating centered columns method
                        }
                        else
                        {

                            S[0][ij] = m_y[0][ij]; // x1
                            S[1][ij] = m_y[1][ij]; // x2
                            S[2][ij] = m_y[2][ij]; // y1
                            S[3][ij] = m_y[3][ij]; // y2
                            S[4][ij] = m_y[4][ij]; // n_x
                            S[5][ij] = m_y[5][ij]; // n_y
                            S[6][ij] = 2.0;        // indicating Youngs' method
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < m + 2; i++)
    {
        free(C_g[i]);
    }

    free(m_i);
    free(C_g);
    free(A);

    return S;
}

double **ELVIRA(double **C, double **S, double dx0, int m, int n)
{

    int ij, p, ecase;
    double dx0i, dx0i2, hgt_x[m][n], hgt_y[m][n], mc_x, mc_y, sgn_mc_x, sgn_mc_y, **C_g;
    double mc_xp, mc_yp, m_i[2], *A0, A1[6][6], err_el[6], m_L;

    ecase = 1; // case for error calculation

    C_g = malloc(sizeof(double *) * (n + 2));
    for (int j = 0; j < n + 2; j++)
    {
        C_g[j] = malloc(sizeof(double *) * (m + 2));
    }

    dx0i = 1.0 / dx0; // inverting dx to prevent repeated division (for computational speed)
    dx0i2 = 1.0 / (2.0 * dx0);

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            C_g[0][i + 1] = C[m - 1][i]; // adding periodic ghost value to color function
            C_g[m + 1][i + 1] = C[0][i];
            C_g[j + 1][0] = C[j][n - 1];
            C_g[j + 1][n + 1] = C[j][0];
            C_g[j + 1][i + 1] = C[j][i];
        }
    }

    C_g[0][0] = 0.5 * (C_g[0][1] + C_g[1][0]); // Corner Ghost values
    C_g[0][n + 1] = 0.5 * (C_g[1][n + 1] + C_g[0][n]);
    C_g[m + 1][n + 1] = 0.5 * (C_g[m][n + 1] + C_g[m + 1][n]);
    C_g[m + 1][0] = 0.5 * (C_g[m][0] + C_g[m + 1][1]);

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            if (C[j][i] > 1e-10 && C[j][i] < 1.0 - 1e-10)
            {

                ij = i + j * m; // linearizing matrix

                for (int k = 0; k < 6; k++)
                { // loop to use produce 6 interface reconstructions

                    if (k == 0)
                    { // (mx)x = (my)y + b	central difference scheme

                        mc_y = 0.5 * ((C_g[j + 2][i] + C_g[j + 2][i + 1] + C_g[j + 2][i + 2]) - (C_g[j][i] + C_g[j][i + 1] + C_g[j][i + 2]));   // mcy = d(hgt_x)dy
                        mc_x = (C_g[j + 1][i + 2] - C_g[j + 1][i]);     // mcx = dC/dx
                    }
                    else if (k == 1)
                    { // (mx)x = (my)y + b 	backwards difference scheme

                        mc_y = ((C_g[j + 1][i] + C_g[j + 1][i + 1] + C_g[j + 1][i + 2]) - (C_g[j][i] + C_g[j][i + 1] + C_g[j][i + 2]));
                        mc_x = (C_g[j + 1][i + 1] - C_g[j + 1][i]);
                    }
                    else if (k == 2)
                    { // (mx)x = (my)y + b 	forward difference scheme

                        mc_y = ((C_g[j + 2][i] + C_g[j + 2][i + 1] + C_g[j + 2][i + 2]) - (C_g[j + 1][i] + C_g[j + 1][i + 1] + C_g[j + 1][i + 2]));
                        mc_x = (C_g[j + 1][i + 2] - C_g[j + 1][i + 1]);
                    }
                    else if (k == 3)
                    { // (my)y = (mx)x + b 	central difference scheme

                        mc_x = 0.5 * ((C_g[j][i + 2] + C_g[j + 1][i + 2] + C_g[j + 2][i + 2]) - (C_g[j][i] + C_g[j + 1][i] + C_g[j + 2][i]));   // mcx = d(hgt_y)/dx
                        mc_y = (C_g[j + 2][i + 1] - C_g[j][i + 1]);     // mcy = dC/dy
                    }
                    else if (k == 4)
                    { // (my)y = (mx)x + b 	backwards difference scheme

                        mc_x = ((C_g[j][i + 1] + C_g[j + 1][i + 1] + C_g[j + 2][i + 1]) - (C_g[j][i] + C_g[j + 1][i] + C_g[j + 2][i]));
                        mc_y = (C_g[j + 1][i + 1] - C_g[j][i + 1]);
                    }
                    else
                    { // (my)y = (mx)x + b 		forward difference schme

                        mc_x = ((C_g[j][i + 2] + C_g[j + 1][i + 2] + C_g[j + 2][i + 2]) - (C_g[j][i + 1] + C_g[j + 1][i + 1] + C_g[j + 2][i + 1]));
                        mc_y = (C_g[j + 2][i + 1] - C_g[j + 1][i + 1]);
                    }

                    if (k < 3)
                    {

                        if (mc_x > 0.0)
                        {
                            sgn_mc_x = 1.0; // calculating sgn of mcx_x
                        }
                        else
                        {
                            sgn_mc_x = -1.0;
                        }

                        mc_xp = sgn_mc_x / (1.0 + fabs(mc_y)); // "normalizing" mcx_x and mcx_y so
                        mc_yp = mc_y / (1.0 + fabs(mc_y));     // ||mcx_x|| + ||mcx_y|| = 1
                    }
                    else
                    {

                        if (mc_y > 0.0)
                        {
                            sgn_mc_y = 1.0; // calculating sgn of mcx_x
                        }
                        else
                        {
                            sgn_mc_y = -1.0;
                        }

                        mc_xp = mc_x / (1.0 + fabs(mc_x));     // "normalizing" mcx_x and mcx_y so
                        mc_yp = sgn_mc_y / (1.0 + fabs(mc_x)); // ||mcx_x|| + ||mcx_y|| = 1
                    }

                    if (fabs(mc_yp) < 1e-10)
                    { // ensuring slope never approaches infinity

                        if (mc_yp > 0.0)
                        {
                            mc_yp = 1.0e-10;
                        }
                        else
                        {
                            mc_yp = -1.0e-10;
                        }

                        mc_xp = 0.5; // making results so off target that it wont be used
                        mc_yp = 0.5;
                    }

                    m_i[0] = mc_xp;
                    m_i[1] = mc_yp;

                    /// printf("mx = %f\t my = %f\n",m_i[0],m_i[1]);

                    A0 = VoF_tools(C[j][i], m_i, dx0, m, n); // interface reconstruction

                    A1[0][k] = A0[0]; // saving interface reconstruction options
                    A1[1][k] = A0[1];
                    A1[2][k] = A0[2];
                    A1[3][k] = A0[3];
                    A1[4][k] = A0[4];
                    A1[5][k] = A0[5];

                    free(A0);

                    err_el[k] = error_el(C_g, A0, dx0, i, j, ecase);
                    //	printf("%i\t %i\t err_el = %f\n",i,j,err_el[k]);
                }

                p = 0; // indexing p as first position of A1

                for (int k = 0; k < 6; k++)
                {
                    if (err_el[k] < err_el[p])
                    {
                        p = k;
                    }
                }

                S[0][ij] = A1[0][p];
                S[1][ij] = A1[1][p];
                S[2][ij] = A1[2][p];
                S[3][ij] = A1[3][p];
                S[4][ij] = A1[4][p];
                S[5][ij] = A1[5][p];
            }
        }
    }

    for (int i = 0; i < m + 2; i++)
    {
        free(C_g[i]);
    }

    free(C_g);

    return S;
}

double **PLIR(double **C, double **S, double dx0, int m, int n, int ICM)
{ // interface reconstruction function

    if (ICM == 1)
    { // centered columns method

        S = CC_method(C, S, m, n, dx0); // row outputs: mci_jp = [x1; x2; y1; y2; n_x; n_y]
    }
    else if (ICM == 2)
    { // Youngs Method

        S = Youngs_method(C, S, m, n, dx0); // row outputs: mci_jp = [x1; x2; y1; y2; n_x; n_y]
    }
    else if (ICM == 3)
    { // Mixed Youngs Centered Method

        S = myc_method(C, S, m, n, dx0);
    }
    else if (ICM == 4)
    { // ELVIRA

        S = ELVIRA(C, S, dx0, m, n); // row outputs: mci_jp = [x1; x2; y1; y2; n_x; n_y]
    }

    return S;
}

double F_xp_le(double **S, double u_bar_l, double u_bar_r, double dx0i, int m, int n, int i, int j, int bc)
{

    int ij;
    double F, n_x, n_y, m_p, x1, x2, y1, y2, x1p, x2p, y1p, y2p, y_i;

    if (bc == 0)
    {
        ij = i - 1 + j * m; // interior indexed cell
    }
    else
    {
        ij = (n - 1) + j * m; // left or right boundary indexed cell
    }

    n_x = S[4][ij];
    n_y = S[5][ij];

    if (S[0][ij] >= S[1][ij])
    { // assigning x1, y1 to donwnind values
        x1 = S[0][ij];
        x2 = S[1][ij];
        y1 = S[2][ij];
        y2 = S[3][ij];
    }
    else
    {
        x1 = S[1][ij];
        x2 = S[0][ij];
        y1 = S[3][ij];
        y2 = S[2][ij];
    }

    x1p = (1.0 + u_bar_r - u_bar_l) * x1 + u_bar_l; // advecting interface in x direction
    x2p = (1.0 + u_bar_r - u_bar_l) * x2 + u_bar_l;
    y1p = y1;
    y2p = y2;

    // calculating right fluxes
    // All fluxes are calculated as a volume fraction with respect to cell area
    if (x1p > 1.0 && x2p > 1.0)
    { // x1 and x2 out of cell
        if (y1p > 1e-10 && y1p < 1.0 - 1e-10)
        { // y1 not top or bottom bndry

            if (n_x > 0.0)
            {
                F = (0.5 * (x1p - x2p) * fabs(y1p - y2p));
            }
            else
            {
                F = (u_bar_r - 0.5 * (x1p - x2p) * fabs(y1p - y2p));
            }
        }
        else
        { // y1 on top or bottom boudary

            if (n_x > 0.0)
            {
                F = 0.5 * ((1.0 + u_bar_r - x1p) + (1.0 + u_bar_r - x2p)); // trapizoid formula
            }
            else
            {
                F = 0.5 * ((x1p - 1.0) + (x2p - 1.0)); // trapizoid fomrula
            }
        }
        //		printf("testing01\n");
    }
    else if (x1p >= 1.0 && x2p <= 1.0)
    { // x1 out of cell; x2 in the  cell

        m_p = (y1p - y2p) / (x1p - x2p); // slope of interface after advection
        y_i = y1p - m_p * (x1p - 1.0);   // y intercept at cell interface; needed for flux

        if (y1p > 1e-10 && y1p < 1.0 - 1e-10)
        { // y1 not top or bottom bndry

            if (n_y > 0.0)
            {

                F = 0.5 * ((1.0 - y1p) + (1.0 - y_i)) * u_bar_r; // trapizoid formula
            }
            else
            {

                F = 0.5 * (y1p + y_i) * u_bar_r; // trapizoid fomrula
            }
        }
        else
        { // y1 on top or bottom boundary

            if (n_x > 0.0)
            {

                F = (u_bar_r - 0.5 * (x1p - 1.0) * fabs(y1p - y_i));
            }
            else
            {

                F = 0.5 * (x1p - 1.0) * fabs(y1p - y_i);
            }
        }
    }
    else
    {

        if (n_x > 0.0 && u_bar_r > 0.0)
        {
            F = u_bar_r;
        }
        else
        {
            F = 0.0;
        }
    }

    return F;
}

double F_xn_le(double **S, double u_bar_l, double u_bar_r, double dx0i, int m, int n, int i, int j, int bc)
{

    int ij;
    double F, n_x, n_y, m_p, x1, x2, y1, y2, x1p, x2p, y1p, y2p, y_i, u_bar_lp, u_bar_rp;

    if (bc == 0)
    {
        ij = i + j * m; // interior indexed cell
    }
    else
    {
        ij = 0 + j * m; // left or right boundary indexed cell
    }

    n_x = S[4][ij];
    n_y = S[5][ij];

    if (S[0][ij] < S[1][ij])
    {
        x1 = S[0][ij];
        x2 = S[1][ij];
        y1 = S[2][ij];
        y2 = S[3][ij];
    }
    else
    {
        x1 = S[1][ij];
        x2 = S[0][ij];
        y1 = S[3][ij];
        y2 = S[2][ij];
    }

    x1p = (1.0 + u_bar_r - u_bar_l) * x1 + u_bar_l; // advecting PLIR in negative x direction
    x2p = (1.0 + u_bar_r - u_bar_l) * x2 + u_bar_l;
    y1p = y1;
    y2p = y2;

    u_bar_lp = fabs(u_bar_l); // converting to posative values to simplify following calculations
    u_bar_rp = fabs(u_bar_r);

    // calculating  left fluxes when one or more points leave the cell
    // All fluxes are calculated as area leaving the cell divided by dx0^2
    if (x1p < 0.0 && x2p < 0.0)
    { // x1 and x2 out of cell
        if (y1p > 1e-10 && y1p < 1.0 - 1e-10)
        { // y1 not top or bottom bndry

            if (n_x > 0.0)
            {
                F = (u_bar_lp - 0.5 * (x2p - x1p) * fabs(y1p - y2p));
            }
            else
            {
                F = 0.5 * (x2p - x1p) * fabs(y1p - y2p);
            }
        }
        else
        { // y1 on top or bottom boudary

            if (n_x > 0.0)
            {
                F = 0.5 * (-1.0 * (x1p + x2p)); // trapizoid formula
            }
            else
            {
                F = 0.5 * ((x1p + u_bar_lp) + (x2p + u_bar_lp)); // trapizoid formula
            }
        }
    }
    else if (x1p <= 0.0 && x2p >= 0.0)
    { // x1 out of cell; x2 in the  cell

        m_p = (y2p - y1p) / (x2p - x1p); // slope of interface after advection
        y_i = y1p + m_p * (-1.0 * x1p);  // y intercept at cell interface; needed for flux

        if (y1p > 1e-10 && y1p < 1.0 - 1e-10)
        { // y1 not top or bottom bndry

            if (n_y > 0.0)
            {
                F = 0.5 * ((1.0 - y1p) + (1.0 - y_i)) * u_bar_lp; // trapazoid
            }
            else
            {
                F = 0.5 * (y1p + y_i) * u_bar_lp; // trazoid formula
            }
        }
        else
        { // y1 on top or bottom boundary

            if (n_x > 0.0)
            {
                F = 0.5 * fabs(x1p * (y1p - y_i));
            }
            else
            {
                F = (u_bar_lp - 0.5 * fabs(x1p * (y1p - y_i)));
            }
        }
    }
    else
    {

        if (n_x < 0.0 && u_bar_l < 0.0)
        {
            F = u_bar_lp;
        }
        else
        {
            F = 0.0;
        }
    }

    return F;
}

double F_yp_le(double **S, double v_bar_b, double v_bar_t, double dx0i, int m, int n, int i, int j, int bc)
{

    int ij;
    double F, n_x, n_y, m_p, x1, x2, y1, y2, x1p, x2p, y1p, y2p, x_i;

    if (bc == 0)
    {
        ij = i + (j - 1) * m; // interior indexed cell
    }
    else
    {
        ij = i + (m - 1) * m; // left or right boundary indexed cell
    }

    n_x = S[4][ij];
    n_y = S[5][ij];

    if (S[2][ij] >= S[3][ij])
    { // assignting (x1,y1) to the downwind cooridnate
        x1 = S[0][ij];
        x2 = S[1][ij];
        y1 = S[2][ij];
        y2 = S[3][ij];
    }
    else
    {
        x1 = S[1][ij];
        x2 = S[0][ij];
        y1 = S[3][ij];
        y2 = S[2][ij];
    }

    x1p = x1;
    x2p = x2;
    y1p = (1.0 + v_bar_t - v_bar_b) * y1 + v_bar_b; // advecting in the posative y direction
    y2p = (1.0 + v_bar_t - v_bar_b) * y2 + v_bar_b;

    // calculating top fluxes
    // All fluxes are calculated as volume fraction with respect to cell area
    if (y1p > 1.0 && y2p > 1.0)
    { // y1 out the cell; y2 out the cell

        if (x1p > 1e-10 && x1p < 1.0 - 1e-10)
        { // x1 not on top or bottom boundary
            if (n_y > 0.0)
            {
                F = (0.5) * (y1p - y2p) * fabs(x1p - x2p); // traiangle formula
            }
            else
            {
                F = (v_bar_t - (0.5) * (y1p - y2p) * fabs(x1p - x2p)); // rectangle minus triangle
            }
        }
        else
        { // x1 on top or bottom boundary
            if (n_y > 0.0)
            {
                F = (0.5) * ((1.0 + v_bar_t - y1p) + (1.0 + v_bar_t - y2p)); // trapizoid
            }
            else
            {
                F = (0.5) * ((y1p - 1.0) + (y2p - 1.0)); // trapizoid
            }
        }
    }
    else if (y1p >= 1.0 && y2p <= 1.0)
    { // y1 out the cell; y2 in the cell

        if (x1p > x2p)
        {
            m_p = (y1p - y2p) / (x1p - x2p); // slope of interface after advection
        }
        else
        {
            m_p = (y2p - y1p) / (x2p - x1p);
        }

        x_i = x1p - ((y1p - 1.0) / m_p); // x intercept at cell face; needed for flux

        if (x1p > 1e-10 && x1p < 1.0 - 1e-10)
        { // x1 not on left or right boudary
            if (n_x > 0.0)
            {
                F = (0.5) * ((1.0 - x1p) + (1.0 - x_i)) * v_bar_t; // trapizoid
            }
            else
            {
                F = (0.5) * (x1p + x_i) * v_bar_t; // trapizoid
            }
        }
        else
        { // x1 on top or bottom boundary

            if (n_y > 0.0)
            {
                F = (v_bar_t - 0.5 * fabs((x1p - x_i) * (y1p - 1.0))); // rectangle minus triangle
            }
            else
            {
                F = (0.5) * fabs((x1p - x_i) * (y1p - 1.0)); // triangle
            }
        }
    }
    else
    { // y1 in cell; y2 in cell

        if (n_y > 0.0 && v_bar_t > 0.0)
        {
            F = v_bar_t; // local CFL number for current cell (j,i)
        }
        else
        {
            F = 0.0;
        }
    }

    return F;
}

double F_yn_le(double **S, double v_bar_b, double v_bar_t, double dx0i, int m, int n, int i, int j, int bc)
{

    int ij;
    double F, n_x, n_y, m_p, x1, x2, y1, y2, x1p, x2p, y1p, y2p, x_i, v_bar_bp, v_bar_tp;

    if (bc == 0)
    {
        ij = i + j * m; // interior indexed cell
    }
    else
    {
        ij = i; // left or right boundary indexed cell
    }

    n_x = S[4][ij];
    n_y = S[5][ij];

    if (S[2][ij] <= S[3][ij])
    { // assignting (x1,y1) to the downwind cooridnate
        x1 = S[0][ij];
        x2 = S[1][ij];
        y1 = S[2][ij];
        y2 = S[3][ij];
    }
    else
    {
        x1 = S[1][ij];
        x2 = S[0][ij];
        y1 = S[3][ij];
        y2 = S[2][ij];
    }

    x1p = x1;
    x2p = x2;
    y1p = (1.0 + v_bar_t - v_bar_b) * y1 + v_bar_b; // advecting in the negative y direction
    y2p = (1.0 + v_bar_t - v_bar_b) * y2 + v_bar_b;

    v_bar_tp = fabs(v_bar_t); // converting v_bar_b and v_bar_t to posative values to simplify calculation
    v_bar_bp = fabs(v_bar_b);

    // calculating bottom fluxes
    // All fluxes are calculated as area fraction leaving cell
    if (y1p <= 0.0 && y2p < 0.0)
    { // y1 out the cell; y2 out the cell

        if (x1p > 1e-10 && x1p < 1.0 - 1e-10)
        { // x1 not on top or bottom boundary
            if (n_y < 0.0)
            {
                F = (0.5) * fabs((y2p - y1p) * (x1p - x2p)); // triangle formula
            }
            else
            {
                F = (v_bar_bp - (0.5) * fabs((y2p - y1p) * (x1p - x2p))); // rectangle minus triangle
            }
        }
        else
        { // x1 on top or bottom boundary
            if (n_y < 0.0)
            {
                F = (0.5) * ((y1p + v_bar_bp) + (y2p + v_bar_bp)); // trapizoid formula
            }
            else
            {
                F = (0.5) * (-1.0 * (y1p + y2p)); // trapizoid formula
            }
        }
    }
    else if (y1p <= 0.0 && y2p >= 0.0)
    { // y1 out the cell; y2 in the cell

        if (x1p > x2p)
        {
            m_p = (y1p - y2p) / (x1p - x2p); // slope of interface after advection
        }
        else
        {
            m_p = (y2p - y1p) / (x2p - x1p);
        }

        x_i = x1p - (y1p / m_p); // x intercept at cell interface; needed for flux

        if (x1p > 1e-10 && x1p < 1.0 - 1e-10)
        { // x1 not on top or bottom boudary
            if (n_x > 0.0)
            {
                F = (0.5) * ((1.0 - x1p) + (1.0 - x_i)) * v_bar_bp; // trapizoid formula
            }
            else
            {
                F = (0.5) * (x1p + x_i) * v_bar_bp; // trapizoid formula
            }
        }
        else
        { // x1 on top or bottom boundary
            if (n_y < 0.0)
            {
                F = (v_bar_bp - (0.5) * fabs(y1p * (x1p - x_i))); // rectangle minus triangle
            }
            else
            {
                F = (0.5) * fabs(y1p * (x1p - x_i)); // triangle
            }
        }
    }
    else
    { // y1 in cell; y2 in cell

        if (n_y < 0.0 && v_bar_b < 0.0)
        {
            F = v_bar_bp;
        }
        else
        {
            F = 0.0;
        }
    }

    return F;
}

double F_xp_ei(double **S, double u_bar_l, double u_bar_r, double dx0i, int m, int n, int i, int j, int bc)
{

    int ij;
    double F, n_x, n_y, m_p, x1, x2, y1, y2, y_i;

    if (bc == 0)
    {
        ij = i - 1 + j * m; // interior indexed cell
    }
    else
    {
        ij = (n - 1) + j * m; // left or right boundary indexed cell
    }

    n_x = S[4][ij];
    n_y = S[5][ij];

    if (S[0][ij] >= S[1][ij])
    {
        x1 = S[0][ij];
        x2 = S[1][ij];
        y1 = S[2][ij];
        y2 = S[3][ij];
    }
    else
    {
        x1 = S[1][ij];
        x2 = S[0][ij];
        y1 = S[3][ij];
        y2 = S[2][ij];
    }

    if (x1 > 1.0 - u_bar_r && x2 > 1.0 - u_bar_r)
    { // x1 in flux boundary; x2 in flux boundary
        if (y1 > 1e-10 && y1 < 1.0 - 1e-10)
        { // y1 not on top or bottom boundary
            if (n_x > 0.0)
            {
                F = 0.5 * fabs((x1 - x2) * (y1 - y2)); // triangle
            }
            else
            {
                F = (u_bar_r - 0.5 * fabs((x1 - x2) * (y1 - y2))); // rectangle minus triangle
            }
        }
        else
        { // y1 on top or bottom boundary
            if (n_x > 0.0)
            {
                F = 0.5 * ((1.0 - x1) + (1.0 - x2)); // trapizoid
            }
            else
            {
                F = 0.5 * ((x1 - (1.0 - u_bar_r)) + (x2 - (1.0 - u_bar_r))); // trapizoid
            }
        }
    }
    else if (x1 > 1.0 - u_bar_r && x2 <= 1.0 - u_bar_r)
    { // x1 in flux boundary; x2 out flux boundary

        m_p = (y1 - y2) / (x1 - x2);             // slope of interface
        y_i = y1 - m_p * (x1 - (1.0 - u_bar_r)); // y intercept at flux boundary

        if (y1 > 1e-10 && y1 < 1.0 - 1e-10)
        { // y1 not on top or bottom boundary

            if (n_y > 0.0)
            {

                F = 0.5 * ((1.0 - y1) + (1.0 - y_i)) * u_bar_r; // trapizoid
            }
            else
            {

                F = 0.5 * (y1 + y_i) * u_bar_r; // trapizoid
            }
        }
        else
        { // y1 on top or bottom boundary
            if (n_x > 0.0)
            {
                F = (u_bar_r - 0.5 * fabs((x1 - (1.0 - u_bar_r)) * (y1 - y_i))); // rectangle minus triangle
            }
            else
            {
                F = 0.5 * fabs((x1 - (1.0 - u_bar_r)) * (y1 - y_i)); // triangle
            }
        }
    }
    else
    { // x1 out flux boundary; x2 out flux boundary

        if (n_x > 0.0)
        {
            F = u_bar_r;
        }
        else
        {
            F = 0.0;
        }
    }

    return F;
}

double F_xn_ei(double **S, double u_bar_l, double u_bar_r, double dx0i, int m, int n, int i, int j, int bc)
{

    int ij;
    double F, n_x, n_y, m_p, x1, x2, y1, y2, y_i, u_bar_lp, u_bar_rp;

    if (bc == 0)
    {
        ij = i + j * m; // interior indexed cell
    }
    else
    {
        ij = 0 + j * m; // left or right boundary indexed cell
    }

    n_x = S[4][ij]; // extractring values from PLIR
    n_y = S[5][ij];

    if (S[0][ij] <= S[1][ij])
    {
        x1 = S[0][ij];
        x2 = S[1][ij];
        y1 = S[2][ij];
        y2 = S[3][ij];
    }
    else
    {
        x1 = S[1][ij];
        x2 = S[0][ij];
        y1 = S[3][ij];
        y2 = S[2][ij];
    }

    u_bar_lp = fabs(u_bar_l); // converting to posative value to simplify calculation
    u_bar_rp = fabs(u_bar_r);

    if (x1 < u_bar_lp && x2 < u_bar_lp)
    { // x1 in flux boundary; x2 in flux boundary

        if (y1 > 1e-10 && y1 < 1.0 - 1e-10)
        { // y1 not on top or bottom boundary of cell
            if (n_x < 0.0)
            {
                F = 0.5 * fabs((x1 - x2) * (y1 - y2)); // Triangle
            }
            else
            {
                F = (u_bar_lp - 0.5 * fabs((x1 - x2) * (y1 - y2))); // rectangle minus triangle
            }
        }
        else
        {
            if (n_x < 0.0)
            {
                F = 0.5 * (x1 + x2); // trapizoid
            }
            else
            {
                F = 0.5 * ((u_bar_lp - x1) + (u_bar_lp - x2)); // trapizoid
            }
        }
    }
    else if (x1 <= u_bar_lp && x2 >= u_bar_lp)
    { // x1 in flux boundary; x2 out flux boundary

        m_p = (y2 - y1) / (x2 - x1);      // slope of interface
        y_i = y1 + m_p * (u_bar_lp - x1); // y intercept at flux boundary

        if (y1 > 1e-10 && y1 < 1.0 - 1e-10)
        {
            if (n_y > 0.0)
            {
                F = 0.5 * ((1.0 - y1) + (1.0 - y_i)) * u_bar_lp; // trapizoid
            }
            else
            {
                F = 0.5 * (y1 + y_i) * u_bar_lp; // trapizoid
            }
        }
        else
        {
            if (n_x < 0.0)
            {
                F = (u_bar_lp - 0.5 * fabs((x1 - u_bar_lp) * (y1 - y_i))); // rectangle minus triangle
            }
            else
            {
                F = 0.5 * fabs((x1 - u_bar_lp) * (y1 - y_i)); // triangle
            }
        }
    }
    else
    { // x1 out flux boundary; x2 out flux boundary

        if (n_x < 0.0)
        {
            F = u_bar_lp;
        }
        else
        {
            F = 0.0;
        }
    }

    return F;
}

double F_yp_ei(double **S, double v_bar_b, double v_bar_t, double dx0i, int m, int n, int i, int j, int bc)
{

    int ij;
    double F, n_x, n_y, m_p, x1, x2, y1, y2, x_i;

    if (bc == 0)
    {
        ij = i + (j - 1) * m; // interior indexed cell
    }
    else
    {
        ij = i + (m - 1) * m; // left or right boundary indexed cell
    }

    n_x = S[4][ij];
    n_y = S[5][ij];

    if (S[2][ij] >= S[3][ij])
    { // assignting (x1,y1) to the downwind cooridnate
        x1 = S[0][ij];
        x2 = S[1][ij];
        y1 = S[2][ij];
        y2 = S[3][ij];
    }
    else
    {
        x1 = S[1][ij];
        x2 = S[0][ij];
        y1 = S[3][ij];
        y2 = S[2][ij];
    }

    // calculating fluxes
    // All fluxes are calcualted as area fraction leaving cell
    if (y1 >= 1.0 - v_bar_t && y2 > 1.0 - v_bar_t)
    { // y1 in flux boundary; y2 in flux boundary

        if (x1 > 1e-10 && x1 < 1.0 - 1e-10)
        { // x1 not on top or bottom boundary of cell
            if (n_y > 0.0)
            {
                F = 0.5 * fabs((x1 - x2) * (y1 - y2)); // triangle
            }
            else
            {
                F = (v_bar_t - 0.5 * fabs((x1 - x2) * (y1 - y2))); // rectangle minus triangle
            }
        }
        else
        { // x1 on top or bottom boundary of cell
            if (n_y > 0.0)
            {
                F = 0.5 * ((1.0 - y1) + (1.0 - y2)); // trapizoid
            }
            else
            {
                F = 0.5 * (y1 - (1.0 - v_bar_t) + (y2 - (1.0 - v_bar_t))); // trapizoid
            }
        }
    }
    else if (y1 >= 1.0 - v_bar_t && y2 <= 1.0 - v_bar_t)
    { // y1 in flux boundary; y2 out flux boundary

        if (x1 > x2)
        {
            m_p = (y1 - y2) / (x1 - x2); // slope of interface
        }
        else
        {
            m_p = (y2 - y1) / (x2 - x1);
        }

        x_i = x1 + ((1.0 - v_bar_t - y1) / m_p); //  x intercept at flux boundary; needed for flux calculation

        if (x1 > 1e-10 && x1 < 1.0 - 1e-10)
        { // x1 not on left or right boundary of cell

            if (n_x > 0.0)
            {
                F = 0.5 * ((1.0 - x1) + (1.0 - x_i)) * v_bar_t; // trapizoid
            }
            else
            {
                F = 0.5 * (x1 + x_i) * v_bar_t; // trapizoid
            }
        }
        else
        { // x1 on top or bottom boundary

            if (n_y > 0.0)
            {
                F = (v_bar_t - 0.5 * fabs((x1 - x_i) * (y1 - (1.0 - v_bar_t)))); // rectangle minus triangle
            }
            else
            {
                F = 0.5 * fabs((x1 - x_i) * (y1 - (1.0 - v_bar_t))); // triangle
            }
        }
    }
    else
    { // y1 out flux boundary; y2 out flux boundary

        if (n_y > 0.0)
        {
            F = v_bar_t;
        }
        else
        {
            F = 0.0;
        }
    }

    return F;
}

double F_yn_ei(double **S, double v_bar_b, double v_bar_t, double dx0i, int m, int n, int i, int j, int bc)
{

    int ij;
    double F, n_x, n_y, m_p, x1, x2, y1, y2, x_i, v_bar_bp;

    if (bc == 0)
    {
        ij = i + j * m; // interior indexed cell
    }
    else
    {
        ij = i; // left or right boundary indexed cell
    }

    n_x = S[4][ij];
    n_y = S[5][ij];

    if (S[2][ij] <= S[3][ij])
    { // assignting (x1,y1) to the downwind cooridnate
        x1 = S[0][ij];
        x2 = S[1][ij];
        y1 = S[2][ij];
        y2 = S[3][ij];
    }
    else
    {
        x1 = S[1][ij];
        x2 = S[0][ij];
        y1 = S[3][ij];
        y2 = S[2][ij];
    }

    v_bar_bp = fabs(v_bar_b); // converting to posative values to simplfy calculation

    // calculating fluxes 		// All fluxes are calculated as area fraction leaving the cell
    if (y1 <= v_bar_bp && y2 < v_bar_bp)
    { // y1 in flux boundary; y2 in flux boundary

        if (x1 > 1e-10 && x1 < 1.0 - 1e-10)
        { // x1 not on left or right boundary of cell
            if (n_y < 0.0)
            {
                F = 0.5 * fabs((x1 - x2) * (y1 - y2)); // triangle
            }
            else
            {
                F = (v_bar_bp - 0.5 * fabs((x1 - x2) * (y1 - y2))); // rectangle minus triangle
            }
        }
        else
        {
            if (n_y < 0.0)
            {
                F = 0.5 * (y1 + y2); // trapizoid
            }
            else
            {
                F = 0.5 * ((v_bar_bp - y1) + (v_bar_bp - y2)); // trapizoid
            }
        }
    }
    else if (y1 <= v_bar_bp && y2 >= v_bar_bp)
    { // y1 in flux boundar; y2 out flux boundary

        if (x1 > x2)
        {
            m_p = (y1 - y2) / (x1 - x2); // slope of interface
        }
        else
        {
            m_p = (y2 - y1) / (x2 - x1);
        }

        x_i = x1 - ((y1 - v_bar_bp) / m_p); // x intercept at flux boundary

        if (x1 > 1e-10 && x1 < 1.0 - 1e-10)
        { // x1 not on left or right boundary of cell
            if (n_x > 0.0)
            {
                F = 0.5 * ((1.0 - x1) + (1.0 - x_i)) * v_bar_bp; // trapizoid
            }
            else
            {
                F = 0.5 * (x1 + x_i) * v_bar_bp; // trapizoid
            }
        }
        else
        { // x1 on left or right boundary
            if (n_y < 0.0)
            {
                F = (v_bar_bp - 0.5 * fabs((x1 - x_i) * (y1 - v_bar_bp))); // rectangle minus triangle
            }
            else
            {
                F = 0.5 * fabs((x1 - x_i) * (y1 - v_bar_bp)); // triangle
            }
        }
    }
    else
    { // y1 out flux boundary; y2 out flux boundary

        if (n_y < 0.0)
        {
            F = v_bar_bp;
        }
        else
        {
            F = 0.0;
        }
    }

    return F;
}

double **advec_x(double **C, double **S, double dx0, int m, int n, double **u, double dt, int t, int sc)
{

    double dx0i, F_l[m][n], F_r[m][n], u_bar, u_bar_l, u_bar_r, bc;

    dx0i = 1.0 / dx0; // recipricol of dx0; used for computational speed to avoid division

    for (int j = 0; j < m; j++)
    { // looping over cell faces
        for (int i = 0; i < n + 1; i++)
        {

            if (i == 0 || i == n)
            {

                u_bar = 0.5 * dt * dx0i * (u[j][0] + u[j][n - 1]); // period boundary conditions
            }
            else
            {

                u_bar = 0.5 * dt * dx0i * (u[j][i - 1] + u[j][i]); // interpolated CFL numbers at cell face
            }

            if (u_bar > 1e-10)
            { // calculate fluxes leaving cell to the left

                if (i == 0)
                { // periodic boundary conditions

                    bc = 1; // left boundary condition

                    u_bar_r = u_bar; // u_bar_r  and u_bar_l needed for mapping interface points
                    u_bar_l = (dt * 0.5 * dx0i) * (u[j][n - 2] + u[j][n - 1]);

                    if (C[j][n - 1] < 1e-10)
                    { // empty cells (C = 0)

                        F_l[j][0] = 0.0;
                    }
                    else if (C[j][n - 1] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_l[j][0] = u_bar_r;
                    }
                    else
                    { // inteface cells

                        if (sc == 1)
                        {
                            F_l[j][0] = F_xp_le(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Lagrangian Explicit
                        }
                        else
                        {
                            F_l[j][0] = F_xp_ei(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Eularian Implicit
                        }
                    }
                }
                else if (i == 1)
                {

                    bc = 0; // interior point

                    u_bar_r = u_bar;
                    u_bar_l = (dt * 0.5 * dx0i) * (u[j][0] + u[j][n - 1]);

                    if (C[j][0] < 1e-10)
                    { // empty cells (C = 0)

                        F_l[j][1] = 0.0;
                        F_r[j][0] = 0.0;
                    }
                    else if (C[j][0] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_l[j][1] = u_bar_r;
                        F_r[j][0] = -u_bar_r;
                    }
                    else
                    { // interface cells

                        if (sc == 1)
                        {
                            F_l[j][1] = F_xp_le(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Lagrangian Explicit
                        }
                        else
                        {
                            F_l[j][1] = F_xp_ei(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Eulerian Implicit
                        }

                        F_r[j][0] = -F_l[j][1];
                    }
                }
                else if (i == n)
                {

                    bc = 1; // right boundary conditon

                    u_bar_r = u_bar; // u_bar_r  and u_bar_l needed for mapping interface points
                    u_bar_l = (dt * 0.5 * dx0i) * (u[j][n - 2] + u[j][n - 1]);

                    if (C[j][n - 1] < 1e-10)
                    { // empty cells (C = 0)

                        F_r[j][n - 1] = 0.0;
                    }
                    else if (C[j][n - 1] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_r[j][n - 1] = -u_bar_r;
                    }
                    else
                    { // inteface cells

                        F_r[j][n - 1] = -F_l[j][0];
                    }
                }
                else
                {

                    bc = 0; // interior point

                    u_bar_r = u_bar; // u_bar_r  and u_bar_l needed for mapping interface points
                    u_bar_l = (dt * 0.5 * dx0i) * (u[j][i - 2] + u[j][i - 1]);

                    if (C[j][i - 1] < 1e-10)
                    { // empty cells (C = 0)

                        F_l[j][i] = 0.0;
                        F_r[j][i - 1] = 0.0;
                    }
                    else if (C[j][i - 1] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_l[j][i] = u_bar_r;
                        F_r[j][i - 1] = -u_bar_r;
                    }
                    else
                    { // inteface cell

                        if (sc == 1)
                        {
                            F_l[j][i] = F_xp_le(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Lagrangian Explicit
                        }
                        else
                        {
                            F_l[j][i] = F_xp_ei(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Eulerian Implicit
                        }

                        F_r[j][i - 1] = -F_l[j][i];
                    }
                }
            }
            else if (u_bar < -1e-10)
            { // calculate fluxes leaving cell to the right

                if (i == 0)
                {

                    bc = 1;

                    u_bar_l = u_bar;
                    u_bar_r = (dt * 0.5 * dx0i) * (u[j][0] + u[j][1]);

                    if (C[j][0] < 1e-10)
                    { // empty cells (C = 0)

                        F_l[j][0] = 0.0;
                    }
                    else if (C[j][0] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_l[j][0] = u_bar_l; // note: this is a negative number
                    }
                    else
                    { // interface cells

                        if (sc == 1)
                        {
                            F_l[j][0] = -F_xn_le(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Lagrangian Explicit
                        }
                        else
                        {
                            F_l[j][0] = -F_xn_ei(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Eulerian Implicit
                        }
                    }
                }
                else if (i == n - 1)
                {

                    bc = 0;

                    u_bar_l = u_bar;
                    u_bar_r = (dt * 0.5 * dx0i) * (u[j][n - 1] + u[j][0]);

                    if (C[j][n - 1] < 1e-10)
                    { // empy cells (C = 0)

                        F_l[j][n - 1] = 0.0;
                        F_r[j][n - 2] = 0.0;
                    }
                    else if (C[j][n - 1] > 1.0 - 1e-10)
                    {

                        F_l[j][n - 1] = u_bar_l;
                        F_r[j][n - 2] = -u_bar_l;
                    }
                    else
                    {

                        if (sc == 1)
                        {
                            F_l[j][n - 1] = -F_xn_le(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Lagrangian Explicit
                        }
                        else
                        {
                            F_l[j][n - 1] = -F_xn_ei(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Eulerian Implicit
                        }

                        F_r[j][n - 2] = -F_l[j][n - 1];
                    }
                }
                else if (i == n)
                {

                    bc = 1;

                    u_bar_l = u_bar;
                    u_bar_r = (dt * 0.5 * dx0i) * (u[j][0] + u[j][1]);

                    if (C[j][0] < 1e-10)
                    { // empty cells (C = 0)

                        F_r[j][n - 1] = 0.0;
                    }
                    else if (C[j][0] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_r[j][n - 1] = -u_bar_l;
                    }
                    else
                    { // interface cells

                        F_r[j][n - 1] = -F_l[j][0];
                    }
                }
                else
                {

                    bc = 0;

                    u_bar_l = u_bar;
                    u_bar_r = (dt * 0.5 * dx0i) * (u[j][i] + u[j][i + 1]);

                    if (C[j][i] < 1e-10)
                    { // empty cells (C = 0)

                        F_l[j][i] = 0.0;
                        F_r[j][i - 1] = 0.0;
                    }
                    else if (C[j][i] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_l[j][i] = u_bar_l;
                        F_r[j][i - 1] = -u_bar_l;
                    }
                    else
                    { // interface cells

                        if (sc == 1)
                        {
                            F_l[j][i] = -F_xn_le(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Lagrangian Explicit
                        }
                        else
                        {
                            F_l[j][i] = -F_xn_ei(S, u_bar_l, u_bar_r, dx0i, m, n, i, j, bc); // Eulerian Implicit
                        }

                        F_r[j][i - 1] = -F_l[j][i];
                    }
                }
            }
            else
            { // u = 0

                if (i == 0)
                {

                    F_l[j][i] = 0.0;
                }
                else if (i == n)
                {

                    F_r[j][n - 1] = 0.0;
                }
                else
                {

                    F_l[j][i] = 0.0;
                    F_r[j][i - 1] = 0.0;
                }
            }
        }
    }

    for (int j = 0; j < m; j++)
    { // updating color function using fluxes
        for (int i = 0; i < n; i++)
        {

            if (i == 0)
            { // periodic boundary conditions

                u_bar_l = (dt * 0.5 * dx0i) * (u[j][i] + u[j][n - 1]);
                u_bar_r = (dt * 0.5 * dx0i) * (u[j][i + 1] + u[j][i]);
            }
            else if (i == n - 1)
            {

                u_bar_l = (dt * 0.5 * dx0i) * (u[j][i] + u[j][i - 1]);
                u_bar_r = (dt * 0.5 * dx0i) * (u[j][0] + u[j][i]);
            }
            else
            {

                u_bar_l = (dt * 0.5 * dx0i) * (u[j][i] + u[j][i - 1]); // interpolated CFL numbers at left and right boundaries
                u_bar_r = (dt * 0.5 * dx0i) * (u[j][i + 1] + u[j][i]);
            }

            if (sc == 1)
            {
                C[j][i] = (1.0 + u_bar_r - u_bar_l) * C[j][i] + F_r[j][i] + F_l[j][i]; // Lagrangian Explicit
            }
            else
            {
                C[j][i] = (C[j][i] + F_r[j][i] + F_l[j][i]) / (1.0 - (u_bar_r - u_bar_l)); // Eulerina Implicit
            }

            if (C[j][i] > -1e-10 && C[j][i] < 0.0)
            { // correcting for rounding errors in machine precision
                C[j][i] = 0.0;
            }
            else if (C[j][i] > 1.0 && C[j][i] < 1.0 + 1e-10)
            {
                C[j][i] = 1.0;
            }
            if (C[j][i] < 0.0 || C[j][i] > 1.0)
            { // terminating program if color function is out of bounds
                if (sc == 1)
                {
                    printf("for t = %i, at i = %i, j = %i, C = %.16f, which is out of bounds. Program needs to be debugged [x_le]\n", t, i, j, C[j][i]);
                    exit(0);
                }
                else
                {
                    printf("for t = %i, at i = %i, j = %i, C = %.16f, which is out of bounds. Program needs to be debugged [x_ei]\n", t, i, j, C[j][i]);
                    exit(0);
                }
            }
        }
    }

    return C;
}

double **advec_y(double **C, double **S, double dx0, int m, int n, double **v, double dt, int t, int sc)
{

    double dx0i, F_t[m][n], F_b[m][n], v_bar, v_bar_t, v_bar_b, bc;

    dx0i = 1.0 / dx0; // recipricol of dx0; used for computational speed to avoid division

    for (int j = 0; j < m + 1; j++)
    {
        for (int i = 0; i < n; i++)
        {

            if (j == 0 || j == m)
            {
                v_bar = 0.5 * dt * dx0i * (v[0][i] + v[m - 1][i]); // periodic boundary condtions
            }
            else
            {
                v_bar = 0.5 * dt * dx0i * (v[j - 1][i] + v[j][i]); // calculating velocites at cell faces
            }

            if (fabs(v_bar) > 1.0)
            {
                printf("|v_bar| > 1.0\t check flow parameters [x_le]\n");
                exit(0);
            }

            if (v_bar > 1e-10)
            { // calculate fluxes leaving bottom cell

                if (j == 0)
                { // periodic boundary conditions

                    bc = 1; // bottom boundary condition

                    v_bar_t = v_bar; // v_bar_r  and v_bar_l needed for mapping interface points
                    v_bar_b = (0.5 * dt * dx0i) * (v[m - 1][i] + v[m - 2][i]);

                    if (C[m - 1][i] < 1e-10)
                    { // empty cells (C = 0)

                        F_b[0][i] = 0.0;
                    }
                    else if (C[m - 1][i] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_b[0][i] = v_bar_t;
                    }
                    else
                    { // interface cells

                        if (sc == 1)
                        {
                            F_b[0][i] = F_yp_le(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc); // Lagranian Explicit
                        }
                        else
                        {
                            F_b[0][i] = F_yp_ei(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc); // Eulerian Implicit
                        }
                    }
                }
                else if (j == 1)
                {

                    bc = 0; // interior cell

                    v_bar_t = v_bar; // u_bar_r  and u_bar_l needed for mapping interface points
                    v_bar_b = (0.5 * dt * dx0i) * (v[0][i] + v[m - 1][i]);

                    if (C[0][i] < 1e-10)
                    { // empty cells (C = 0)

                        F_b[j][i] = 0.0;
                        F_t[j - 1][i] = 0.0;
                    }
                    else if (C[0][i] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_b[j][i] = v_bar_t;
                        F_t[j - 1][i] = -v_bar_t;
                    }
                    else
                    { // interface cells

                        if (sc == 1)
                        {
                            F_b[j][i] = F_yp_le(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }
                        else
                        {
                            F_b[j][i] = F_yp_ei(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }

                        F_t[j - 1][i] = -F_b[j][i];
                    }
                }
                else if (j == m)
                {

                    bc = 1; // top boundry condition

                    v_bar_t = v_bar;
                    v_bar_b = (0.5 * dt * dx0i) * (v[0][i] + v[m - 1][i]);

                    if (C[m - 1][i] < 1e-10)
                    { // empty cells (C = 0)

                        F_t[m - 1][i] = 0.0;
                    }
                    else if (C[m - 1][i] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_t[m - 1][i] = -v_bar_t;
                    }
                    else
                    { // interface cells

                        F_t[m - 1][i] = -F_b[0][i];
                    }
                }
                else
                {

                    bc = 0; // interior point

                    v_bar_t = v_bar;
                    v_bar_b = (0.5 * dt * dx0i) * (v[j - 1][i] + v[j - 2][i]);

                    if (C[j - 1][i] < 1e-10)
                    { // empty cells (C = 0)

                        F_b[j][i] = 0.0;
                        F_t[j - 1][i] = 0.0;
                    }
                    else if (C[j - 1][i] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_b[j][i] = v_bar_t;
                        F_t[j - 1][i] = -v_bar_t;
                    }
                    else
                    { // interface cells

                        if (sc == 1)
                        {
                            F_b[j][i] = F_yp_le(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }
                        else
                        {
                            F_b[j][i] = F_yp_ei(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }

                        F_t[j - 1][i] = -F_b[j][i];
                    }
                }
            }
            else if (v_bar < 1e-10)
            {

                if (j == 0)
                { // Periodic boundary conditons

                    bc = 1; // bottom boundary conditon

                    v_bar_b = v_bar;
                    v_bar_t = (0.5 * dx0i * dt) * (v[0][i] + v[1][i]);

                    if (C[0][i] < 1e-10)
                    { // empty cells (C = 0)

                        F_b[0][i] = 0.0;
                    }
                    else if (C[0][i] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_b[0][i] = v_bar_b; // note: this is a negative value
                    }
                    else
                    { // interface cells
                        if (sc == 1)
                        {
                            F_b[0][i] = -F_yn_le(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }
                        else
                        {
                            F_b[0][i] = -F_yn_ei(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }
                    }
                }
                else if (j == m - 1)
                {

                    bc = 0; // interior cell

                    v_bar_b = v_bar;
                    v_bar_t = (0.5 * dt * dx0i) * (v[m - 1][i] + v[0][i]);

                    if (C[j][i] < 1e-10)
                    { // empty cells (C = 0)

                        F_b[j][i] = 0.0;
                        F_t[j - 1][i] = 0.0;
                    }
                    else if (C[j][i] > 1.0 - 1e-10)
                    {

                        F_b[j][i] = v_bar_b; // note: this is a negative value
                        F_t[j - 1][i] = -v_bar_b;
                    }
                    else
                    { // interface cells

                        if (sc == 1)
                        {
                            F_b[j][i] = -F_yn_le(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }
                        else
                        {
                            F_b[j][i] = -F_yn_ei(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }
                        F_t[j - 1][i] = -F_b[j][i];
                    }
                }
                else if (j == m)
                {

                    bc = 1; // top cell boundary

                    v_bar_b = v_bar_b;
                    v_bar_t = (0.5 * dt * dx0i) * (v[0][i] + v[1][i]);

                    if (C[0][i] < 1e-10)
                    { // empty cells (C = 0)

                        F_t[m - 1][i] = 0.0;
                    }
                    else if (C[0][i] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_t[m - 1][i] = -v_bar_b; // note: this is a posative number
                    }
                    else
                    {
                        if (sc == 1)
                        {
                            F_t[m - 1][i] = F_yn_le(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }
                        else
                        {
                            F_t[m - 1][i] = F_yn_ei(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }
                    }
                }
                else
                {

                    bc = 0; // interior point

                    v_bar_b = v_bar;
                    v_bar_t = (0.5 * dt * dx0i) * (v[j][i] + v[j + 1][i]);

                    if (C[j][i] < 1e-10)
                    { // empty cells (C = 0)

                        F_b[j][i] = 0.0;
                        F_t[j - 1][i] = 0.0;
                    }
                    else if (C[j][i] > 1.0 - 1e-10)
                    { // full cells (C = 1)

                        F_b[j][i] = v_bar_b; // note: this is a negative number
                        F_t[j - 1][i] = -v_bar_b;
                    }
                    else
                    {

                        if (sc == 1)
                        {
                            F_b[j][i] = -F_yn_le(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }
                        else
                        {
                            F_b[j][i] = -F_yn_ei(S, v_bar_b, v_bar_t, dx0i, m, n, i, j, bc);
                        }

                        F_t[j - 1][i] = -F_b[j][i];
                    }
                }
            }
            else
            { // v = 0
                if (j == 0)
                {

                    F_b[j][i] = 0.0;
                }
                else if (i == n)
                {

                    F_t[j - 1][i] = 0.0;
                }
                else
                {

                    F_b[j][i] = 0.0;
                    F_t[j - 1][i] = 0.0;
                }
            }
        }
    }

    for (int j = 0; j < m; j++)
    { // updating color function using fluxes
        for (int i = 0; i < n; i++)
        {

            if (j == 0)
            { // period boundary conditons

                v_bar_b = (dt * 0.5 * dx0i) * (v[j][i] + v[m - 1][i]);
                v_bar_t = (dt * 0.5 * dx0i) * (v[j + 1][i] + v[j][i]);
            }
            else if (j == m - 1)
            {

                v_bar_b = (dt * 0.5 * dx0i) * (v[j][i] + v[j - 1][i]);
                v_bar_t = (dt * 0.5 * dx0i) * (v[0][i] + v[j][i]);
            }
            else
            {

                v_bar_b = (dt * 0.5 * dx0i) * (v[j][i] + v[j - 1][i]);
                v_bar_t = (dt * 0.5 * dx0i) * (v[j + 1][i] + v[j][i]);
            }

            if (sc == 1)
            {
                C[j][i] = C[j][i] * (1.0 + v_bar_t - v_bar_b) + F_t[j][i] + F_b[j][i];
            }
            else
            {
                C[j][i] = (C[j][i] + F_t[j][i] + F_b[j][i]) / (1.0 - v_bar_t + v_bar_b);
            }

            if (C[j][i] < 0.0 && C[j][i] > -1e-10)
            { // correcting for machine precision

                C[j][i] = 0.0;
            }
            else if (C[j][i] > 1.0 && C[j][i] < 1.0 + 1e-10)
            {

                C[j][i] = 1.0;
            }

            if (C[j][i] < 0.0 || C[j][i] > 1.0)
            { // terminating program if color function is out of bounds

                if (sc == 1)
                {
                    printf("for t = %i, at i = %i, j = %i, C = %.16f, which is out of bounds. Program needs to be debugged [y_le]\n", t, i, j, C[j][i]);
                    exit(0);
                }
                else
                {
                    printf("for t = %i, at i = %i, j = %i, C = %.16f, which is out of bounds. Program needs to be debugged [y_ei]\n", t, i, j, C[j][i]);
                    exit(0);
                }
            }
        }
    }

    return C;
}

double *error_m_g(double **C, double **C_0, double mass_0, double dx0, int m, int n)
{

    double *Err_a;
    Err_a = malloc(sizeof(double) * 2);

    double mass_g = 0.0;
    double mass_t = 0.0;

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            mass_g = mass_g + fabs(C[j][i] - C_0[j][i]) * dx0 * dx0;
            mass_t = mass_t + C[j][i] * dx0 * dx0;
        }
    }

    //	Err[0] = mass_g;
    Err_a[1] = (mass_0 - mass_t) / mass_0;

    Err_a[0] = mass_g / mass_t;

    return Err_a;
}

double *error_d(double **S, double **C, int n, int m, double dx0, double x0, double y0, double r, int fcase)
{

    int ij;
    double *Err_d, x_val_a[8], y_val_a[8], dif_x1, dif_x2, dif_y1, dif_y2, dif_1, dif_2, dif;

    Err_d = malloc(sizeof(double *) * 2);

    Err_d[0] = 0.0; // initializing value
    Err_d[1] = 0.0;
    for (int j = 1; j < m - 1; j++)
    {
        for (int i = 1; i < n - 1; i++)
        {
            if (C[j][i] > 1e-10 && C[j][i] < 1.0 - 1e-10)
            {

                ij = i + j * m;

                if (fcase == 1)
                { // interface is line

                    x_val_a[2] = dx0 * i; // coordinate value of cell face
                    x_val_a[3] = dx0 * (i + 1);
                    y_val_a[0] = dx0 * j;
                    y_val_a[1] = dx0 * (j + 1);

                    y_val_a[2] = x0 * x_val_a[2] + y0;   // values of anlaytical functions at cell faces
                    y_val_a[3] = x0 * x_val_a[3] + y0;   // x0 = slope, m, from y = mx + b
                    x_val_a[0] = (y_val_a[0] - y0) / x0; // y0 = y intercept, b, from y = mx + b
                    x_val_a[1] = (y_val_a[1] - y0) / x0;

                    for (int k = 0; k < 4; k++)
                    {
                        if (k < 2)
                        {
                            if (x_val_a[k] > dx0 * i && x_val_a[k] < dx0 * (i + 1))
                            { // checking if analytical y value is within cell

                                dif_x1 = x_val_a[k] - dx0 * (i + S[0][ij]); // diference between analytical coordinates and 1st PLIC point
                                dif_y1 = y_val_a[k] - dx0 * (j + S[2][ij]);

                                dif_x2 = x_val_a[k] - dx0 * (i + S[1][ij]); // difference between analytical cooridnates and 2nd PLIC point
                                dif_y2 = y_val_a[k] - dx0 * (j + S[3][ij]);

                                dif_1 = sqrt(dif_x1 * dif_x1 + dif_y1 * dif_y1); // distance vector between analytical point and PLIC point
                                dif_2 = sqrt(dif_x2 * dif_x2 + dif_y2 * dif_y2);

                                if (dif_1 < dif_2)
                                {
                                    dif = dif_1;
                                }
                                else
                                {
                                    dif = dif_2;
                                }

                                Err_d[0] = Err_d[0] + dif;
                                if (i == 0 || i == n - 1 || j == 0 || j == m - 1)
                                {
                                    Err_d[1] = Err_d[1] + dif;
                                }
                            }
                        }
                        else
                        {
                            if (y_val_a[k] > dx0 * j && y_val_a[k] < dx0 * (j + 1))
                            {

                                dif_x1 = x_val_a[k] - dx0 * (i + S[0][ij]); // diference between analytical coordinates and 1st PLIC point
                                dif_y1 = y_val_a[k] - dx0 * (j + S[2][ij]);

                                dif_x2 = x_val_a[k] - dx0 * (i + S[1][ij]); // difference between analytical cooridnates and 2nd PLIC point
                                dif_y2 = y_val_a[k] - dx0 * (j + S[3][ij]);

                                dif_1 = sqrt(dif_x1 * dif_x1 + dif_y1 * dif_y1); // distance vector between analytical point and PLIC point
                                dif_2 = sqrt(dif_x2 * dif_x2 + dif_y2 * dif_y2);

                                if (dif_1 < dif_2)
                                {
                                    dif = dif_1;
                                }
                                else
                                {
                                    dif = dif_2;
                                }

                                Err_d[0] = Err_d[0] + dif;
                                if (i == 0 || i == n - 1 || j == 0 || j == m - 1)
                                {
                                    Err_d[1] = Err_d[1] + dif;
                                }
                            }
                        }
                    }
                }
                else if (fcase == 2)
                { // interface is circle

                    x_val_a[4] = dx0 * i; // coordinate value of cell face
                    x_val_a[5] = dx0 * i;
                    x_val_a[6] = dx0 * (i + 1);
                    x_val_a[7] = dx0 * (i + 1);
                    y_val_a[0] = dx0 * j;
                    y_val_a[1] = dx0 * j;
                    y_val_a[2] = dx0 * (j + 1);
                    y_val_a[3] = dx0 * (j + 1);

                    y_val_a[4] = y0 + sqrt(r * r - (x_val_a[4] - x0) * (x_val_a[4] - x0)); // values of analytical functions at cell faces
                    y_val_a[5] = y0 - sqrt(r * r - (x_val_a[5] - x0) * (x_val_a[5] - x0));

                    y_val_a[6] = y0 + sqrt(r * r - (x_val_a[6] - x0) * (x_val_a[6] - x0));
                    y_val_a[7] = y0 - sqrt(r * r - (x_val_a[7] - x0) * (x_val_a[7] - x0));

                    x_val_a[0] = x0 + sqrt(r * r - (y_val_a[0] - y0) * (y_val_a[0] - y0));
                    x_val_a[1] = x0 - sqrt(r * r - (y_val_a[1] - y0) * (y_val_a[1] - y0));

                    x_val_a[2] = x0 + sqrt(r * r - (y_val_a[2] - y0) * (y_val_a[2] - y0));
                    x_val_a[3] = x0 - sqrt(r * r - (y_val_a[3] - y0) * (y_val_a[3] - y0));

                    for (int k = 0; k < 8; k++)
                    {
                        if (k < 4)
                        {
                            if (x_val_a[k] > dx0 * i && x_val_a[k] < dx0 * (i + 1))
                            { // checking if analytical y value is within cell

                                dif_x1 = x_val_a[k] - dx0 * (i + S[0][ij]); // diference between analytical coordinates and 1st PLIC point
                                dif_y1 = y_val_a[k] - dx0 * (j + S[2][ij]);

                                dif_x2 = x_val_a[k] - dx0 * (i + S[1][ij]); // difference between analytical cooridnates and 2nd PLIC point
                                dif_y2 = y_val_a[k] - dx0 * (j + S[3][ij]);

                                dif_1 = sqrt(dif_x1 * dif_x1 + dif_y1 * dif_y1); // distance vector between analytical point and PLIC point
                                dif_2 = sqrt(dif_x2 * dif_x2 + dif_y2 * dif_y2);

                                if (dif_1 < dif_2)
                                {
                                    dif = dif_1;
                                }
                                else
                                {
                                    dif = dif_2;
                                }

                                Err_d[0] = Err_d[0] + dif;
                                if (i == 0 || i == n - 1 || j == 0 || j == m - 1)
                                {
                                    Err_d[1] = Err_d[1] + dif;
                                }
                            }
                        }
                        else
                        {
                            if (y_val_a[k] > dx0 * j && y_val_a[k] < dx0 * (j + 1))
                            {

                                dif_x1 = x_val_a[k] - dx0 * (i + S[0][ij]); // diference between analytical coordinates and 1st PLIC point
                                dif_y1 = y_val_a[k] - dx0 * (j + S[2][ij]);

                                dif_x2 = x_val_a[k] - dx0 * (i + S[1][ij]); // difference between analytical cooridnates and 2nd PLIC point
                                dif_y2 = y_val_a[k] - dx0 * (j + S[3][ij]);

                                dif_1 = sqrt(dif_x1 * dif_x1 + dif_y1 * dif_y1); // distance vector between analytical point and PLIC point
                                dif_2 = sqrt(dif_x2 * dif_x2 + dif_y2 * dif_y2);

                                if (dif_1 < dif_2)
                                {
                                    dif = dif_1;
                                }
                                else
                                {
                                    dif = dif_2;
                                }

                                Err_d[0] = Err_d[0] + dif;
                                if (i == 0 || i == n - 1 || j == 0 || j == m - 1)
                                {
                                    Err_d[1] = Err_d[1] + dif;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return Err_d;
}

double *error_ic3(double **S, double **C, int n, int m, double dx0, double Lx, double Ly, double r, double x0, double y0, int fcase)
{

    int ij, nx_L;
    double *Err_ic, dx_L, x1[m][n], y1[m][n], n_x[m][n], n_y[m][n], a_x[m][n], a_y[m][n], error_L[m * n], m0, b0;

    Err_ic = malloc(sizeof(double *) * 2);
    dx_L = 1.0 / pow(2.0, 10.0);

    nx_L = dx0 / dx_L;

    double x_L[nx_L + 1], y_L_p[nx_L + 1], y_L_a1[nx_L + 1], y_L_a2[nx_L + 1], dif[nx_L + 1];

    Err_ic[0] = 0.0; // initializing value
    Err_ic[1] = 0.0;

    for (int j = 1; j < m - 1; j++)
    {
        for (int i = 1; i < n - 1; i++)
        {

            ij = i + j * m;

            if (C[j][i] > 1e-10 && C[j][i] < 1.0 - 1e-10)
            {

                x1[j][i] = S[0][ij]; // extracting values from input arguments
                y1[j][i] = S[2][ij];
                n_x[j][i] = S[4][ij];
                n_y[j][i] = S[5][ij];

                a_x[j][i] = dx0 * (i + x1[j][i]); // converting line segment point to global coordinates
                a_y[j][i] = dx0 * (j + y1[j][i]);

                m0 = -n_x[j][i] / n_y[j][i]; // solving for y = mx + b constants
                b0 = a_y[j][i] - m0 * a_x[j][i];

                for (int k = 0; k < nx_L + 1; k++)
                {

                    x_L[k] = dx0 * i + dx_L * k; // local discretization

                    y_L_p[k] = m0 * x_L[k] + b0; // PLIC

                    if (fcase == 1)
                    { // if interface is line
                        y_L_a1[k] = x0 * x_L[k] + y0;
                    }
                    else if (fcase == 2)
                    { // if interface is circle
                        if (x_L[k] < x0 - r + 1.0e-10 || x_L[k] > x0 + r - 1.0e-10)
                        {
                            y_L_a1[k] = 0.0;
                            y_L_a2[k] = 0.0;
                        }
                        else
                        {
                            y_L_a1[k] = y0 + sqrt(r * r - pow((x_L[k] - x0), 2.0)); // top half of circle
                            y_L_a2[k] = y0 - sqrt(r * r - pow((x_L[k] - x0), 2.0)); // bottom half of circle
                        }
                    }

                    if (y_L_p[k] > dx0 * (j + 1))
                    {

                        y_L_p[k] = dx0 * (j + 1); // boundary condition if PLIC is above cell
                    }
                    else if (y_L_p[k] < dx0 * j)
                    {

                        y_L_p[k] = dx0 * j; // boundary condition if PLIC is below cell
                    }

                    if (y_L_a1[k] > dx0 * (j + 1))
                    {

                        y_L_a1[k] = dx0 * (j + 1); // boundary condition if function is above cell
                    }
                    else if (y_L_a1[k] < dx0 * j)
                    {

                        y_L_a1[k] = dx0 * j; // boundary condition if function is below cell
                    }

                    if (fcase == 2)
                    {
                        if (y_L_a2[k] > dx0 * (j + 1))
                        {

                            y_L_a2[k] = dx0 * (j + 1); // boundary condition if circle is above cell
                        }
                        else if (y_L_a2[k] < dx0 * j)
                        {

                            y_L_a2[k] = dx0 * j; // boundary condition if circle is below cell
                        }
                    }

                    if (fcase == 1)
                    {
                        dif[k] = fabs(y_L_a1[k] - y_L_p[k]); // difference between PLIC and analyitcal line
                    }
                    else if (fcase == 2)
                    {
                        if (y0 > dx0 * (j + 1) || y0 < dx0 * j) // if horizontal edge of circle outside of cell
                        {                                       // because the equation of a circle is in y = f(x) form, two functions will be needed to construct the circle in this cell
                            if (y_L_p[k] > y0)
                            {
                                dif[k] = fabs(fabs(y_L_a1[k] - y_L_a2[k]) - fabs(y_L_p[k] - dx0 * j)); // difference between PLIC and Circle
                            }
                            else
                            {
                                dif[k] = fabs(fabs(y_L_a1[k] - y_L_a2[k]) - fabs(y_L_p[k] - dx0 * (j + 1))); // difference between PLIC and Circle
                            }
                        }
                        else
                        {
                            if (S[5][ij] > 0.0)
                            {
                                dif[k] = fabs(fabs(y_L_a1[k] - y_L_a2[k]) - fabs(y_L_p[k] - dx0 * (j + 1))); // difference between PLIC and Circle
                            }
                            else
                            {
                                dif[k] = fabs(fabs(y_L_a1[k] - y_L_a2[k]) - fabs(y_L_p[k] - dx0 * j)); // difference between PLIC and Circle
                            }
                        }
                    }
                }

                error_L[ij] = (0.5 * dx_L) * dif[0]; // first term of trapizoid integral

                for (int k = 1; k < nx_L; k++)
                {
                    error_L[ij] = error_L[ij] + dx_L * dif[k];
                }

                error_L[ij] = error_L[ij] + (0.5 * dx_L) * dif[nx_L];

                Err_ic[0] = Err_ic[0] + error_L[ij];

                if (i == 0 || i == n - 1 || j == 0 || j == m - 1)
                {
                    Err_ic[1] = Err_ic[1] + error_L[ij];
                }
            }
        }
    }

    return Err_ic;
}

int main()
{
    // Volume of Fluids Method
    // Drake Jones
    // M.S. in Mechanical Engineering Thesis
    // San Diego State University
    // Fall 2023

    clock_t t0_t;
    t0_t = clock();

    // initialzing variables
    double r, x0, y0, dx0, Lx, dt, t_f, CFL, u0, v0, w, *ipv_x, *ipv_y, t_ic, t_ad, t_t;
    double vertp0[2 * nv], vf, mass_0, mass_a, A, test;
    int n, m, nt, ccase, ICM, fcase, acase, ntv1, sc1, sc2, ij, glbvt[4], ipv0[nv], grd[6], nm_tm[6], numtm0[6], nt_8;

    int nc = 50;       // local number of divisions along each coordinate
    double tol = 10.0; //	posative tolerance
    int ntv2 = 4;      // number fo total vertices in a single cell
    int ntp0 = 4;      // last local vertex (following VoF tools test format)
    int t = 0;

    double **C, **C_0, **S, **u, **v;

    ICM = 3; // interface construction method: 1 = centered columns method 2 = Youngs method
             // 3 = Mixed Youngs-centered	4 = ELVIRA

    ccase = 5; // 1 = interface reconstruction; 2 = Linear Translation; 3 = Rotational Translation;  5 = reverse vortex flow  

    acase = 4; // Advection Scheme: 1 = EI-EI;     2 = LE-Le;      3 = LE-EI;      4 = EI-LE

    fcase = 2; // 1 = line;	2 = cirlce;

    n = 8; // numer of cells across the diameter; Must obey 2^N
    m = n;
    nt = 641; // number of time steps

    if (ccase == 1)
    {
        printf("Test case: Interface Reconstruction\n");
        if (fcase == 1)
        {
            x0 = -1.5; // SLOPE; wanted to use same varibes or x0, y0 and m,b
            y0 = 5.1;  // Y INTERCEPT; wanted to use same variables as x0,y0 and m,b
            r = 1.0;   // need this for this to run
        }
        else if (fcase == 2)
        {
            r = 1.0;   // radius of cirlce
            x0 = 1.51; // x coordinate of cirlce orgin
            y0 = 1.6;  // y coordinate of cirlce orgin
        }

        Lx = 4.0; // domain in x and y direction
        dx0 = Lx / n;
        nt = 1;
    }
    else if (ccase == 2)
    { // Linear translation

        printf("Test Case: Linear Translation\n");

        r = 0.25;  // radius of cirlce
        x0 = 0.37; // x coordinate of cirlce orgin
        y0 = 0.38; // y coordinate of cirlce orgin

        u0 = 1.0; // initial velocity in x direction
        v0 = 1.0; // initial velocity in y direction

        Lx = 1.0;     // domain in x and y direction
        dx0 = Lx / n; // cell size in x and y direction

        t_f = fabs(2 * sqrt(2 * (0.62 - 0.38) * (0.62 - 0.38)) / sqrt(u0 * u0 + v0 * v0));
        dt = t_f / nt;
        // printf("d = %f  tf = %f  dt = %f\n",sqrt(2*(2.49-1.51)*(2.49-1.51)),t_f,dt);
        CFL = ((fabs(u0) + fabs(v0)) * dt) / dx0;
    }
    else if (ccase == 3)
    { // rotational translation

        printf("Test Case: Rotational Translaiton\n");

        r = 0.15;  // radius fo cirlce
        x0 = 0.26; // x coordinate of cirlce orgin
        y0 = 0.48; // y coordinate of cirlce orgin

        Lx = 1.0;     // domain in x and y direction
        dx0 = Lx / n; // cell size in x and y direction

        w = M_PI / 2.0; // Angluar Speed

        t_f = 4.0;
        dt = t_f / nt;

        CFL = (w * dt) / dx0;
    }
    else if (ccase == 5)
    { // Reverse Vortex flow

        printf("Test Case: Reverse Vortex Flow\n");

        w = M_PI; // Angular frequency [rad/s]

        A = 2.0; // amplitude

        r = 0.15;  // radius of circle
        x0 = 0.5;  // x coordinate of cirlce orgin
        y0 = 0.75; // y coordinate of cirlce orgin

        Lx = 1.0;
        dx0 = Lx / n;

        t_f = 8.0;
        dt = t_f / (nt - 1);

        CFL = (t_f / (dx0 * nt)) * (0.755 + 0.5295);
    }

    printf("N = %i\t nt = %i\t dt = %f\t CFL = %f\n", n, nt, dt, CFL);

    u = malloc(sizeof(double *) * m);
    v = malloc(sizeof(double *) * m);

    for (int j = 0; j < m; j++)
    {
        u[j] = malloc(sizeof(double *) * n);
        v[j] = malloc(sizeof(double *) * n);
    }

    ntv1 = (m + 1) * (n + 1); // number of total  vertices
                              // printf("Testing01\n");

    ipv_x = malloc(sizeof(double *) * ntv1);
    ipv_y = malloc(sizeof(double *) * ntv1);

    mass_0 = 0.0; // initializing initial mass
    for (int j = 0; j < m + 1; j++)
    {
        for (int i = 0; i < n + 1; i++)
        {
            ij = i + j * (n + 1); // turing 2d grid into 1d array

            ipv_x[ij] = dx0 * i; // x cooridnate of each vertice
            ipv_y[ij] = dx0 * j; // y coordinate of each vertice
        }
    }

    glbvt[0] = -1; // intialzing so indexing loop works
    glbvt[1] = 0;
    glbvt[2] = n + 1;
    glbvt[3] = n;

    for (int i = 0; i < 4; i++)
    {
        ipv0[i] = i + 1; // initializing ipv for rectangular grid
    }

    C = malloc(sizeof(double *) * m);
    C_0 = malloc(sizeof(double *) * m);

    for (int j = 0; j < m; j++)
    {
        C[j] = malloc(sizeof(double *) * n);
        C_0[j] = malloc(sizeof(double *) * n);
    }

    S = malloc(sizeof(double *) * 7);
    for (int j = 0; j < 7; j++)
    {
        S[j] = malloc(sizeof(double *) * m * n);
    }

    int cnt = 0;
    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            if (i == 0)
            {
                glbvt[0] = glbvt[0] + 2; // global vertices for each cell
                glbvt[1] = glbvt[1] + 2;
                glbvt[2] = glbvt[2] + 2;
                glbvt[3] = glbvt[3] + 2;
            }
            else
            {
                glbvt[0] = glbvt[0] + 1;
                glbvt[1] = glbvt[1] + 1;
                glbvt[2] = glbvt[2] + 1;
                glbvt[3] = glbvt[3] + 1;
            }
            vertp0[0] = ipv_x[glbvt[0] - 1]; // global vertex coordinates
            vertp0[1] = ipv_x[glbvt[1] - 1];
            vertp0[2] = ipv_x[glbvt[2] - 1];
            vertp0[3] = ipv_x[glbvt[3] - 1];
            vertp0[120] = ipv_y[glbvt[0] - 1];
            vertp0[121] = ipv_y[glbvt[1] - 1];
            vertp0[122] = ipv_y[glbvt[2] - 1];
            vertp0[123] = ipv_y[glbvt[3] - 1];

            initf2d_(func2d2_, ipv0, &nc, &ntp0, &ntv2, &tol, vertp0, &vf); // VoF tools initialization function for egg

            C[j][i] = fabs(1.0 - vf);   // color function
            C_0[j][i] = fabs(1.0 - vf); // storing original color function for error

            mass_0 = mass_0 + C[j][i] * dx0 * dx0; // Calculating mass of original color function

        }
    }

    printf("m_0 = %f\n", mass_0);

    double *Err_ic, *Err_d, *Err_a;

    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            C[j][i] = C_0[j][i];
        }
    }

    if (ICM == 1)
    {
        printf("ICM: Centerd Columns Method\n");
    }
    else if (ICM == 2)
    {
        printf("ICM: Youngs' Method\n");
    }
    else if (ICM == 3)
    {
        printf("ICM: Mixed Youngs' Centered\n");
    }
    else if (ICM == 4)
    {
        printf("ICM: ELVIRA\n");
    }

    if (ccase != 1)
    {
        if (acase == 1)
        {
            printf("Advection Scheme: EI-EI\n");
        }
        else if (acase == 2)
        {
            printf("Advetion Scheme: LE-LE\n");
        }
        else if (acase == 3)
        {
            printf("Advection Scheme: LE-EI\n");
        }
        else if (acase == 4)
        {
            printf("Advection Scheme: EI-LE\n");
        }
    }

    clock_t t0_ic;
    clock_t t0_ad; // initializing clock
    t0_ic = clock();
    t0_ad = clock();

    S = PLIR(C, S, dx0, m, n, ICM); // PLIR at t = 0;

    t0_ic = clock() - t0_ic;
    if (ccase == 1)
    {
        t_ic = ((double)t0_ic) / CLOCKS_PER_SEC;
        printf("t_ic = %f\n", t_ic);
    }

    if (ccase == 1)
    {
        Err_ic = error_ic3(S, C, n, m, dx0, Lx, Lx, r, x0, y0, fcase);
           printf("E_ic = %e\n", Err_ic[0]);
        Err_d = error_d(S, C, n, m, dx0, x0, y0, r, fcase);
            printf("E_ic_t = %e E_ic_b = %e  E_d_t = %e  E_d_b = %e\n", Err_ic[0], Err_ic[1], Err_d[0], Err_d[1]);
    }
    test = 0.0;
    if (ccase != 1)
    { // interface reconstruction error

        for (int t = 1; t < nt; t++)
        { // time step

            for (int j = 0; j < m; j++)
            { // test velocity profiles
                for (int i = 0; i < n; i++)
                {
                    if (ccase == 2)
                    {

                        if (t <= nt / 2)
                        {
                            u[j][i] = u0;
                            v[j][i] = u0;
                        }
                        else
                        {
                            u[j][i] = -u0;
                            v[j][i] = -u0;
                        }
                    }
                    else if (ccase == 3)
                    {
                        u[j][i] = -w * (dx0 * (j + 0.5) - 0.5);
                        v[j][i] = w * (dx0 * (i +0.5) - 0.5);
                    }
                    else if (ccase == 5)
                    {
                        u[j][i] = A * sin(w * dx0 * (i+0.5)) * sin(w * dx0 * (i+0.5)) * sin(w * dx0 * (j+0.5)) * cos(w * dx0 * (j+0.5)) * cos((M_PI / 4.0) * t * dt); // reverse vortex flow
                        v[j][i] = -A * sin(w * dx0 * (i+0.5)) * cos(w * dx0 * (i+0.5)) * sin(w * dx0 * (j+0.5)) * sin(w * dx0 * (j+0.5)) * cos((M_PI / 4.0) * t * dt);
                    }

                    if (dt * u[j][i] >= dx0 || dt * v[j][i] >= dx0)
                    { // exiting program if lengthstep is too large
                        printf("ALERT!\n");
                        printf("i = %i\t j = %i\t dx = %f\t u*dt = %f\t vdt = %f\n", i, j, dx0, u[j][i] * dt, v[j][i] * dt);
                        printf("u_i*dt > dx0. Adjust parameters so v_i*dt is less than the grid size.\n");
                        exit(0);
                    }
                }
            }

            if (t % 2 == 0)
            {

                if (acase == 1)
                { // EI-EI

                    sc1 = 2; // Eulerian Implicit
                    sc2 = 2; // Eulerian Implicit
                }
                else if (acase == 2)
                { // LE-LE

                    sc1 = 1; // Lagrangian Explicit
                    sc2 = 1; // Lagrangian Explicit
                }
                else if (acase == 3)
                { // LE-EI

                    sc1 = 1; // Lagrangian Explict
                    sc2 = 2; // Eularian Implicit
                }
                else
                {

                    sc1 = 2; // Eulerian Implicit
                    sc2 = 1; // Lagranian Explicit
                }

                C = advec_x(C, S, dx0, m, n, u, dt, t, sc1);

                S = PLIR(C, S, dx0, m, n, ICM);

                C = advec_y(C, S, dx0, m, n, v, dt, t, sc2);

                S = PLIR(C, S, dx0, m, n, ICM);
            }
            else if (t % 2 == 1)
            {

                if (acase == 1)
                { // EI-EI

                    sc1 = 2; // Eulerian Implicit
                    sc2 = 2; // Eulerian Implicit
                }
                else if (acase == 2)
                { // LE-LE

                    sc1 = 1; // Lagrangian Explicit
                    sc2 = 1; // Lagrangian Explicit
                }
                else if (acase == 3)
                { // LE-EI

                    sc1 = 1; // Lagrangian Explict
                    sc2 = 2; // Eularian Implicit
                }
                else
                {

                    sc1 = 2; // Eulerian Implicit
                    sc2 = 1; // Lagranian Explicit
                }

                C = advec_y(C, S, dx0, m, n, v, dt, t, sc1);

                S = PLIR(C, S, dx0, m, n, ICM);

                C = advec_x(C, S, dx0, m, n, u, dt, t, sc2);

                S = PLIR(C, S, dx0, m, n, ICM);
            }
        }

        t0_ad = clock() - t0_ad;
        if (ccase != 1)
        {
            t_ad = ((double)t0_ad) / CLOCKS_PER_SEC;
            printf("t_ad = %f\n", t_ad);
        }

        // advection mass conservation and geometric error
        Err_a = error_m_g(C, C_0, mass_0, dx0, m, n);
        printf("E_m = %e\t E_g = %e\n", Err_a[1], Err_a[0]);
        Err_ic = error_ic3(S, C, n, m, dx0, Lx, Lx, r, x0, y0, fcase);
        Err_d = error_d(S, C, n, m, dx0, x0, y0, r, fcase);
        printf("E_ic = %e\t E_d = %e\n", Err_ic[0], Err_d[0]);
    }
    printf("----------\n");

    for (int i = 0; i < m; i++)
    {
        free(u[i]);
        free(v[i]);
        free(C[i]);
        free(C_0[i]);
    }

    for (int i = 0; i < 7; i++)
    {
        free(S[i]);
    }

    free(u);
    free(v);
    free(C);
    free(C_0);
    free(S);

    t0_t = clock() - t0_t;

    t_t = ((double)t0_t) / CLOCKS_PER_SEC;
    printf("t_t = %f\n", t_t);

    return 0;
}
