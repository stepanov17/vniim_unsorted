clear;
clc;

// F crit. points

n1 = 100;
n2 = 200;

P0 = 0.95;

dP = 0.5 * (1. - P0);

c1 = 1. / cdff("F", n2 - 1, n1 - 1, 1. - dP, dP);
c2 =      cdff("F", n1 - 1, n2 - 1, 1. - dP, dP);

printf("c1 = %.4f, c2 = %.4f\n", c1, c2)
