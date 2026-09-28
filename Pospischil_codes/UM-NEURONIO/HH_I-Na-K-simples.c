//~ gcc -O3 -o x.x HH_I-Na-K-simples.c -lm; time ./x.x &

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define n_n 20000 // numero de passos
#define NR_END 1
#define FREE_ARG char *
#define N 5             // numero de equacoes
#define NN 1            // quantidade de HH
#define neuron_classe 1 // number of neurons classes

FILE *input;

void derivs(double y[], double df[], double *Iext);
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void nrerror(char error_text[]);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
int delta(int a, int b);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);

int main(void)
{

   int i, t;
   double *x, tempo, h, *a, *b, *c, *df, *y, *Iext;

   y = dvector(1, N * NN + 1);
   df = dvector(1, N * NN + 1);
   x = dvector(1, N * NN + 1);
   a = dvector(1, N * NN + 1);
   b = dvector(1, N * NN + 1);
   c = dvector(1, N * NN + 1);

   Iext = dvector(1, N * NN + 1);

   h = 0.05; // passo de integracao

   input = fopen("HH_N1.dat", "wt");

   // condicoes iniciais
   for (i = 1; i <= N * NN; i = i + 5)
   {
      x[i] = -70.0; //*ran1(& idum);
      x[i + 1] = 0.0;
      x[i + 2] = 0.0;
      x[i + 3] = 0.0;
      x[i + 4] = 0.0;

      Iext[i] = 0.0;
   }

   // loop do tempo
   tempo = 0.0;
   for (t = 1; t <= n_n; t++)
   {
      tempo = tempo + h;

      for (i = 1; i <= N * NN; i++)
         y[i] = x[i];

      // correte externa aplicada entre 1000 ms e 3000 ms
      if (tempo > 300.0 && tempo < 700.00)
         for (i = 1; i <= N * NN; i = i + N) // depolarizing current pulse in uA
            Iext[i] = 0.0007;

      if (tempo > 700.01 && tempo < 3000.02)
         for (i = 1; i <= N * NN; i = i + N)
            Iext[i] = 0.0;

      // ------------ Runge-Kutta --------------------
      derivs(y, df, Iext);
      for (i = 1; i <= N * NN; i++)
      {
         a[i] = h * df[i];
         y[i] = x[i] + a[i] / 2.0;
      }
      derivs(y, df, Iext);
      for (i = 1; i <= N * NN; i++)
      {
         b[i] = h * df[i];
         y[i] = x[i] + b[i] / 2.0;
      }
      derivs(y, df, Iext);
      for (i = 1; i <= N * NN; i++)
      {
         c[i] = h * df[i];
         y[i] = x[i] + c[i];
      }
      derivs(y, df, Iext);
      for (i = 1; i <= N * NN; i++)
         x[i] = x[i] + (a[i] + h * df[i]) / 6.0 + (b[i] + c[i]) / 3.0;
      // ------------

      fprintf(input, "%.2f %.3f %.3f\n", tempo, x[1], x[5]);

   } // fim do loop do tempo

   free_dvector(y, 1, N * NN + 1);
   free_dvector(df, 1, N * NN + 1);
   free_dvector(x, 1, N * NN + 1);
   free_dvector(a, 1, N * NN + 1);
   free_dvector(b, 1, N * NN + 1);
   free_dvector(c, 1, N * NN + 1);

   free_dvector(Iext, 1, N * NN + 1);

   return 0;

   fclose(input);
}

////////// equacoes
void derivs(double y[], double df[], double *Iext)
{
   int i, ii, neurontype;
   double C_m, E_Na, E_K;
   double Area[neuron_classe + 1], V_T[neuron_classe + 1], E_leak[neuron_classe + 1], g_leak[neuron_classe + 1], g_Na[neuron_classe + 1], g_Kd[neuron_classe + 1], g_M[neuron_classe + 1], tau_max[neuron_classe + 1];
   double V, m, h, n, p, I_Na, alpha_m, beta_m, alpha_h, beta_h, dm, dh, I_Kd, alpha_n, beta_n, dn, I_M, p_inf, tau_p, dp;

   C_m = 1.0;   // uF/cm^2
   E_Na = 50.0; // mV
   E_K = -90.0; // mV

   // Regular-spiking pyramidal neuron - Fig. 1
   neurontype = 1;
   Area[neurontype] = 0.0096 * 0.0096 * 3.14159; // cm^2
   V_T[neurontype] = -55.0;                      // mV
   E_leak[neurontype] = -70.0;                   // mV
   g_leak[neurontype] = 0.1;                     // mS/cm^2
   g_Na[neurontype] = 50.0;                      // mS/cm^2
   g_Kd[neurontype] = 5.0;                       // mS/cm^2
   g_M[neurontype] = 0.07;                       // mS/cm^2
   tau_max[neurontype] = 1000.0;                 // ms

   //////////////////////////EQUAÇOES ACOPLADAS////////////////////

   for (i = 1; i <= N * NN; i = i + 5)
   {
      ii = 1;

      V = y[i];
      m = y[i + 1];
      h = y[i + 2];
      n = y[i + 3];
      p = y[i + 4];

      // sodium current responsible for action potentials
      I_Na = g_Na[ii] * pow(m, 3) * h * (V - E_Na);
      alpha_m = (-0.32 * (V - V_T[ii] - 13.0)) / (exp(-(V - V_T[ii] - 13.0) / 4.0) - 1.0);
      beta_m = (0.28 * (V - V_T[ii] - 40.0)) / (exp((V - V_T[ii] - 40.0) / 5.0) - 1.0);
      alpha_h = 0.128 * exp(-(V - V_T[ii] - 17.0) / 18.0);
      beta_h = 4.0 / (1.0 + exp(-(V - V_T[ii] - 40.0) / 5.0));

      dm = alpha_m * (1 - m) - beta_m * m;
      dh = alpha_h * (1 - h) - beta_h * h;

      // potassium current responsible for action potentials
      I_Kd = g_Kd[ii] * (pow(n, 4) * (V - E_K));
      alpha_n = (-0.032 * (V - V_T[ii] - 15.0)) / (exp(-(V - V_T[ii] - 15.0) / 5.0) - 1.0);
      beta_n = 0.5 * exp(-(V - V_T[ii] - 10.0) / 40.0);

      dn = alpha_n * (1 - n) - beta_n * n;

      // slow voltage-dependent potassium current responsible for spike-frequency adaptation
      I_M = g_M[ii] * p * (V - E_K);
      p_inf = 1.0 / (1.0 + exp(-(V + 35.0) / 10.0));
      tau_p = tau_max[ii] / (3.3 * exp((V + 35.0) / 20.0) + exp(-(V + 35.0) / 20.0));

      dp = (p_inf - p) / tau_p;

      df[i] = (1.0 / C_m) * (-g_leak[ii] * (V - E_leak[ii]) - I_Na - I_Kd - I_M + (Iext[i] / Area[ii]));
      df[i + 1] = dm;
      df[i + 2] = dh;
      df[i + 3] = dn;
      df[i + 4] = dp;
   }
}
/////////////////////////////////////////

double *dvector(long nl, long nh)
{
   double *v;

   v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
   if (!v)
      nrerror("allocation failure in dvector()");
   return v - nl + NR_END;
}

void free_dvector(double *v, long nl, long nh)
{
   free((FREE_ARG)(v + nl - NR_END));
}

void nrerror(char error_text[])
{
   fprintf(stderr, "Numerical Recipes run-time error...\n");
   fprintf(stderr, "%s\n", error_text);
   fprintf(stderr, "...now exiting to system...\n");
   exit(1);
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
   double **m;

   m = (double **)malloc((size_t)((nrow + NR_END) * sizeof(double *)));
   if (!m)
      nrerror("allocation failure 1 in matrix()");
   m += NR_END;
   m -= nrl;

   m[nrl] = (double *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
   if (!m[nrl])
      nrerror("allocation failure 2 in matrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for (i = nrl + 1; i <= nrh; i++)
      m[i] = m[i - 1] + ncol;

   return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
   free((FREE_ARG)(m[nrl] + ncl - NR_END));
   free((FREE_ARG)(m + nrl - NR_END));
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
   int **m;

   m = (int **)malloc((size_t)((nrow + NR_END) * sizeof(int *)));
   if (!m)
      nrerror("allocation failure 1 in matrix()");
   m += NR_END;
   m -= nrl;

   m[nrl] = (int *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));
   if (!m[nrl])
      nrerror("allocation failure 2 in matrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for (i = nrl + 1; i <= nrh; i++)
      m[i] = m[i - 1] + ncol;

   return m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
   free((FREE_ARG)(m[nrl] + ncl - NR_END));
   free((FREE_ARG)(m + nrl - NR_END));
}
