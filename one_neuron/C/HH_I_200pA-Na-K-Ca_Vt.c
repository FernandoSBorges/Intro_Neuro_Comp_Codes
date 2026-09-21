//~ Biol Cybern (2008) 99:427–441
//~ DOI 10.1007/s00422-008-0263-8
//~ Minimal Hodgkin–Huxley type models for different classes
//~ of cortical and thalamic neurons
//~ Martin Pospischil · Maria Toledo-Rodriguez ·
//~ Cyril Monier · Zuzanna Piwkowska · Thierry Bal ·
//~ Yves Frégnac · Henry Markram · Alain Destexhe

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define PI acos(-1.0)

#define n_n 300000
#define NR_END 1
#define FREE_ARG char *
#define N 9             // eq number
#define NN 1            // HH number
#define neuron_classe 4 // number of neurons classes

#define stim_delay 500.0 // in ms
#define stim_dur 2000.0

float ran1(long *idum);
float gasdev(long *idum);

FILE *input;
FILE *input2;

void derivs(double y[], double df[], double *Iext, double gM_bar, double gL_bar, double gT_bar);
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

  int i, j, t, spikeID[NN + 2], ii;
  double *x, time, step, *a, *b, *c, *df, *y, *Iext, I_uA, Eca, pico[NN + 2], spikeFreq[NN + 2][11], spikeTime[NN + 2][11], tud[NN + 2], time10;
  double gM_bar, gL_bar, gT_bar;
  char label[100];
  long idum = -254321547251;

  y = dvector(1, N * NN + 1);
  df = dvector(1, N * NN + 1);
  x = dvector(1, N * NN + 1);
  a = dvector(1, N * NN + 1);
  b = dvector(1, N * NN + 1);
  c = dvector(1, N * NN + 1);

  Iext = dvector(1, NN + 1);

  step = 0.01; // RK4 integration step

  for (gM_bar = 0.03; gM_bar <= 0.07001; gM_bar = gM_bar + 0.02) // slow voltage-dependent potassium current
    for (gL_bar = 0.0; gL_bar <= 0.4000001; gL_bar = gL_bar + 0.2) // low-threshold Calcium current
      for (gT_bar = 0.0; gT_bar <= 0.8000001; gT_bar = gT_bar + 0.4) // high-threshold	Calcium current
      {      
        
        printf(" gM_bar = %.3f      gL_bar = %.3f      gT_bar = %.3f \n", gM_bar, gL_bar, gT_bar); 

        sprintf(label,"HH_V_RS_I200pA_gM%.3fe-3_gL%.3fe-3_gT%.3fe-3.dat",gM_bar,gL_bar,gT_bar);
        input=fopen(label,"wt");
        sprintf(label,"HH_Freq_RS_I200pA_gM%.3fe-3_gL%.3fe-3_gT%.3fe-3.dat",gM_bar,gL_bar,gT_bar);
        input2=fopen(label,"wt");

        ii = 0;

        for (I_uA = 0.000200; I_uA <= 0.0002001; I_uA = I_uA + 0.000005) // depolarizing current pulse in uA
        {
          // initial conditions
          for (i = 1; i <= NN; i++)
          {
            x[1 + (i - 1) * N] = -84.0; //*ran1(& idum);
            x[2 + (i - 1) * N] = 0.0;
            x[3 + (i - 1) * N] = 0.0;
            x[4 + (i - 1) * N] = 0.0;
            x[5 + (i - 1) * N] = 0.0;
            x[6 + (i - 1) * N] = 0.0;
            x[7 + (i - 1) * N] = 0.0;
            x[8 + (i - 1) * N] = 0.00024;
            x[9 + (i - 1) * N] = 0.0;

            Iext[i] = 0.0; // depolarizing current pulse in uA
            pico[i] = 0.0;

            spikeID[i] = 0;
            for (j = 1; j <= 10; j++)
              spikeFreq[i][j] = 0;
              spikeTime[i][j] = 0;
            tud[i] = 0;
          }

          time = 0.0;
          time10 = 0;

          ////////////////////////////////////////////////////////////loop do time
          for (t = 1; t <= n_n; t++)
          {
            time = time + step;

            for (i = 1; i <= N * NN; i++)
              y[i] = x[i];

            // depolarizing current pulse in uA
            if (time > stim_delay && time < stim_delay + step)
              for (i = 1; i <= NN; i++)
                Iext[i] = I_uA;

            if (time > stim_delay + stim_dur && time < stim_delay + stim_dur + step)
              for (i = 1; i <= NN; i++)
                Iext[i] = 0.0;

            // ------------ Runge-Kutta --------------------
            derivs(y, df, Iext, gM_bar, gL_bar, gT_bar);
            for (i = 1; i <= N * NN; i++)
            {
              a[i] = step * df[i];
              y[i] = x[i] + a[i] / 2.0;
            }
            derivs(y, df, Iext, gM_bar, gL_bar, gT_bar);
            for (i = 1; i <= N * NN; i++)
            {
              b[i] = step * df[i];
              y[i] = x[i] + b[i] / 2.0;
            }
            derivs(y, df, Iext, gM_bar, gL_bar, gT_bar);
            for (i = 1; i <= N * NN; i++)
            {
              c[i] = step * df[i];
              y[i] = x[i] + c[i];
            }
            derivs(y, df, Iext, gM_bar, gL_bar, gT_bar);
            for (i = 1; i <= N * NN; i++)
              x[i] = x[i] + (a[i] + step * df[i]) / 6.0 + (b[i] + c[i]) / 3.0;
            // -------------------------------------------

            time10 = time10 + step;

            if (ii == 0 || ii == 10 || ii == 20 || ii == 30)
              if (time10 > 0.01 && time>300 && time<2800)
              {
                time10 = 0;
                fprintf(input, "%.2f", time - 0.0); // nA
                for (i = 1; i <= NN; i++)
                  fprintf(input, " %.2f", x[1 + (i - 1) * N] - ((ii/10.0) - 0) * 140.0);
                fprintf(input, "\n");
              }

            for (i = 1; i <= NN; i++)
            {
              if (x[1 + (i - 1) * N] < -20.0)
                pico[i] = 0.0;

              if (x[1 + (i - 1) * N] > pico[i]) // save the spike when potential V_i reaches 0 mV
              {
                pico[i] = 1000;
                spikeID[i] = spikeID[i] + 1;

                if (spikeID[i] <= 10)
                  spikeTime[i][spikeID[i]] = time;

                if (tud[i] > 0.0 && spikeID[i] <= 11)
                  spikeFreq[i][spikeID[i] - 1] = 1000.0 / (time - tud[i]);

                tud[i] = time;
              }
            }

          } // time loop end          

          if (ii == 0 || ii == 10 || ii == 20 || ii == 30)
          {
            fprintf(input, "\n");
            printf("  %f\n", 1000.0 * I_uA);        // nA
          }

          ii = ii + 1;
          
          printf("%f\n", 1000.0 * I_uA);        // nA
          fprintf(input2, "%.2f", 1000.0 * 1000.0 * I_uA); // pA
          for (i = 1; i <= NN; i++)
            for (j = 1; j <= 10; j++)
              fprintf(input2, " %.2f", spikeFreq[i][j]);    
          for (i = 1; i <= NN; i++)    
            fprintf(input2, " %.2f %.2f %d %.2f %.2f",spikeID[i]/2.0,spikeID[i]*1000.0/(tud[i]-stim_delay),spikeID[i],spikeTime[i][1]-stim_delay,spikeTime[i][10]-stim_delay);
          fprintf(input2, "\n");

        } // I_uA loop end

        fclose(input);
        fclose(input2);

      } // g_bars loop end

  free_dvector(y, 1, N * NN + 1);
  free_dvector(df, 1, N * NN + 1);
  free_dvector(x, 1, N * NN + 1);
  free_dvector(a, 1, N * NN + 1);
  free_dvector(b, 1, N * NN + 1);
  free_dvector(c, 1, N * NN + 1);
  free_dvector(Iext, 1, N * NN + 1);

  return 0;
}

////////// equations
void derivs(double y[], double df[], double *Iext, double gM_bar, double gL_bar, double gT_bar)
{
  int i, ii, neurontype;
  double C_m, E_Na, E_K, E_Ca;
  double V_T[neuron_classe + 1], E_leak[neuron_classe + 1], g_leak[neuron_classe + 1], g_Na[neuron_classe + 1], g_Kd[neuron_classe + 1], g_M[neuron_classe + 1], tau_max[neuron_classe + 1];
  double V, m, h, n, p, I_Na, alpha_m, beta_m, alpha_h, beta_h, dm, dh, I_Kd, alpha_n, beta_n, dn, I_M, p_inf, tau_p, dp;
  double Area[neuron_classe + 1];
  double V_x[neuron_classe + 1], g_L[neuron_classe + 1], g_T[neuron_classe + 1];
  double q, r, I_L, alpha_q, beta_q, alpha_r, beta_r, dq, dr, s_inf, I_T, u, u_inf, tau_u, du;
  double cai, dcai, drivechannel, depth, cai_inf, tau_r;

  C_m = 1.0;    // uF/cm^2
  E_Na = 50.0;  // mV
  E_K = -100.0; //-90.0 mV no paper // -100.0 mV no modelDB
  depth = 1.0;
  cai_inf = 0.00024;
  tau_r = 5.0;

  // ------------------------------------------------------------- //
  // Regular-spiking pyramidal neuron
  neurontype = 1;
  Area[neurontype] = 0.0096 * 0.0096 * PI; // cm^2
  V_T[neurontype] = -55.0;                      // mV
  E_leak[neurontype] = -85.0;                   // mV
  g_leak[neurontype] = 0.01;                     // mS/cm^2
  g_Na[neurontype] = 50.0;                      // mS/cm^2
  g_Kd[neurontype] = 5.0;                       // mS/cm^2
  g_M[neurontype] = gM_bar;                       // mS/cm^2
  tau_max[neurontype] = 1000.0;                 // ms
  V_x[neurontype] = 2.0;                        // mV
  g_L[neurontype] = gL_bar;                        // mS/cm^2
  g_T[neurontype] = gT_bar;                        // mS/cm^2

  //////////////////////////EQUAÇOES ACOPLADAS////////////////////
  for (i = 1; i <= NN; i++)
  {
    ii = 1; // neuron type

    V = y[1 + (i - 1) * N];
    m = y[2 + (i - 1) * N];
    h = y[3 + (i - 1) * N];
    n = y[4 + (i - 1) * N];
    p = y[5 + (i - 1) * N];
    q = y[6 + (i - 1) * N];
    r = y[7 + (i - 1) * N];
    cai = y[8 + (i - 1) * N]; // concentration of the Calcium (Ca2+) in the intracellular fluid in mM
    u = y[9 + (i - 1) * N];

    // Calcium reversal potential
    E_Ca = 13.320242889960909 * log(2.0 / cai); // mV

    // sodium current responsible for action potentials
    I_Na = g_Na[ii] * pow(m, 3) * h * (V - E_Na);
    alpha_m = (-0.32 * (V - V_T[ii] - 13.0)) / (exp(-(V - V_T[ii] - 13.0) / 4.0) - 1.0);
    beta_m = (0.28 * (V - V_T[ii] - 40.0)) / (exp((V - V_T[ii] - 40.0) / 5.0) - 1.0);
    alpha_h = 0.128 * exp(-(V - V_T[ii] - 17.0) / 18.0);
    beta_h = 4.0 / (1.0 + exp(-(V - V_T[ii] - 40.0) / 5.0));

    dm = alpha_m * (1.0 - m) - beta_m * m;
    dh = alpha_h * (1.0 - h) - beta_h * h;

    // potassium current responsible for action potentials
    I_Kd = g_Kd[ii] * (pow(n, 4) * (V - E_K));
    alpha_n = (-0.032 * (V - V_T[ii] - 15.0)) / (exp(-(V - V_T[ii] - 15.0) / 5.0) - 1.0);
    beta_n = 0.5 * exp(-(V - V_T[ii] - 10.0) / 40.0);

    dn = alpha_n * (1.0 - n) - beta_n * n;

    // slow voltage-dependent potassium current responsible for spike-frequency adaptation
    I_M = g_M[ii] * p * (V - E_K);
    p_inf = 1.0 / (1.0 + exp(-(V + 35.0) / 10.0));
    tau_p = tau_max[ii] / (3.3 * exp((V + 35.0) / 20.0) + exp(-(V + 35.0) / 20.0));

    dp = (p_inf - p) / tau_p;

    // low-threshold Calcium currents to generate bursting
    I_L = g_L[ii] * q * q * r * (V - E_Ca); // not causing [Ca2+] influx

    alpha_q = (0.055 * (-27.0 - V)) / (exp((-27.0 - V) / 3.8) - 1.0);
    beta_q = 0.94 * exp((-75.0 - V) / 17.0);
    alpha_r = 0.000457 * exp((-13.0 - V) / 50.0);
    beta_r = 0.0065 / (exp((-15.0 - V) / 28.0) + 1.0);

    dq = alpha_q * (1.0 - q) - beta_q * q;
    dr = alpha_r * (1.0 - r) - beta_r * r;

    // high-threshold	Calcium currents to generate bursting
    s_inf = 1.0 / (1.0 + exp(-(V + V_x[ii] + 57.0) / 6.2));

    I_T = g_T[ii] * pow(s_inf, 2) * u * (V - E_Ca);

    u_inf = 1.0 / (1.0 + exp((V + V_x[ii] + 81.0) / 4.0));

    if (V < 0.0) // approximation to "exp()<10^32"
      tau_u = 30.8 / 3.7372 + (211.4 + exp((V + V_x[ii] + 113.2) / 5.0)) / (3.7372 * (1.0 + exp((V + V_x[ii] + 84.0) / 3.2)));
    else
      tau_u = 30.8 / 3.7372 + ((exp(((V + V_x[ii] + 113.2) / 5.0) - ((V + V_x[ii] + 84.0) / 3.2))) / 3.7372);

    du = (u_inf - u) / tau_u;

    // Fast mechanism for submembranal Ca++ concentration (cai)
    drivechannel = -10.0 * (I_L + I_T) / (2.0 * 96489 * depth);
    if (drivechannel <= 0.0)
      drivechannel = 0.0;

    dcai = drivechannel + (cai_inf - cai) / tau_r;

    df[1 + (i - 1) * N] = (1.0 / C_m) * (-g_leak[ii] * (V - E_leak[ii]) - I_Na - I_Kd - I_M - I_L - I_T + (Iext[i] / Area[ii]));
    df[2 + (i - 1) * N] = dm;
    df[3 + (i - 1) * N] = dh;
    df[4 + (i - 1) * N] = dn;
    df[5 + (i - 1) * N] = dp;
    df[6 + (i - 1) * N] = dq;
    df[7 + (i - 1) * N] = dr;
    df[8 + (i - 1) * N] = dcai;
    df[9 + (i - 1) * N] = du;
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

float ran1(long *idum)
{
  int j;
  long k;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy)
  {
    if (-(*idum) < 1)
      *idum = 1;
    else
      *idum = -(*idum);
    for (j = NTAB + 7; j >= 0; j--)
    {
      k = (*idum) / IQ;
      *idum = IA * (*idum - k * IQ) - IR * k;
      if (*idum < 0)
        *idum += IM;
      if (j < NTAB)
        iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0)
    *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
}

float gasdev(long *idum)
{
  float ran1(long *idum);
  static int iset = 0;
  static float gset;
  float fac, rsq, v1, v2;

  if (*idum < 0)
    iset = 0;
  if (iset == 0)
  {
    do
    {
      v1 = 2.0 * ran1(idum) - 1.0;
      v2 = 2.0 * ran1(idum) - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  }
  else
  {
    iset = 0;
    return gset;
  }
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
