
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI acos(-1.0)

#define n_n 550000
#define NR_END 1
#define FREE_ARG char*
#define N 9           // numero de equacoes
#define NN 1000    // quantidade de HH
#define conexMAX  1.0*pow(10,3)   // each neuron max conections number
#define fe 0.8    // % excitatory neurons (80%)
#define neuron_classe 2   // number of neurons classes
#define transient 2500 //in ms

#define stim_delay 500.0
#define stim_dur 2000.0

float ran1(long *idum);
float gasdev(long *idum);

FILE *input;
FILE *input2;

 void derivs(double y[],double df[],double *Iext,double epsilon, double *Gexc,double *Gini);
 double *dvector(long nl,long nh);
 void free_dvector(double *v, long nl, long nh);
 void nrerror(char error_text[]);
 double **dmatrix(long nrl, long nrh, long ncl, long nch);
 int **imatrix(long nrl, long nrh, long ncl, long nch);
 void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
 int delta(int a, int b);
 void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
 int *vector(long nl,long nh);
 void free_vector(int *v, long nl, long nh);

 int main(void)
{

 int i,j,t,auxISI,sum,KM,numerodeconexoes,auxx,jj;
 int **listaexc,*conextotalex,**listaini,*conextotalin;
 double *x,tempo,step,*a,*b,*c,*df,*y,*Iext,I_uA,epsilon,pico[NN+2],somaISI,somaISI2,tpeak[NN+2],*Gexc,*Gini,delay,g_exc,razao,g_ini;
 double desviopadrao,CV;
 long idum = -12345678987;


 y=dvector(1,N*NN+1);
 df=dvector(1,N*NN+1);
 x=dvector(1,N*NN+1);
 a=dvector(1,N*NN+1);
 b=dvector(1,N*NN+1);
 c=dvector(1,N*NN+1);

 Iext=dvector(1,NN+1);

  Gexc=dvector(1,NN+1);
  Gini=dvector(1,NN+1);
   
  conextotalex=vector(1,NN+2);
  listaexc=imatrix(1,NN+1,1,conexMAX+2);

  conextotalin=vector(1,NN+2);
  listaini=imatrix(1,NN+1,1,conexMAX+2);

 step=0.01; // passo de integracao

  input=fopen("raster_HH_RS_g1_I0.00008_gexc0.0001_gM0.0_gL0.0_gT0.0.dat","wt");
  input2=fopen("HH.dat","wt");

 //*****************************************************************
for(g_exc=0.0001;g_exc<=0.00010000000002;g_exc=g_exc+0.0001) 		// 
for(I_uA=0.00008;I_uA<=0.0000800000000002;I_uA=I_uA+0.00003) 		// 
    {

epsilon=0;
//~ I_uA=0.00029; // depolarizing current pulse in uA
//~ g_exc=0.0009;
KM=100;
razao=1.0;
g_ini=razao*g_exc;
numerodeconexoes=NN*KM;//rede aleatória com KM //se KM=(NN-1); //global network	 	

 //****************** Criando lista com as conexoes *********************
//-----------zerando tudo------------//  
for(i=1;i<=NN;i++)
   {
   conextotalex[i]=0.0;
   conextotalin[i]=0.0;
   }

//--------pelo menos uma conexao para cada neuronio------//
sum=0.0; 
for(i=1;i<=fe*NN;i++)
   {
   conextotalex[i]=1.0;
   listaexc[i][1]=(int)NN*ran1(& idum)+1;
   sum=sum+1;   
   } 
for(i=fe*NN+1;i<=NN;i++)
   {
   conextotalin[i]=1.0;
   listaini[i][1]=(int)NN*ran1(& idum)+1;
   } 

//------- criando lista de conexoes exc---------// 
  while(sum<fe*numerodeconexoes) 
    {
    i=(int)NN*ran1(& idum)+1; //i neurônio pós sinaptico
    j=(int)NN*fe*ran1(& idum)+1;// j neurônio pré sinaptico
  
    auxx=0.0;      
    for(jj=1;jj<=conextotalex[j];jj++) //conta apenas novas conexoes
	if(listaexc[j][jj]==i) 
	auxx=1.0;

    if(auxx==0.0) 
	if(conextotalex[j]<=2*KM) 
	if(j!=i) 
	 {	
	 conextotalex[j]=conextotalex[j]+1;
	 listaexc[j][conextotalex[j]]=i; 
	 sum=sum+1;   
     }
    }
//------- criando lista de conexoes inh---------//
for(i=1;i<=NN;i++)
sum=sum+conextotalin[i];

  while(sum<numerodeconexoes) 
    {
    i=(int)NN*ran1(& idum)+1; //i neurônio pós sinaptico
    j=fe*NN+(int)NN*(1.0-fe)*ran1(& idum)+1;// j neurônio pré sinaptico
 
    auxx=0.0;
    for(jj=1;jj<=conextotalin[j];jj++) //conta apenas novas conexoes
	if(listaini[j][jj]==i) 
	auxx=1.0;

    if(auxx==0.0) 
    if(conextotalin[j]<=2*KM) 
	if(j!=i) 
	 {	
	 conextotalin[j]=conextotalin[j]+1;
	 listaini[j][conextotalin[j]]=i;
	 sum=sum+1;   
     }
     }
 /////////////////////////////////////////////
 // condicoes iniciais
    for(i=1;i<=NN;i++)
	  {
		x[1+(i-1)*N] = -20.0*ran1(& idum)-50.0;
		x[2+(i-1)*N] = 0.1*ran1(& idum);
		x[3+(i-1)*N] = 0.1*ran1(& idum);
		x[4+(i-1)*N] = 0.1*ran1(& idum);
		x[5+(i-1)*N] = 0.0;
		x[6+(i-1)*N] = 0.0;
		x[7+(i-1)*N] = 0.5;
		x[8+(i-1)*N] = 0.00024;
		x[9+(i-1)*N] = 0.0;
		
		Iext[i]=0.0; // depolarizing current pulse in uA
	    Iext[i]=I_uA+0.05*I_uA*gasdev(& idum);    
		pico[i]=0.0;
		
		Gexc[i]=0.0;
		Gini[i]=0.0;
		
		auxISI = 0;
		somaISI = 0;
		somaISI2 = 0;
		tpeak[i] = -10.0;
	  }


  tempo=0.0;

////////////////////////////////////////////////////////////loop do tempo
 for(t=1;t<=n_n;t++)
    {   
     tempo=tempo+step;

      for(i=1;i<=N*NN;i++)
        y[i]=x[i];
                
     if(tempo<=1000.0)  //para dessinctonizar
		delay =  2.0*ran1(& idum); 
     if(tempo>1000.0 && tempo<=1000.0+step)  //para dessinctonizar
        delay=0.01;
     
	  for(i=1;i<=NN;i++)  // tempo anterior
		{         
			//~ if(i<=fe*NN) 
			  //~ delay=0.01;
			//~ else 
			  //~ delay=0.01;	
				
			if(tempo>=tpeak[i]+delay && tempo<tpeak[i]+delay+step)  //1.5  0.8 ms
			{
			if(i<=fe*NN) 
			  for(j=1;j<=conextotalex[i];j++) 
				Gexc[listaexc[i][j]]=Gexc[listaexc[i][j]]+g_exc; // condutancia exc recebida
			else 
			  for(j=1;j<=conextotalin[i];j++) 
				Gini[listaini[i][j]]=Gini[listaini[i][j]]+g_ini; // condutancia inhi recebida
			}
		 }
             
 // ------------ Runge-Kutta --------------------         
     derivs(y,df,Iext,epsilon,Gexc,Gini);  
     for(i=1;i<=N*NN;i++)
           {
            a[i]=step*df[i];
            y[i]=x[i]+a[i]/2.0;
           }
     derivs(y,df,Iext,epsilon,Gexc,Gini);
     for(i=1;i<=N*NN;i++)
           {
            b[i]=step*df[i];
            y[i]=x[i]+b[i]/2.0;
	   }
     derivs(y,df,Iext,epsilon,Gexc,Gini);
     for(i=1;i<=N*NN;i++)
           {
            c[i]=step*df[i];
            y[i]=x[i]+c[i]; 
	   }
     derivs(y,df,Iext,epsilon,Gexc,Gini);
     for(i=1;i<=N*NN;i++)
        x[i]=x[i]+(a[i]+step*df[i])/6.0+(b[i]+c[i])/3.0;  
  // ----------------------------------------------
         
	  //~ if(tempo>0.0)
	  //~ {	
		//~ fprintf(input,"%.2f",tempo-0.0); //nA
		//~ for(i=1;i<=NN;i++)  
        //~ fprintf(input," %.2f",x[1+(i-1)*N]);
		//~ fprintf(input,"\n");
	 //~ }

     for(i=1;i<=NN;i++)  
       {		   
		 if(x[1+(i-1)*N]<-20.0)
		   pico[i]=20.0;	 
	
		 if(x[1+(i-1)*N]>pico[i])  //marca o disparo quando o potencial V_i passa de 0 mV
		   {
		     pico[i]=1000;     
	
		     if(tempo>transient && tpeak[i]>0.0)
		       {
				  //~ if(auxISI[i]==0)
	                somaISI=somaISI+(tempo-tpeak[i]);
	                somaISI2=somaISI2+(tempo-tpeak[i])*(tempo-tpeak[i]);
	                
	                auxISI=auxISI+1;
		       }
	            tpeak[i]=tempo;    
	            fprintf(input,"%.2f %d\n",tempo,i);//raster plot
	  	    }	  	    
	  	    
		  Gexc[i]=Gexc[i]*exp(-step/5.0);  // atualiza a condutancia recebida	
		  Gini[i]=Gini[i]*exp(-step/5.0);  // 
       }   

   }//fim do loop do tempo

 
	//~ fprintf(input,"\n");//raster plot linha

// CALCULO DO CV 
  desviopadrao=sqrt(somaISI2/auxISI - (somaISI/auxISI)*(somaISI/auxISI));
  CV=desviopadrao/(somaISI/auxISI);           


//~ fprintf(input,"%f %f %f %.2f %.3f %.3f\n",I_uA,g_ini,g_exc,razao,CV,1000.0*auxISI/somaISI);
printf("%f %f %f %.2f %.3f %.3f\n",I_uA,g_ini,g_exc,razao,CV,1000.0*auxISI/somaISI);

}
   free_dvector(y,1,N*NN+1);
   free_dvector(df,1,N*NN+1); 
   free_dvector(x,1,N*NN+1);
   free_dvector(a,1,N*NN+1);
   free_dvector(b,1,N*NN+1);
   free_dvector(c,1,N*NN+1);

   free_dvector(Iext,1,N*NN+1);
   
  free_dvector(Gexc,1,NN+1); 
  free_dvector(Gini,1,NN+1); 

  free_vector(conextotalex,1,NN+2);
  free_vector(conextotalin,1,NN+2);
  free_imatrix(listaexc,1,NN+1,1,conexMAX+2);
  free_imatrix(listaini,1,NN+1,1,conexMAX+2);

   return 0;

    
  fclose(input);
  fclose(input2);
}

////////// equacoes 
void derivs(double y[],double df[],double *Iext,double epsilon, double *Gexc,double *Gini)
{	
 int i,ii,neurontype;
 double C_m,E_Na,E_K,E_Ca;
 double V_T[neuron_classe+1],E_leak[neuron_classe+1],g_leak[neuron_classe+1],g_Na[neuron_classe+1],g_Kd[neuron_classe+1],g_M[neuron_classe+1],tau_max[neuron_classe+1];
 double V,m,h,n,p,I_Na,alpha_m,beta_m,alpha_h,beta_h,dm,dh,I_Kd,alpha_n,beta_n,dn,I_M,p_inf,tau_p,dp;
 double Area[neuron_classe+1]; 
 double V_x[neuron_classe+1],g_L[neuron_classe+1],g_T[neuron_classe+1];	   	  
 double q,r,I_L,alpha_q,beta_q,alpha_r,beta_r,dq,dr,s_inf,I_T,u,u_inf,tau_u,du;	   	  
 double cai,dcai,drivechannel,depth,cai_inf,tau_r;
 double Vr_exc,Vr_ini;
 
      C_m = 1.0; // uF/cm^2
	  E_Na = 50.0; // mV
	  E_K = -100.0; //-90.0 mV no paper // -100.0 mV no modelDB
	  depth = 1.0;
	  cai_inf = 0.00024;
	  tau_r = 50.0;  
  Vr_exc=0.0;      // mV Potencial reverso exc 
  Vr_ini=-80.0;      // mV Potencial reverso inh
	  
  //RS - adaptado
      neurontype = 1;      
      Area[neurontype] = 0.0096*0.0096*PI; // cm^2
	  V_T[neurontype] = -55.0; // mV 
      E_leak[neurontype] = -70.0; // mV            
      g_leak[neurontype] = 0.02; // mS/cm^2
	  g_Na[neurontype] = 50.0; // mS/cm^2	  
	  g_Kd[neurontype] = 5.0; // mS/cm^2		  
	  g_M[neurontype] = 0.0; // mS/cm^2	 
	  tau_max[neurontype] = 1000.0; // ms	
      V_x[neurontype] = 2.0; // mV            
      g_L[neurontype] = 0.0; // mS/cm^2
	  g_T[neurontype] = 0.0; // mS/cm^2	
	    
  //FS - adaptado
      neurontype = 2;      
      Area[neurontype] = 0.0096*0.0096*PI; // cm^2
	  V_T[neurontype] = -55.0; // mV 
      E_leak[neurontype] = -70.0; // mV            
      g_leak[neurontype] = 0.02; // mS/cm^2
	  g_Na[neurontype] = 50.0; // mS/cm^2	  
	  g_Kd[neurontype] = 5.0; // mS/cm^2		  
	  g_M[neurontype] = 0.01; // mS/cm^2	 
	  tau_max[neurontype] = 1000.0; // ms	
      V_x[neurontype] = 2.0; // mV            
      g_L[neurontype] = 0.1; // mS/cm^2
	  g_T[neurontype] = 0.1; // mS/cm^2		  
	  

//////////////////////////EQUAÇOES ACOPLADAS////////////////////  
  for(i=1;i<=NN;i++)
    {
	  if(i<=fe*NN)  
        ii=1;  //neuron type
	  else  
        ii=1;  
      
	  V = y[1+(i-1)*N];
	  m = y[2+(i-1)*N];
	  h = y[3+(i-1)*N];	
	  n = y[4+(i-1)*N];
	  p = y[5+(i-1)*N];
	  q = y[6+(i-1)*N];
	  r = y[7+(i-1)*N];
	  cai = y[8+(i-1)*N]; // concentration of the Calcium (Ca2+) in the intracellular fluid in mM
	  u = y[9+(i-1)*N];
	  	   	   
	//Calcium reversal potential		
	  //~ E_Ca = 120.0; // mV         
      E_Ca = 13.31854*log(2.0/cai);         
	  	    
	//sodium current responsible for action potentials		
	  I_Na = g_Na[ii]*pow(m,3)*h*(V-E_Na);   	  
      alpha_m = (-0.32*(V-V_T[ii]-13.0))/(exp(-(V-V_T[ii]-13.0)/4.0)-1.0);  
      beta_m = (0.28*(V-V_T[ii]-40.0))/(exp((V-V_T[ii]-40.0)/5.0)-1.0);   	  
      alpha_h = 0.128*exp(-(V-V_T[ii]-17.0)/18.0);                    
      beta_h = 4.0/(1.0+exp(-(V-V_T[ii]-40.0)/5.0));      
      
      dm = alpha_m*(1.0-m)-beta_m*m;
      dh = alpha_h*(1.0-h)-beta_h*h;	  
	  
	 //potassium current responsible for action potentials	  
	  I_Kd = g_Kd[ii]*(pow(n,4)*(V-E_K));                      
      alpha_n = (-0.032*(V-V_T[ii]-15.0))/(exp(-(V-V_T[ii]-15.0)/5.0)-1.0); 
      beta_n = 0.5*exp(-(V-V_T[ii]-10.0)/40.0);    
      
      dn = alpha_n*(1.0-n)-beta_n*n;      
      
	 //slow voltage-dependent potassium current responsible for spike-frequency adaptation	   
	  I_M = g_M[ii]*p*(V-E_K);    		     
      p_inf = 1.0/(1.0+exp(-(V+35.0)/10.0));                     
      tau_p = tau_max[ii]/(3.3*exp((V+35.0)/20.0)+exp(-(V+35.0)/20.0));   
	  
	  dp = (p_inf-p)/tau_p;	 	  
            
	  // low-threshold Calcium currents to generate bursting  	  
	  I_L = g_L[ii]*pow(q,2)*r*(V-E_Ca);    //E_Ca = 120.0 no code // not causing [Ca2+] influx

	  if(V==-27.0) //evitar ()/0.0
	  V=-27.001;
	  	  	  
      alpha_q = (0.055*(-27.0-V))/(exp((-27.0-V)/3.8)-1.0);     
      beta_q = 0.94*exp((-75.0-V)/17.0);   	     
      alpha_r = 0.000457*exp((-13.0-V)/50.0);                                
      beta_r = 0.0065/(exp((-15.0-V)/28.0)+1.0);   
                        
      dq = alpha_q*(1.0-q)-beta_q*q;
      dr = alpha_r*(1.0-r)-beta_r*r;	
      
	  
	  // high-threshold	Calcium currents to generate bursting  
	  s_inf = 1.0/(1.0+exp(-(V+V_x[ii]+57.0)/6.2));     
	  	   
	  I_T = g_T[ii]*pow(s_inf,2)*u*(V-E_Ca);    
	  
	  u_inf = 1.0/(1.0+exp((V+V_x[ii]+81.0)/4.0));	
	  	  
	  if(V<0.0) // approximation to "exp()<10^32"
	  tau_u = 30.8/3.7372+(211.4+exp((V+V_x[ii]+113.2)/5.0))/(3.7372*(1.0+exp((V+V_x[ii]+84.0)/3.2)));  
	  else
	  tau_u = 30.8/3.7372+((exp(((V+V_x[ii]+113.2)/5.0)-((V+V_x[ii]+84.0)/3.2)))/3.7372);  
	  
      du = (u_inf-u)/tau_u;
	  
	  
	  // Fast mechanism for submembranal Ca++ concentration (cai) 	        
      drivechannel = -10.0*(I_L+I_T)/(2.0*96489*depth); //
      //~ drivechannel = 0.0;
      if(drivechannel<=0.0)
       drivechannel = 0.0;
       
      dcai = drivechannel+(cai_inf-cai)/tau_r;    
	  	  

      df[1+(i-1)*N] = (1.0/C_m)*(-g_leak[ii]*(V-E_leak[ii]) - I_Na - I_Kd - I_M - I_L - I_T + (Iext[i]/Area[ii])+ Gexc[i]*(Vr_exc-V) + Gini[i]*(Vr_ini-V));       
      df[2+(i-1)*N] = dm;
      df[3+(i-1)*N] = dh;
      df[4+(i-1)*N] = dn;      
      df[5+(i-1)*N] = dp;      
      df[6+(i-1)*N] = dq;
      df[7+(i-1)*N] = dr;      
      df[8+(i-1)*N] = dcai;      
      df[9+(i-1)*N] = du;
   }
}
/////////////////////////////////////////


double *dvector(long nl,long nh)
{
   double *v;
   
   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}


void free_dvector(double *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}


float ran1(long *idum)
{
 int j;
 long k;
 static long iy=0;
 static long iv[NTAB];
 float temp;
 
 if(*idum<=0 || !iy)
   {
     if(-(*idum)<1) *idum=1;
     else *idum = -(*idum);
    for(j=NTAB+7;j>=0;j--)
      {
       k=(*idum)/IQ;
       *idum=IA*(*idum-k*IQ)-IR*k;
       if(*idum<0) *idum +=IM;
       if(j<NTAB) iv[j]=*idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if(*idum<0) *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j]=*idum;
   if((temp=AM*iy)>RNMX) return RNMX;
   else return temp;
}

float gasdev(long *idum)
{
  float ran1(long *idum);
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;  
  
  if(*idum<0) iset=0;
  if(iset==0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq>=1.0 || rsq==0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
} 

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

   m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m += NR_END;
   m -= nrl;

   m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
   if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}


int **imatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   int **m;

   m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m += NR_END;
   m -= nrl;

   m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
   if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   return m;
}


void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}


int *vector(long nl,long nh)
{
   int *v;
   
   v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}

void free_vector(int *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}

