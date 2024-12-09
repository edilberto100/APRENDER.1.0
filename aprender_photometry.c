/**************************************************************
 *                           APPHI                            *
 * AN AUTOMATED PHOTOMETRY PIPELINE FOR HIGH CADENCE,         *
 *                     LARGE VOLUME DATA                      *
 * Edilberto Sanchez Moreno                                   *
 **************************************************************/
 #include "aprender.h"

 void get_centroids_list(int x, int y, MASTER_DATA * night, int N, int D, long int *naxes, float * xy, float * PX, float * PY, float * xyabs, double exptime)
 {
   //printf("|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|\n");
   //printf("Iniciando función get_centroids_list\n");
   //printf("Paramentros de entrada: x=%d, y=%d, N=%d, D=%d, naxes[0]=%ld, naxes[1]=%ld, naxes[2]=%ld\n",x,y,N,D,naxes[0],naxes[1],naxes[2]);
   long int i=0, j=0, k=0 , l=0, t=0, m=0;
   long int px=0, py=0;
   //printf("Paramentros de reserva: x=%d, y=%d EXPTIME=%lf\n", x, y, exptime);
   PX[0] = x;
   PY[0] = y;

   float *aux =  (float*) calloc(D*D,sizeof(float));
   for(k=0 ; k< N ; k++)
   {
        for(l=0 ; l< night[k].naxes[2] ; l++)
        {
        px = x - ((D-1)/2);
        py = y - ((D-1)/2);
        //printf("(px,py)=[%d,%d] ", x,y);
        if (px > night[k].naxes[0]-(D/2) || py > night[k].naxes[1]-(D/2) || px < (D/2) || py < (D/2))
        {
                //printf("Condición--->\n");
                px = (long int)x;
                py = (long int)y;
        }
        else
        {
            for(m=0 ; m<10 ; m++)
            {
                //P = 1.0;
                for(j=0 ; j<D ; j++)
                {
                    for(i=0 ; i<D ; i++)
                    {
                        if( (px+i) > night[k].naxes[0] || (py+j) > night[k].naxes[1] || (px+i < 0) || (py+j < 0))
                        {
                            aux[(j*D) + i] = 0.0;
                        }
                        else
                        {
                            aux[(j*D) + i] = night[k].data[(l*night[k].naxes[0]*night[k].naxes[1]) + (( (py) +j)*night[k].naxes[0]) + ( (px) +i)];
                        }
                    }
                }
                get_centroids_bin(aux,D,&px,&py);
                px = px - ((D-1)/2);
                py = py - ((D-1)/2);
            }
            get_centroids_bin(aux,D,&px,&py);
        }
        PX[t] = px;
        PY[t] = py;
        xy[t] = sqrt(pow(px-x,2)+pow(py-y,2));
        if(xy[t] > (D/2))
        {
            PX[t] = x;
            PY[t] = y;
        }
        //printf("t=%ld, (px,py)=[%d,%d] - (PX,PY)=[%ld,%ld] EXPTIME=%lf\n", t,x, y, px, py,exptime);
        }
   }
free(aux);
//printf("|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|\n");
return;
 }

 void get_naxes(MASTER_DATA * night, int N, long int *naxes)
{
    long int i;
    for(i=0 ; i<N ; i++)
    {
        naxes[2] += night[i].naxes[2];
    }

return;
}

float * get_st(MASTER_DATA * night, int N, int D, long int *naxes, float * PX, float * PY)
{

    long int i, j, k , l, t=0, m;
    long int x,y;
    float *star;
    long int px=0, py=0;
    naxes[0] = D;
    naxes[1] = D;

    for(i=0 ; i<N ; i++)
    {
        naxes[2] += night[i].naxes[2];
    }
   star = (float*) calloc(D*D,sizeof(float));
   for(k=0 ; k< N ; k++)
   {
        for(l=0 ; l< night[k].naxes[2] ; l++)
        {
            x=(long int)PX[t];
            y=(long int)PY[t];
            //printf ("x=%ld, y=%ld\n",x,y);
            px = x - ((D-1)/2);
            py = y - ((D-1)/2);
            //printf ("px=%ld, py=%ld\n",px,py);
            for(j=0 ; j<D ; j++)
            {
                for(i=0 ; i<D ; i++)
                {
                    //printf("%ld, %f\n",(j*D) + i, night[k].data[(l*night[k].naxes[0]*night[k].naxes[1]) + (( (py) +j)*night[k].naxes[0]) + ( (px) +i)]);
                    if( (px+i) > night[k].naxes[0] || (py+j) > night[k].naxes[1] || (px+i < 0) || (py+j < 0))
                    {
                        star[(j*D) + i] = 0.0;
                    }
                    else
                    {
                        star[(j*D) + i] = night[k].data[(l*night[k].naxes[0]*night[k].naxes[1]) + (( (py) +j)*night[k].naxes[0]) + ( (px) +i)];
                    }

                }
            }
        }
   }
return star;
}

void get_centroids_bin(float *bin, int B, long int *x, long int *y)
{
int i,j;
float vec_x[B];
float vec_y[B];
float mx = 0.0, my = 0.0;
long int dx = 0;
long int dy = 0;

/*for(j=0 ; j<B ; j++)
{
    for(i=0 ; i<B ; i++)
    {
       printf("%0.2f ",bin[(j*B)+i]);
    }
printf("\n");
}
printf("\n");//*/
//printf(" -- X -- | -- Y --\n");
for(j=0 ; j<B ; j++)
{
  vec_x[j] = 0;
  vec_y[j] = 0;
    for(i=0 ; i<B ; i++)
    {
      vec_x[j] += bin[(i*B)+j];
      vec_y[j] += bin[(j*B)+i];
    }
//printf("  %lf  |  %lf  \n",vec_x[j],vec_y[j]);
}
//printf(" ------- | -------\n");
float max_x=0.0, max_y=0.0;
for(j=0 ; j<B ; j++)
{
    if (vec_x[j] > max_x){
        max_x = vec_x[j];
        dx=j;
    }
    if(vec_y[j] > max_y){
        max_y=vec_y[j] ;
        dy=j;
    }
}

x[0] = x[0] + dx;
y[0] = y[0] + dy;

return ;
}

float * get_curve(float * star, float * sky , long int *naxes_star, long int *naxes_sky, long int D3)
{
float * curve = (float*) calloc(naxes_star[2],sizeof(float));
float *A= (float*) calloc(naxes_star[2],sizeof(float));
float *B= (float*) calloc(naxes_star[2],sizeof(float));
long int i, j, k,l=0;
//printf("Dim star= %ld, Dim Sky=%ld,D3=%ld\n",naxes_star[0],naxes_sky[0],D3);
int * circle = create_mask_circle(naxes_star[0]);
int * ring = create_mask_ring(naxes_sky[0],D3);

float a_R= 0.0;
float a_S= 0.0;

for(i=0 ; i<naxes_sky[0]; i++){
    for(j=0 ; j<naxes_sky[1]; j++){
        if(ring[(i*naxes_sky[1])+j] == 1)
            a_R= a_R + 1.0;
    }
}
for(i=0 ; i<naxes_star[0]; i++){
    for(j=0 ; j<naxes_star[1]; j++){
        if(circle[(i*naxes_star[1])+j] == 1)
            a_S= a_S + 1.0;
    }
}
//printf("Area_Ring = %0.2f, Area_Star = %0.2f\n",a_R, a_S);
for(k=0 ; k<naxes_star[2]; k++)
{
    A[k]=0.0;
    for(i=0 ; i<naxes_star[0]; i++){
        for(j=0 ; j<naxes_star[1]; j++){
            A[k]+= (star[(k*naxes_star[0] * naxes_star[1])+(i*naxes_star[1])+j] * circle[(i*naxes_star[1])+j]);
            //printf("%0.2f ",(star[(k*naxes_star[0] * naxes_star[1])+(i*naxes_star[1])+j] * circle[(i*naxes_star[1])+j]));
        }
        //printf("\n");
    }
    //printf("A[%ld] = %f \n",k,A[k]);
}
//printf("\n");
for(k=0 ; k<naxes_sky[2]; k++)
{
    B[k]=0.0;
    for(i=0 ; i<naxes_sky[0]; i++){
        for(j=0 ; j<naxes_sky[1]; j++){
            B[k]+= (sky[(k*naxes_sky[0] * naxes_sky[1])+(i*naxes_sky[1])+j] * ring[(i*naxes_sky[1])+j]);
            //printf("%0.2f ",(sky[(k*naxes_sky[0] * naxes_sky[1])+(i*naxes_sky[1])+j] * ring[(i*naxes_sky[1])+j]));
        }
        //printf("\n");
    }
    //printf("B[%ld] = %f \n",k,B[k]);
}
//printf("\n");
for(k=0 ; k<naxes_star[2]; k++)
{

// Para tiempo de integración
curve[k] =  (25.0-2.5*log10(A[k]-(a_S*(B[k]/a_R))))+2.5*log10(100);

// Para flujo
//curve[k] =  (A[k]-(a_S*(B[k]/a_R)));
//printf("curve[%ld] = %f \n",k,curve[k]);
}

//supr_outfl(curve,naxes_star[2]);
 free(circle);
 free(ring);
 free(A);
 free(B);
return curve;
}

float * get_mag(float * star, float * sky , long int *naxes_star, long int *naxes_sky, long int D3,int ti, float * err, float * errfluxes)
{
float * curve = (float*) calloc(naxes_star[2],sizeof(float));
float *A= (float*) calloc(naxes_star[2],sizeof(float));   // Fluxes of stars
float *B= (float*) calloc(naxes_star[2],sizeof(float));   // Fluxes of Sky
float *SD= (float*) calloc(naxes_star[2],sizeof(float));  // Standar Desviation
long int i, j, k,l=0;
//printf("Dim star= %ld, Dim Sky=%ld,D3=%ld\n",naxes_star[0],naxes_sky[0],D3);
int * circle = create_mask_circle(naxes_star[0]);
int * ring = create_mask_ring(naxes_sky[0],D3);

float a_R= 0.0;
float a_S= 0.0;

for(i=0 ; i<naxes_sky[0]; i++){
    for(j=0 ; j<naxes_sky[1]; j++){
        if(ring[(i*naxes_sky[1])+j] == 1)
            a_R= a_R + 1.0;
    }
}
for(i=0 ; i<naxes_star[0]; i++){
    for(j=0 ; j<naxes_star[1]; j++){
        if(circle[(i*naxes_star[1])+j] == 1)
            a_S= a_S + 1.0;
    }
}

for(k=0 ; k<naxes_star[2]; k++)
{
    A[k]=0.0;
    for(i=0 ; i<naxes_star[0]; i++){
        for(j=0 ; j<naxes_star[1]; j++){
            A[k]+= (star[(k*naxes_star[0] * naxes_star[1])+(i*naxes_star[1])+j] * circle[(i*naxes_star[1])+j]);
        }
    }
}
for(k=0 ; k<naxes_sky[2]; k++)
{
    B[k]=0.0;
    for(i=0 ; i<naxes_sky[0]; i++){
        for(j=0 ; j<naxes_sky[1]; j++){
            B[k]+= (sky[(k*naxes_sky[0] * naxes_sky[1])+(i*naxes_sky[1])+j] * ring[(i*naxes_sky[1])+j]);
        }
    }
}
// To Standar desviation
for(k=0 ; k<naxes_sky[2]; k++)
{
    SD[k]=0.0;
    float mean=B[k]/a_R;

    for(i=0 ; i<naxes_sky[0]; i++){
        for(j=0 ; j<naxes_sky[1]; j++){
            if(ring[(i*naxes_sky[1])+j] == 1)
                SD[k] += pow((((sky[(k*naxes_sky[0] * naxes_sky[1])+(i*naxes_sky[1])+j] * ring[(i*naxes_sky[1])+j])) - mean),2) ;
        }
    }
    SD[k] = pow((SD[k]/a_R),(0.5));
}

//printf("Fluxes   SD\n");
for(k=0 ; k<naxes_star[2]; k++)
{
float fluxes = (A[k]-(a_S*(B[k]/a_R)));
float epadu = 1.0;
curve[k] =  (25.0-2.5*log10(A[k]-(a_S*(B[k]/a_R))))+2.5*log10(ti);
err[k] = sqrt((fluxes/epadu) + (a_S*(SD[k]*SD[k])) + (a_S*a_S*(SD[k]*SD[k])/a_R));
err[k] =  1.0857 * err[k] / fluxes;
errfluxes[k] = fluxes/(SD[k]*sqrt(a_R));
//printf("errfluxes= %lf\n",errfluxes[k]);
}
//printf("\n");


 free(circle);
 free(ring);
 free(A);
 free(B);
 free(SD);
return curve;
}

float * get_fluxes(float * star, float * sky , long int *naxes_star, long int *naxes_sky, long int D3)
{
float * curve = (float*) calloc(naxes_star[2],sizeof(float));
float *A= (float*) calloc(naxes_star[2],sizeof(float));
float *B= (float*) calloc(naxes_star[2],sizeof(float));
long int i, j, k,l=0;
//printf("Dim star= %ld, Dim Sky=%ld,D3=%ld\n",naxes_star[0],naxes_sky[0],D3);
int * circle = create_mask_circle(naxes_star[0]);
int * ring = create_mask_ring(naxes_sky[0],D3);

float a_R= 0.0;
float a_S= 0.0;

for(i=0 ; i<naxes_sky[0]; i++){
    for(j=0 ; j<naxes_sky[1]; j++){
        if(ring[(i*naxes_sky[1])+j] == 1)
            a_R= a_R + 1.0;
    }
}
for(i=0 ; i<naxes_star[0]; i++){
    for(j=0 ; j<naxes_star[1]; j++){
        if(circle[(i*naxes_star[1])+j] == 1)
            a_S= a_S + 1.0;
    }
}
//printf("Area_Ring = %0.2f, Area_Star = %0.2f\n",a_R, a_S);
for(k=0 ; k<naxes_star[2]; k++)
{
    A[k]=0.0;
    for(i=0 ; i<naxes_star[0]; i++){
        for(j=0 ; j<naxes_star[1]; j++){
            A[k]+= (star[(k*naxes_star[0] * naxes_star[1])+(i*naxes_star[1])+j] * circle[(i*naxes_star[1])+j]);
            //printf("%0.2f ",(star[(k*naxes_star[0] * naxes_star[1])+(i*naxes_star[1])+j] * circle[(i*naxes_star[1])+j]));
        }
        //printf("\n");
    }
    //printf("A[%ld] = %f \n",k,A[k]);
}
//printf("\n");
for(k=0 ; k<naxes_sky[2]; k++)
{
    B[k]=0.0;
    for(i=0 ; i<naxes_sky[0]; i++){
        for(j=0 ; j<naxes_sky[1]; j++){
            B[k]+= (sky[(k*naxes_sky[0] * naxes_sky[1])+(i*naxes_sky[1])+j] * ring[(i*naxes_sky[1])+j]);
            //printf("%0.2f ",(sky[(k*naxes_sky[0] * naxes_sky[1])+(i*naxes_sky[1])+j] * ring[(i*naxes_sky[1])+j]));
        }
        //printf("\n");
    }
    //printf("B[%ld] = %f \n",k,B[k]);
}
//printf("\n");
for(k=0 ; k<naxes_star[2]; k++)
{

// Para tiempo de integración
//curve[k] =  (25.0-2.5*log10(A[k]-(a_S*(B[k]/a_R))))+2.5*log10(100);

// Para flujo
curve[k] =  (A[k]-(a_S*(B[k]/a_R)));
//printf("curve[%ld] = %f \n",k,curve[k]);
}

//supr_outfl(curve,naxes_star[2]);
 free(circle);
 free(ring);
 free(A);
 free(B);
return curve;
}

void supr_outfl(float *curva, long int naxes)
{
long int i, j;
float p=0.0, s=0.0;
float *aux;
aux =  (float*) calloc(naxes,sizeof(float));
float a,b;

for (i=0; i< naxes ; i++)
 {
   p+= curva[i];
   aux[i] = curva[i];
 }

p = p/(naxes );

for (i=0; i< naxes ; i++)
 {
   s+= pow((curva[i]-p),2) ;
 }
s = sqrt(s/(naxes -1));
a = p - (3.5*s);
b = p + (3.5*s);

for (i=0; i< naxes ; i++)
 {
       if(aux[i] == 0 )
       {
        curva[i] = p;
        aux[i] = p;
       }
 }

for (i=0; i< naxes ; i++)
 {



    if( aux[i] < a || aux[i] > b)
    {
        if(i > 4)
        {
          curva[i] = 0.0;
          for (j=1; j<6 ; j++)
          {
          curva[i]+= aux[i-j];
          }
          curva[i] = curva[i]/5;
        }
        else
        {
          curva[i] = 0.0;
          for (j=1; j<6 ; j++)
          {
          curva[i]+= (aux[i+j]);
          }
          curva[i] = curva[i]/5;
        }
    }

  }
free(aux);
return ;
}

float * get_diff_curve(float * star, float * comp, long int *naxes)
{
long int i;
float *curve;
curve =  (float*) calloc(naxes[2],sizeof(float));
for(i=0 ; i<naxes[2] ; i++)
    {
        curve[i] = (star[i]) - (comp[i]);
    }
return curve;
}

float * copy_vec(float * vec, long int *naxes)
{
long int i;
float *output;
output =  (float*) calloc(naxes[2],sizeof(float));
for(i=0 ; i<naxes[2] ; i++)
    {
    output[i] = vec[i];
    }

return output;
}

long int * size_snr(long int z, long int lz)
{
long int i;
float p, l, s;
long int *t_snr = (long int*) calloc((lz-1),sizeof(long int));
p = log10(z);
l=(6/(float)lz);
s=2;
for(i=0 ; i<(lz-1) ; i++)
  {
    t_snr[i] =(i+1)*10;// (long int)pow(s,p);
    //printf("%ld\n",(int)z/t_snr[i]);
    s= s+l;
  }
return t_snr ;
}

float ** sign_noise(float *star, float *sky, long int *naxes, long int D, long int *z, long int lz, long int D3)
{
    long int i,j,k,l,p;
    float **curves = (float**) calloc((lz-1),sizeof(float*));
    float *aux_star;
    float *aux_sky;
    long int naxes_sky[3];
    naxes_sky[0] = D;
    naxes_sky[1] = D;
    long int naxes_star[3];
    naxes_star[0] = naxes[0];
    naxes_star[1] = naxes[1];
   for (i=0 ; i<(lz-1) ; i++)
   {
      //printf("%ld\n",i);
       p = (long int)(naxes[2]/z[i]);
       naxes_sky[2] = p;
       naxes_star[2]=p;
       aux_star = (float*) calloc(naxes[0]*naxes[1]*p,sizeof(float));
       aux_sky = (float*) calloc(D*D*p,sizeof(float));
       for (l=0 ; l< p ; l++)
       {
            for (k=0 ; k< z[i] ; k++)
            {
                for (j=0 ; j< (naxes_star[0]*naxes_star[1]) ; j++)
                {
                aux_star[(l*(naxes_star[0]*naxes_star[1]))+j] += star[(l*z[i]*(naxes_star[0]*naxes_star[1]))+(k*(naxes_star[0]*naxes_star[1]))+j];
                }
            }
            for (k=0 ; k< z[i] ; k++)
            {
                for (j=0 ; j< (naxes_sky[0]*naxes_sky[1]) ; j++)
                {
                aux_sky[(l*(naxes_sky[0]*naxes_sky[1]))+j] += sky[(l*z[i]*(naxes_sky[0]*naxes_sky[1]))+(k*(naxes_sky[0]*naxes_sky[1]))+j];
                }
            }

       }
       curves[i] = get_curve(aux_star,aux_sky,naxes_star,naxes_sky,D3);
       free(aux_star);
       free(aux_sky);
   }
return curves;
}

float sign_curve(float *datos, long int z)
{
int j,i,k;
float prom = 0, std = 0, signal;
//printf("Z=%ld\n",z);

for (i=0 ; i<z   ; i++)
{
    prom += datos[i];
}
prom = prom/z;
//printf("prom %f\n",prom);

    for (i=0 ; i<z ; i++)
    {
     //printf("     i=%ld\n",i);
     std += pow(datos[i]-prom, 2);
    }

    std = pow((std/(z-1)),0.5);
    signal = prom/(2.5*std);
    //printf("%f\n",signal);
    return signal;
}

void free_STAR_DATA(STAR_DATA star,long int lz)
{
int i;
for (i=0 ; i<(lz-1) ; i++)
   {
   free(star.C_SNR[i]);
   }
free(star.err);
free(star.C_SNR);
free(star.star);
free(star.curve);
free(star.fluxes);
free(star.diff);
free(star.z);
free (star.SNR);
free (star.xy);
free (star.xyabs);
free(star.px);
free(star.py);
free(star.errfluxes);
}
