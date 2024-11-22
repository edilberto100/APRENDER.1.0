/**************************************************************
 *                           APPHI                            *
 * AN AUTOMATED PHOTOMETRY PIPELINE FOR HIGH CADENCE,         *
 *                     LARGE VOLUME DATA                      *
 * Edilberto Sanchez Moreno                                   *
 **************************************************************/
 #include "aprender.h"

STAR_BIN read_position(struct arguments arguments,int * nobj, int *naper)
 {
 STAR_BIN comp, data;
 int nist, ndiam,i,j;
 comp = read_lst(arguments.ist,&nist);
 long int * diameter= get_diameter(arguments.d1,&ndiam);
 /*for(i=0; i<ndiam;i++)
   printf("%ld \n",diameter[i]);*/
 naper[0]=ndiam;
 nobj[0]=comp.cant;

 data.cant=ndiam*comp.cant;
 data.bin=(long int*) calloc(data.cant,sizeof(long int));
 data.x=(long int*) calloc(data.cant,sizeof(long int));
 data.y=(long int*) calloc(data.cant,sizeof(long int));
 data.xe=(long int*) calloc(data.cant,sizeof(long int));
 data.ye=(long int*) calloc(data.cant,sizeof(long int));
 data.x0=(long int*) calloc(data.cant,sizeof(long int));
 data.y0=(long int*) calloc(data.cant,sizeof(long int));
 data.D2=(long int*) calloc(data.cant,sizeof(long int));
 data.D3=(long int*) calloc(data.cant,sizeof(long int));
 data.dt=(double*) calloc(data.cant,sizeof(double));

 for(i=0; i<comp.cant;i++)
   {
     for(j=0; j<ndiam;j++)
       {
           data.x[(j*comp.cant)+i]=comp.x[i];
           data.y[(j*comp.cant)+i]=comp.y[i];
           data.xe[(j*comp.cant)+i]=comp.xe[i];
           data.ye[(j*comp.cant)+i]=comp.ye[i];
           data.x0[(j*comp.cant)+i]=comp.x[i];
           data.y0[(j*comp.cant)+i]=comp.y[i];
           data.bin[(j*comp.cant)+i]=diameter[j];
           data.D2[(j*comp.cant)+i]=arguments.d3;
           data.D3[(j*comp.cant)+i]=arguments.d2;
           data.dt[(j*comp.cant)+i]=comp.dt[i];
       }
   //  printf("%ld \n",diameter[i]);
   }
  free_STAR_BIN(comp);
  free(diameter);

 return data;
 }

STAR_BIN conv_img(long int *data, float umbrall, long int *naxes, int D , float *comp)
{
float *FrameTotal = (float*) calloc(naxes[0]*naxes[1],sizeof(float));
int i,j;
FILE_FITS Dato;
Dato.data = (long int*) calloc(naxes[0]*naxes[1],sizeof(long int));
Dato.naxes[0] = naxes[0];
Dato.naxes[1] = naxes[1];
Dato.naxes[2] = 1;
Dato.dim = 3;
long int max = FrameTotal[(1*naxes[0])+1], min = FrameTotal[(1*naxes[0])+1];

for(j=0  ; j< naxes[1] ; j++)
{
    for(i=0 ; i< naxes[0] ; i++)
    {
        FrameTotal[(j*naxes[0])+i] = (long int)((log(data[(j*naxes[0])+i]))+0.5) ;
        Dato.data[(j*naxes[0])+i] = (long int)((log(data[(j*naxes[0])+i]))+0.5) ;
    }
}

float maximo= get_maximum(FrameTotal,naxes[0]*naxes[1]);
float minimo= get_minimum(FrameTotal,naxes[0]*naxes[1]);
//printf("\nMaximo: %f, Minimo: %f, \n",maximo,minimo);

for(j=0  ; j< naxes[1] ; j++)
{
    for(i=0 ; i< naxes[0] ; i++)
    {
       if(FrameTotal[(j*naxes[0])+i] > (minimo + (maximo - minimo)/2) )
       {
        FrameTotal[(j*naxes[0])+i] = 1.0 ;
        Dato.data[(j*naxes[0])+i] = 1.0 ;
       }
       else
       {
        FrameTotal[(j*naxes[0])+i] = 0.0 ;
        Dato.data[(j*naxes[0])+i] = 0.0;
       }
    }
}

save_fits_ld(Dato,"segmentation.fits");
free(Dato.data);

CORD ftp;
ftp=get_perfil(FrameTotal,naxes);
STAR_BIN pstars, stars, salida;
long int cantx, canty;
pstars.x=filter_positions(ftp.x,naxes[0],0.9,&cantx,1);
pstars.y=filter_positions(ftp.y,naxes[1],0.9,&canty,2);
pstars.cant=cantx*canty;
PROFILE tama;
tama.x=cantx;
tama.y=canty;
stars=cut_star(pstars,tama,FrameTotal,naxes);
salida=obtain_stars(stars,FrameTotal,naxes,15);
free(stars.x);
free(stars.y);
free(FrameTotal);
free_CORD(ftp);
return salida;
}

/*********************** Obtener perfiles ***************************/
//Funcion
CORD get_perfil(float* frame,long int* naxes)
{
CORD perfil;
perfil.x=(float*) calloc(naxes[0], sizeof(float));
perfil.y=(float*) calloc(naxes[1], sizeof(float));
int i, j;
for(j=0  ; j< naxes[1] ; j++)
{
    for(i=0 ; i< naxes[0] ; i++)
    {
    perfil.x[i]+= frame[i+(j*naxes[0])];
    perfil.y[j]+= frame[i+(j*naxes[0])];
    }
    //for(i=0 ; i< naxes[0] ; i++)
    //{
    //perfil.y[j]= perfil.y[j]+frame[i+(j*naxes[0])];
    //}
}
//for(j=0  ; j< naxes[1] ; j++)
//printf("%f, ",perfil.y[j]);
//printf("\n\n");
//for(j=0  ; j< naxes[0] ; j++)
//printf("%f, ",perfil.x[j]);
//printf("\n\n");
return perfil;
}
/*********************** Máximo, mínimo y ancho libre medio ***************************/
//Funcion
float get_maximum(float* P,long int N)
{
long int i=0;
float M=P[0];
for (i=1; i<N;i++)
{
    if(P[i] > M)
     M = P[i];
}
return M;
}
/*********************** Máximo, mínimo y ancho libre medio ***************************/
//Funcion
float get_minimum(float* P,long int N)
{
long int i;
float M=P[0];
for (i=1; i<N;i++)
{
    if(P[i] < M)
     M = P[i];
}
return M;
}
/*********************** Promedios x, y, Total y Desv. Est.***************************/
float get_avg(float *V,long int N)
{
long int i;
float avg=0;

for(i=0 ; i< N ; i++)
    {
        avg+=(V[i]/(N));
    }
    return avg;
}

/***********************************************************************/
float get_moda(float *V , long int N)
{
/*long int i,j, v, p;
long int *aux = (long int *) calloc(N,sizeof(long int));
    for(i=0;i<N;i++)
    {
        v = (long int)V[i];
        p = i;
        for(j=i;j<N;j++) {
            if(V[j]==v)
            aux[p]++;
        }
    }

    long int max = aux[0];
    long int pm = 0;
    for(i=0;i<N;i++) {
        if(aux[i]>max) {
            pm=i;
            max=aux[i];
        }
    }
    // Visualizar el elemento con mas frecuencia de aparicion
    //printf("\nMODA : %f",V[pm]);
    free(aux);*/
    V[0] = 7.0;
    return V[0];

}

/*********************** Normalizar Perfiles ***************************/
//Función
void normal_vec(float* P, long int N)
{
long int i, j;
float prom=get_avg(P,N);
float max=get_maximum(P,N);
float moda = get_moda(P,N);

for(i=0 ; i< N ; i++)
    {
        //P[i]=((P[i]/prom)-1)/((max/prom)-1);
        P[i]= P[i]-moda; //((P[i]/prom)-1)/((max/prom)-1);
    }

}

/*********************** Filtrado Pasa altas ***************************/
void filter_hp(float* P,long int N, float Umb)
{
long int i;
normal_vec(P,N);
for (i=0; i<N;i++)
{
    //if (P[i]>Umb)
    if (P[i]>0)
    {
     P[i]=P[i];//(P[i]-Umb);
    }
    else
    {
     P[i]=0;
    }
}
}

/*********************** Obtener perfiles ***************************/
long int* filter_positions(float* P,long int N, float Umb, long int *Cant, long int A)
{
//filter_hp(P,N,Umb);
long int* Peaks=(long int*) calloc(N,sizeof(long int));
long int* Ps=(long int*) calloc(N,sizeof(long int));
long int h, i, k, C=0;
float g,l,v;
h=0;
//printf("condicion\n\n");
for(i=1; i<N-1; i++)
{
    g=P[i-1]-P[i];
    l=P[i+1]-P[i];
    v=P[i-1]-P[i+1];
    //printf("%ld  -- g=%f, l=%f, v=%f\n",i,g,l,v);
    if(g <= 0.0 || l == 0 )
        if(P[i] > 0 && i> 16 && i< (N-17))
        {
            Ps[i] = i;
            Peaks[h]=i;
            h++;
        }
}
//printf("\n\n");
//for(i=0  ; i< N ; i++)
//printf("%f, %ld\n",P[i],Ps[i]);
//printf("\n\n");

free(Ps);
C=0;
for (i=0; i<N-1; i++)
{
    float dif = Peaks[i+1]-Peaks[i] ;
    if((dif) >= 2 ||  (dif) < 0 )
    {
        C ++;
    }
}
long int * maximos=(long int *) calloc(C,sizeof(long int));
Cant[0] = C;
k=-1;
for (i=0; i<N-1; i++)
{
    float dif = Peaks[i+1]-Peaks[i] ;
    if(dif >= 2.0 ||  dif < 0.0 )
    {
        k++;
        maximos[k]=Peaks[i];
    }
}
free(Peaks);
//printf("maximos %ld\n",A);
//for(i=0  ; i<C ; i++)
//printf("%ld\n",maximos[i]+1);
return maximos;
}

/*********************** Obtener perfiles ***************************/
STAR_BIN cut_star(STAR_BIN pst, PROFILE tam, float* FrameTotal, long int* naxes)
{

long int h, i, j;
STAR_BIN star;
star.x=(long int*) calloc(pst.cant,sizeof(long int));
star.y=(long int*) calloc(pst.cant,sizeof(long int));
star.cant=0;
//Ubicar estrellas
float moda;
moda= 0.0; //get_moda(FrameTotal, (naxes[0])*(naxes[1]));
for(j=0 ; j< tam.y ; j++)
{
    for(i=0 ; i< tam.x ; i++)
    {
       float FT=FrameTotal[pst.y[j]*naxes[0] + pst.x[i]];
       if(FT>moda)
       {
        star.x[star.cant]=pst.x[i];
        star.y[star.cant]=pst.y[j];
        star.cant++;
       }
    }
}

free(pst.x);
free(pst.y);
return star;
}

/*********************** Obtener perfiles ***************************/
STAR_BIN obtain_stars(STAR_BIN stars, float* FrameTotal, long int* naxes, long int box)
{
long int i, j, k, lim;

float **cajas;
cajas=(float**) calloc(stars.cant,sizeof(float*));

for(i=0  ; i<stars.cant ; i++)
{
cajas[i]=(float*) calloc(box*box,sizeof(float));
}

lim=box/2;

for(k=0; k<stars.cant ; k++)
{
for(j=0; j<box ; j++)
{
    for(i=0 ; i<box; i++)
    {
        cajas[k][i+(j*box)]=FrameTotal[(stars.x[k]+i-lim)+((stars.y[k]+j-lim)*naxes[0])];
    }
}
}

float **bol;
bol=(float**) calloc(stars.cant,sizeof(float*));

for(i=0  ; i<stars.cant ; i++)
{

bol[i]=(float*) calloc(box*box,sizeof(float));

}

for(k=0  ; k<stars.cant ; k++)
{
bol[k][0] = (cajas[k][1] + cajas[k][box] + cajas[k][box+1])/3;
bol[k][box-1] = (cajas[k][box-2] + cajas[k][(box*2)-1] + cajas[k][(box*2)-2])/3;
bol[k][(box*(box-1))] = (cajas[k][(box*(box-2))] + cajas[k][(box*(box-2))+1] + cajas[k][(box*(box-1))+1] )/3;
bol[k][((box*box)-1)] = (cajas[k][(box*box)-2] + cajas[k][(box*(box-1))-1] + cajas[k][(box*(box-1))-2] )/3;

for(i=1  ; i<box-1 ; i++)
{
j=1;
bol[k][i] = (cajas[k][i-1] + cajas[k][i+1] + cajas[k][box*j+(i-1)] + cajas[k][box*j+(i)] + cajas[k][box*j+(i+1)])/5;
}

for(j=1  ; j<box-1 ; j++)
{
i=1;
bol[k][j*box] = (cajas[k][(j-1)*box] + cajas[k][(j+1)*box] + cajas[k][((j-1)*box)+i] + cajas[k][((j)*box)+i] + cajas[k][((j+1)*box)+i])/5;
}

for(i=1  ; i<box-1 ; i++)
{
j=box-2;
bol[k][((j+1)*box)+i] = (cajas[k][((j+1)*box)+(i-1)] + cajas[k][((j+1)*box)+(i+1)]  + cajas[k][(j*box)+(i-1)]  + cajas[k][(j*box)+i]  + cajas[k][(j*box)+(i+1)] )/5;
}

for(j=1  ; j<box-1 ; j++)
{
i=box-2;
bol[k][(j*box)+(i+1)] = (cajas[k][((j-1)*box)+(i+1)] + cajas[k][((j+1)*box)+(i+1)]  + cajas[k][((j-1)*box)+i]  + cajas[k][(j*box)+i]  + cajas[k][((j+1)*box)+i] )/5;
}

for(j=1  ; j<box-1 ; j++)
{

    for(i=1 ; i<box-1; i++)
    {
        bol[k][i+(j*box)]=(cajas[k][(i-1)+(j-1)*box]+cajas[k][(i)+(j-1)*box]+cajas[k][(i+1)+(j-1)*box]+cajas[k][(i-1)+(j)*box]+cajas[k][(i+1)+(j)*box]+cajas[k][(i-1)+(j+1)*box]+cajas[k][(i)+(j+1)*box]+cajas[k][(i+1)+(j+1)*box])/8;

    }
}
}
for (k=0; k<stars.cant; k++)
{
free(cajas[k]);
}

free(cajas);

long int dx,dy;
STAR_BIN Salida;
Salida.cant=stars.cant;
Salida.x=(long int*) calloc(stars.cant,sizeof(long int));
Salida.y=(long int*) calloc(stars.cant,sizeof(long int));
Salida.bin=(long int*) calloc(stars.cant,sizeof(long int));

//Agregado Edilberto
Salida.D2= (long int*) calloc(stars.cant,sizeof(long int));
Salida.D3= (long int*) calloc(stars.cant,sizeof(long int));
/**/

for(k=0;k<stars.cant;k++)
{
dx=0;
dy=0;
get_centroids_bin(bol[k],box,&dx,&dy );
//printf("dx = %ld, dy = %ld, val %f\n",dx,dy,bol[k][(dy*box)+dx]);
Salida.x[k]=(stars.x[k]-lim)+dx;
Salida.y[k]=(stars.y[k]-lim)+dy;
}

CORD* PerfilCajas;
PerfilCajas=(CORD*) calloc(stars.cant,sizeof(CORD));
long int dimcaja[2];
dimcaja[0]=box;
dimcaja[1]=box;
long int c=0;

for(k=0  ; k<stars.cant ; k++)
    {
    PerfilCajas[k]=get_perfil(bol[k],dimcaja);
    }
for (k=0; k<stars.cant; k++)
{
free(bol[k]);
}
free(bol);
PROFIL* maxicajas;
maxicajas = (PROFIL*) calloc(stars.cant, sizeof(PROFIL));
PROFIL* minicajas;
minicajas = (PROFIL*) calloc(stars.cant, sizeof(PROFIL));
for (k=0; k<stars.cant; k++)
{
maxicajas[k].x=get_maximum(PerfilCajas[k].x, dimcaja[0]);
maxicajas[k].y=get_maximum(PerfilCajas[k].y, dimcaja[1]);
minicajas[k].x=get_minimum(PerfilCajas[k].x, dimcaja[0]);
minicajas[k].y=get_minimum(PerfilCajas[k].y, dimcaja[1]);
float hx,hy;
hx= maxicajas[k].x-((maxicajas[k].x-minicajas[k].x)/2);
hy= maxicajas[k].y-((maxicajas[k].y-minicajas[k].y)/2);

int Cx=0, Cy=0;
for (i=0; i<box;i++)
{
    if(PerfilCajas[k].x[i]> hx)
    {
     Cx++;
    }
    if(PerfilCajas[k].y[i] > hy)
    {
     Cy++;
    }
}

free_CORD(PerfilCajas[k]);

if (Cx>=Cy)
{
Salida.bin[k]=2*Cx;
Salida.D2[k]= Salida.bin[k] *2.28;
Salida.D3[k]= Salida.bin[k] *1.55;
}
else
{
Salida.bin[k]=2*Cy;
Salida.D2[k]= Salida.bin[k] *2.28;
Salida.D3[k]= Salida.bin[k] *1.55;
}
}




free(PerfilCajas);
free(maxicajas);
free(minicajas);
return Salida;
}

void free_CORD(CORD v)
{
free(v.x);
free(v.y);
}

float * create_gaussean(int d1)
{
  int i,j;
  float *Y = (float*) calloc(d1*d1,sizeof(float));
  float *X = (float*) calloc(d1,sizeof(float));
  float e=0.4342944819;
  float N=d1/(d1);
  float M = d1/2;
  float V=pow((2*0.85),2);
  for (j=0; j<d1; j++)
  {
  X[j] = j;
  //printf("%f \n",  X[j] );
  }

  for (j=0; j<d1; j++)
  {
    for (i=0; i<d1; i++)
    {
    float a = pow((-X[i]+M),2)/V;
    float b = pow((-X[j]+M),2)/V;
    Y[(j*d1) + i] = pow(e,(a)+(b));
    //printf("%f ", Y[(j*d1) + i]);
    }
    //printf("\n");
  }
free(X);
return Y;
}

float * create_mask_circle(int d1)
{
  int i,j;
  float *maskC = (float*) calloc(d1*d1,sizeof(float));
  float ii;
  int area=d1;
  float r2=(float)d1/2;
  float ax=r2, bx=r2;
  float dspace=0.00001;
  float ndata=(float)area/dspace;
  int sum=0;
  int arr=area*ndata;

  MASK_CIRCLE cir,qa,qb;

  cir.dim=ndata*2;
  cir.x = (float*) calloc((cir.dim),sizeof(float));
  cir.y = (float*) calloc((cir.dim),sizeof(float));

  //printf("area:%d,r2:%f,ax:%f,bx:%f,dspace:%f,arr:,cir.dim:%d, ndata:%f\n",area,r2,ax,bx,dspace,arr,cir.dim,ndata);

  int pa=0,pb=0;

  for (i=0; i<ndata; i++)
  {
        float px=(i*dspace);
        cir.x[i]= px;
        cir.x[(int)(2*ndata)-1 - i]= px;

        cir.y[i] = (bx+sqrt((r2*r2)-((px-ax)*(px-ax))));
        cir.y[(int)(2*ndata)-1 - i]= (bx-sqrt((r2*r2)-((px-ax)*(px-ax))));
  }


  for(i=0; i<cir.dim ;i++)
  {
    if(isnan(cir.y[i]) == 1)
    {
        cir.x[i] = 0.0;
        cir.y[i] = 0.0;
    }
  }

  qa=supr_zeros(cir.x,cir.y,cir.dim);
  free_MASK_CIRCLE(cir);
  delete_equals(qa);
  cir=supr_zeros(qa.x,qa.y,qa.dim);
  free_MASK_CIRCLE(qa);

  maskC=map_area(cir,dspace,area,d1);
  free_MASK_CIRCLE(cir);
/*
for(j=0; j<d1 ; j++)
{
    for(i=0; i<d1 ; i++)
    {
        printf("%f ",maskC[(j*d1)+i]);
    }
    printf("\n");
}
printf("\n");*/

return maskC;
}

int * create_mask_ring(int d1, int d2)
{
  float r1 = (d1+0.2)/2;
  float r2 = (d2+0.2)/2;
  long int i,j;
  int *mask = (int*) calloc(d1*d1,sizeof(int));

  float  a,b,y1,y2, y1_b,y2_b;
  //printf("d1= %d , d2=%d aqui \n",d1,d2);
  a=d1/2;
  b=a;

  for (i=0; i<d1; i++)
  {
  y2= (b+sqrt((r1*r1)-((i-a)*(i-a)))) ;
  y1= (b-sqrt((r1*r1)-((i-a)*(i-a)))) ;

  y2_b= (b+sqrt((r2*r2)-((i-a)*(i-a)))) ;
  y1_b= (b-sqrt((r2*r2)-((i-a)*(i-a)))) ;
    for (j=0; j<d1; j++)
    {
        if(j <= (y2) && j >= (y1) )
        {
        mask[(i*d1)+j]=1;
        }
        else
        {
        mask[(i*d1)+j]=0;
        }
        if(j <= (y2_b) && j >= (y1_b) )
        {
        mask[(i*d1)+j]=0;
        }
    }
  }

  /*for (i=0; i<d1; i++)
  {
    for (j=0; j<d1; j++)
    {
        printf("%d ",mask[(i*d1)+j]);
    }
     printf("\n");
  }printf("\n");*/
return mask;
}

void free_STAR_BIN(STAR_BIN datos)
{
free(datos.bin);
free(datos.x);
free(datos.y);
free(datos.xe);
free(datos.ye);
free(datos.x0);
free(datos.y0);
free(datos.dt);
}

STAR_BIN select_stars(STAR_BIN comp, long int * naxes , float *night)
{
int j,i,k=0,l;
long int x[comp.cant], y[comp.cant], v[comp.cant], pos = 0;
float max, val[comp.cant];
for(j=0 ; j< naxes[1] ; j++)
{
    for(i=0 ; i< naxes[0] ; i++)
    {
       if(comp.bin[(j*naxes[0])+i] == 1)
       {
          x[pos] = i;
          y[pos] = j;
          v[pos] = 0;
          pos ++;
       }
    }
}

for (i=0 ; i<comp.cant ; i++)
{
   val[i] = 0.0;
   val[i] = get_max(comp.bin,&x[i],&y[i],night, naxes);
   if(i>0 && sqrt(pow((x[i]-x[i-1]),2)+pow((y[i]-y[i-1]),2)) > 25 && x[i] > 25  &&  x[i] < naxes[0] -25 && y[i] > 25  &&  y[i] < naxes[1] -25)
   {
   k++;
   }
   v[i] = k;
}
comp.x = (long int*) calloc(k+1,sizeof(long int));
comp.y = (long int*) calloc(k+1,sizeof(long int));

j=0;
l=0;
for(i=0 ; i<comp.cant ; i++)
{
    if(v[i] == j+1 )
    {
    comp.x[j] = (int)comp.x[j]/((i-l));
    comp.y[j] = (int)comp.y[j]/(i-l);
    l= i;
    j++;
    }

    comp.x[j] += x[i];
    comp.y[j] += y[i];
}
    comp.x[j] = (int)comp.x[j]/((i-l));
    comp.y[j] = (int)comp.y[j]/(i-l);

comp.cant = k+1;
return comp;
}

float get_max(long int *bin, long int *px, long int *py, float *night, long int * naxes)
{
int i,j, k = 3;

int x = px[0]-1;
int y = py[0]-1;

float max = 0.0;

for(j=y ; j< y+k ; j++)
{
    for(i=x ; i< x+k ; i++)
    {
        if(night[(j*naxes[0])+i] > max)
        {
        max = night[(j*naxes[0])+i];
        px[0] = i;
        py[0] = j;
        }
    }
}
return max;
}

void sort_two_arr(float *X, float *Y, int ndata, int d)
{
int i,j;
if (d==0)
{
for (i = 1; i < ndata; ++i)
  {
    float key = Y[i];
    float kex = X[i];
    j = i - 1;

    while (j >= 0 && Y[j] > key)
    {
        Y[j + 1] = Y[j];
        X[j + 1] = X[j];
        j = j - 1;
    }
    Y[j + 1] = key;
    X[j + 1] = kex;
  }
}
else if(d==1)
{
for (i = 1; i < ndata; ++i)
  {
    float key = Y[i];
    float kex = X[i];
    j = i - 1;

    while (j >= 0 && Y[j] < key)
    {
        Y[j + 1] = Y[j];
        X[j + 1] = X[j];
        j = j - 1;
    }
    Y[j + 1] = key;
    X[j + 1] = kex;
  }
}
return;
}

MASK_CIRCLE supr_zeros(float *x, float *y, int din)
{
int i;
int d=0;
MASK_CIRCLE out;

for (i=0;i<din;i++)
{
    if((x[i] + y[i])!=0.0)
    {
        d++;
    }
}
out.dim = d;
d=0;
//printf("OUT=%d\n",out.dim);
out.x = (float*) calloc((out.dim),sizeof(float));
out.y = (float*) calloc((out.dim),sizeof(float));
for (i=0;i<din;i++)
{
    if((x[i] + y[i])!=0.0)
    {
        out.x[d]=x[i];
        out.y[d]=y[i];
        d++;
    }
}
return out;
}

void free_MASK_CIRCLE(MASK_CIRCLE data)
{
free(data.x);
free(data.y);
return;
}

MASK_CIRCLE concat_arr(MASK_CIRCLE c1, MASK_CIRCLE c2, MASK_CIRCLE c3, MASK_CIRCLE c4)
{
int i;
int d=0;
MASK_CIRCLE out;
out.dim = c1.dim+c2.dim+c3.dim+c4.dim;

out.x = (float*) calloc((out.dim),sizeof(float));
out.y = (float*) calloc((out.dim),sizeof(float));

for (i=0;i<c1.dim;i++)
{
   out.x[i] = c1.x[i];
   out.y[i] = c1.y[i];
}
for (i=0;i<c3.dim;i++)
{
   out.x[i+c1.dim] = c3.x[i];
   out.y[i+c1.dim] = c3.y[i];
}
for (i=0;i<c4.dim;i++)
{
   out.x[i+c1.dim+c3.dim] = c4.x[i];
   out.y[i+c1.dim+c3.dim] = c4.y[i];
}
for (i=0;i<c2.dim;i++)
{
   out.x[i+c1.dim+c3.dim+c4.dim] = c2.x[i];
   out.y[i+c1.dim+c3.dim+c4.dim] = c2.y[i];
}

return out;
}

void delete_equals(MASK_CIRCLE data)
{
int i;
int d=0;
for(i=0; i<data.dim -1 ; i++)
{
    if (data.x[i] == data.x[i+1] && data.y[i] == data.y[i+1] )
    {
     data.x[i]=0.0;
     data.y[i]=0.0;
    }
}

if(data.x[i] == data.x[data.dim - 1])
{
    data.x[i] = 0.0;
    data.y[i] = 0.0;
}

return;
}

float * map_area(MASK_CIRCLE in, float dspace, int area, int d)
{
int i,j,k;
MASK_CIRCLE out;
out.dim = area;
float sp=(float)(area/2.0);
out.x = (float*) calloc((out.dim*out.dim),sizeof(float));
out.y = (float*) calloc((out.dim*out.dim),sizeof(float));
float * mask=(float*) calloc((out.dim*out.dim),sizeof(float));
float * ones=(float*) calloc((out.dim*out.dim),sizeof(float));
float * maskC=(float*) calloc((out.dim*out.dim),sizeof(float));

int p;
for(j=0; j<out.dim ; j++)
{
    for(i=0; i<out.dim ; i++)
    {
        p=0;
        for(k=0; k<in.dim;k++)
        {
            if (in.x[k] >= i && in.x[k] < (i+1))
            {
                if (in.y[k] > j && in.y[k] < (j+1))
                {
                    out.x[(j*out.dim )+i] += (in.y[k] - (float)j);
                    p++;
                }
            }
        }
        if(p!=0)
        {
            out.x[(j*out.dim )+i] = out.x[(j*out.dim )+i]/p;
            if(j< sp)
            {
                out.x[(j*out.dim )+i] = 1.0 - out.x[(j*out.dim )+i] ;
            }
            out.y[(i*out.dim )+j] = out.x[(j*out.dim )+i];
        }
        else
        {
            out.x[(j*out.dim )+i] = 1.0;
        }
    }
}

for(j=0; j<out.dim ; j++)
{
    for(i=0; i<out.dim ; i++)
    {
        if(out.x[(j*out.dim )+i] > out.y[(j*out.dim )+i])
        {
            mask[(j*out.dim)+i] = out.x[(j*out.dim)+i];
        }
        else if(out.x[(j*out.dim )+i] < out.y[(j*out.dim )+i])
        {
            mask[(j*out.dim)+i] =  out.y[(j*out.dim )+i];
        }
        else
        {
            mask[(j*out.dim)+i] =  out.x[(j*out.dim )+i];
        }
    }
}

for(j=0; j<out.dim ; j++)
{
    int dm1=0;
    for(i=0; i<(int)out.dim/2 ; i++)
    {
       if (mask[(j*out.dim)+i] < 1.0 && dm1==0)
            dm1=1;

       if(dm1==1)
       {
            ones[(j*out.dim)+i] = 1.0;
            ones[(j*out.dim)+(out.dim-1)-i] = 1.0;
       }

    }
}

for(j=0; j<out.dim ; j++)
{
    for(i=0; i<out.dim ; i++)
    {
        maskC[(j*out.dim)+i] = ((mask[(j*out.dim)+i]) * ones[(j*out.dim)+i]);
  //      printf("%f ",maskC[(j*out.dim)+i]);
    }
//    printf("\n");
}
//printf("\n");

free(mask);
free(ones);
free_MASK_CIRCLE(out);
return maskC;
}
