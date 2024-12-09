/**************************************************************
 *                           APPHI                            *
 * AN AUTOMATED PHOTOMETRY PIPELINE FOR HIGH CADENCE,         *
 *                     LARGE VOLUME DATA                      *
 * Edilberto Sanchez Moreno                                   *
 **************************************************************/
#include "aprender.h"

int apphi_field(struct arguments arguments) {
  /********************* Review of the parameters into *************************/
    clock_t start, stop;
    start = clock();
    printf("\n|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|\n");
    printf(" APRENDER has been initialized successfully \n");
    printf(" Starting to measure time \n");
  /*****************************************************************************/

  /************************* Data condition ************************************/
  FILE_FITS *bias_files, *flats_files, *night_files, *darks_files;
  STAR_BIN comp;
  float Umb;
  long int lz;
  Umb = 0.4;
  int nbias, ndarks, nflats, nnight;
  long int i,j,k;
  if (arguments.np < 1 || arguments.np > 99)
  lz = 1;
  else
  lz = arguments.np;
  char name[500];
  file_name(arguments.science,name);
  MASTER_DATA darks, flats, bias;
  /*****************************************************************************/

  /********************************* Get list **********************************/
  LIST_TXT *files;
  long int nx_z = 0, *nx;
  int mef_science, mef_bias, mef_darks, mef_flats;
  LIST_TXT **filemef_science;
  LIST_TXT **filemef_bias;
  LIST_TXT ** filemef_darks;
  LIST_TXT ** filemef_flats;
  int B =1, Da=1, F=1, m;
  FILE * file;

  file = fopen(arguments.science,"r");
  if (file == NULL)
  {
    printf(" ERROR; APRENDER don't have Science Image\n");
    exit(0);
  }
  else
  {
    fclose(file);
    files = read_list(arguments.science,&nnight);
    mef_science = get_mef_size(files[0].name);
    printf(" Image MEF = %d\n",mef_science);
    filemef_science=get_file_mef(arguments.science,mef_science);
    free(files);
  }
  file = fopen(arguments.bias,"r");
  if (file == NULL)
  {
    printf(" Warning; APRENDER don't have BIAS Image\n");
    B = 0;
  }
  else
  {
    fclose(file);
    files = read_list(arguments.bias,&nbias);
    mef_bias = get_mef_size(files[0].name);
    printf(" Bias MEF = %d\n",mef_bias);
    if (mef_bias == mef_science)
    filemef_bias=get_file_mef(arguments.bias,mef_bias);
    else
    B = 0;
    free(files);
  }
  file = fopen(arguments.dark,"r");
  if (file == NULL)
  {
    printf(" Warning; APRENDER don't have Darks Image\n");
    Da= 0;
  }
  else
  {
    fclose(file);
    files = read_list(arguments.dark,&ndarks);
    mef_darks = get_mef_size(files[0].name);
    printf(" Dark MEF=%d\n",mef_darks);
    if (mef_darks == mef_science)
        filemef_darks=get_file_mef(arguments.dark,mef_darks);
    else
        Da = 0;
    free(files);
  }
  file = fopen(arguments.flat,"r");
  if (file == NULL)
  {
    printf(" Warning; APRENDER don't have Flat Image\n");
    F = 0;
  }
  else
  {
    fclose(file);
    files = read_list(arguments.flat,&nflats);
    mef_flats = get_mef_size(files[0].name);
    printf(" Flat MEF = %d\n",mef_flats);
    if (mef_flats == mef_science)
        filemef_flats = get_file_mef(arguments.flat,mef_flats);
    else
        F = 0;
    free(files);
  }

  printf(" Header is reading\n");
  HEADER header=get_Header(filemef_science[0][0].name,0);

  for(m=0; m<mef_science; m++)
  {
  printf(" APRENDER is working \n");
  nx = (long int *) calloc(nnight,sizeof(long int));
  long int dimx=0, dimy=0;
  nx_z = 0;
  double *JD = (double *) calloc(nnight,sizeof(double));
  double *EXPTIME= (double *) calloc(nnight,sizeof(double));
  double *xtmp=(double *) calloc(nnight,sizeof(double));
  double *ytmp=(double *) calloc(nnight,sizeof(double));
  double TexpTotal=0.0;
  for (i=0; i<nnight; i++)
  {
    int nf=1;
    JD[i]= get_JD(filemef_science[m][i].name,i);
    EXPTIME[i]=get_EXPTIME(filemef_science[m][i].name,i);
    TexpTotal=TexpTotal+EXPTIME[i];
    night_files = get_data_fits(&filemef_science[m][i], &nf);
    nx[i] = night_files[0].naxes[2];
    dimx=night_files[0].naxes[0];
    dimy=night_files[0].naxes[1];
    nx_z += night_files[0].naxes[2];
    free_FILE_FITS(night_files,1);
  }

   if(nx_z > nnight)
   {
      for (i=0; i<nnight; i++)
        {
            JD[i]=1.0;
        }
   }

  /*****************************************************************************/
  printf(" Read Bias data\n");

  /***************************** Read data bias ********************************/
  if (B == 0)
  {
  bias.naxes[0] = dimx;
  bias.naxes[1] = dimy;
  bias.naxes[2] = 1;
  bias.data = (float*) calloc(bias.naxes[0]*bias.naxes[1],sizeof(float));
  printf(" bias ------------ [OK]       \n");
  }
  else
  {
  bias_files = get_data_fits(filemef_bias[m],&nbias);
  free(filemef_bias[m]);
  bias = master_bias(bias_files,nbias);
  free_FILE_FITS(bias_files,nbias);
  printf(" bias ------------ [OK]       \n\n");
  }
  /*****************************************************************************/

  printf(" Read Darks data\n");

  /****************************** Read data darks ******************************/
  if (Da == 0)
  {
  darks.naxes[0] = dimx;
  darks.naxes[1] = dimy;
  darks.naxes[2] = 1;
  darks.data = (float*) calloc(darks.naxes[0]*darks.naxes[1],sizeof(float));
  printf(" darks ----------- [OK]       \n\n");
  }
  else
  {
  darks_files = get_data_fits(filemef_darks[m],&ndarks);
  free(filemef_darks[m]);
  darks = master_darks(bias,darks_files,ndarks);
  free_FILE_FITS(darks_files,ndarks);
  printf("darks ----------- [OK]       \n\n");
  }
  /*****************************************************************************/
 printf(" Read data flats data\n");

  /**************************** Read data flats ********************************/
  if (F == 0)
  {
  flats.naxes[0] = dimx;
  flats.naxes[1] = dimy;
  flats.naxes[2] = 1;
  flats.data = (float*) calloc(flats.naxes[0]*flats.naxes[1],sizeof(float));
	for(i=0 ; i<(flats.naxes[0]*flats.naxes[1]) ; i++)
	flats.data[i] = 1.0;
  printf(" flats ----------- [OK]       \n");
  }
  else
  {
  flats_files = get_data_fits(filemef_flats[m],&nflats);
  free(filemef_flats[m]);
  flats = master_flats(bias, darks, flats_files, nflats);
  free_FILE_FITS(flats_files,nflats);
  printf(" flats ----------- [OK]       \n");
  }
  /*****************************************************************************/


  /************************* Select the stars **********************************/
  int nf=1, nist, nobj, naper;
  night_files = get_data_fits(&filemef_science[m][0], &nf);
  file = fopen(arguments.ist,"r");
  if (file == NULL)
  {
  printf(" The automatic search algorithm has been activated");
  comp = conv_img(night_files[0].data,Umb,night_files[0].naxes,15,bias.data);
  }
  else
  {
  fclose(file);
  comp = read_position(arguments,&nobj,&naper);
  }
  MASTER_DATA * night;
  if (comp.cant == 0)
  {
  free(files);
  free_MASTER_DATA(&bias,1);
  free_MASTER_DATA(&darks,1);
  free_MASTER_DATA(&flats,1);
  free(nx);
  printf("\n APRENDER found %ld stars  \n", comp.cant);
  exit(0);
  }
  else
  printf("\n APRENDER is working with %ld stars\n", comp.cant);
  /*****************************************************************************/

  /*****************************************************************************/
  STAR_DATA *st, st_aux;
  MASTER_DATA *sky, sky_aux;
  int *D, *D2;
  st = (STAR_DATA*) calloc(comp.cant,sizeof(STAR_DATA));
  sky = (MASTER_DATA*) calloc(comp.cant,sizeof(MASTER_DATA));
  D  = (int*) calloc(comp.cant,sizeof(int));
  D2 = (int*) calloc(comp.cant,sizeof(int));

  for(i=0 ; i<comp.cant ; i++)
  {

   D[i]=comp.bin[i];
   D2[i]=comp.D2[i];
   int D3=comp.D3[i];

   //printf(" star(%ld), [%ld,%ld]. Diameters:[%ld,%ld,%ld]\n",i,comp.x[i],comp.y[i],comp.bin[i],comp.D3[i], comp.D2[i]);
   star_name(name,comp.x[i],comp.y[i],D[i],m,st[i].name);
   txt_name(name,comp.x[i],comp.y[i],D[i],m,st[i].txt);
   st[i].SNR=(float*) calloc(lz,sizeof(float));
   st[i].naxes[0] = D[i];
   st[i].naxes[1] = D[i];
   st[i].naxes[2] = nx_z;
   st[i].xy    = (float*) calloc(nx_z,sizeof(float));
   st[i].px    = (float*) calloc(nx_z,sizeof(float));
   st[i].py    = (float*) calloc(nx_z,sizeof(float));
   st[i].xyabs = (float*) calloc(nx_z,sizeof(float));
   st[i].star  = (float*) calloc(nx_z*D[i]*D[i],sizeof(float));
   sky[i].naxes[0] = D2[i];
   sky[i].naxes[1] = D2[i];
   sky[i].naxes[2] = nx_z;
   sky[i].data = (float*) calloc(nx_z*D2[i]*D2[i],sizeof(float));
  }
  /*****************************************************************************/

  /*****************************************************************************/
  long int data = 0, l;
  for(l=0 ; l<nnight ; l++)
  {
    night_files = get_data_fits(&filemef_science[m][l], &nf);
    night = image_reduction(bias,darks,flats,night_files,1);
    //printf("Empezando la Fotometría %ld \n",l);
    for(i=0 ; i<comp.cant ; i++)
    {
        //printf("Estrella en proceso: %ld \n",i);
        sky_aux.naxes[0]=0;
        sky_aux.naxes[1]=0;
        sky_aux.naxes[2]=0;
        get_naxes(night,1,sky_aux.naxes);
        st_aux.xy    = (float*) calloc(sky_aux.naxes[2],sizeof(float));
        st_aux.px    = (float*) calloc(sky_aux.naxes[2],sizeof(float));
        st_aux.py    = (float*) calloc(sky_aux.naxes[2],sizeof(float));
        st_aux.xyabs = (float*) calloc(sky_aux.naxes[2],sizeof(float));
        st_aux.naxes[0] = st[i].naxes[0];
        st_aux.naxes[1] = st[i].naxes[1];
        st_aux.naxes[2] = 1;
        //printf("Paramentros de entrada: %d, %d, %d, %d, %ld, %ld, %ld\n",comp.x[i],comp.y[i],1,D[i],st_aux.naxes[0],st_aux.naxes[1],st_aux.naxes[2]);
        //double TexpPixel=comp.dt[i]/TexpTotal;

        if(comp.x0[i]-comp.xe[i] != 0)
        {
            double TexpPixel=TexpTotal/nnight;
            double pxx=((double)comp.xe[i]-(double)comp.x0[i])/(double)nnight;
            double pyy=((double)comp.ye[i]-(double)comp.y0[i])/(double)nnight;
            comp.x[i]=(long int) comp.x0[i]+(l*pxx);
            comp.y[i]=(long int) comp.y0[i]+(l*pyy);
            //printf("nnight=%d, TexpPix=%lf, Xo[%ld,%ld], Xe[%ld,%ld], pxx,pyy=[%lf,%lf], Recorido[%ld,%ld]  \n",nnight,TexpPixel,comp.x0[i],comp.y0[i],comp.xe[i],comp.ye[i],pxx,pyy,comp.x[i],comp.y[i]);
            //
            get_centroids_list(comp.x[i],comp.y[i],night,1,arguments.d2,st_aux.naxes, st_aux.xy, st_aux.px, st_aux.py, st_aux.xyabs,EXPTIME[l]); //Hice más grande el radio de busqueda
        }
        get_centroids_list(comp.x[i],comp.y[i],night,1,arguments.d2,st_aux.naxes, st_aux.xy, st_aux.px, st_aux.py, st_aux.xyabs,EXPTIME[l]); //Hice más grande el radio de busqueda

        st_aux.star  = get_st(night,1,D[i],st_aux.naxes,st_aux.px,st_aux.py);
        sky_aux.data = get_st(night,1,D2[i],sky_aux.naxes,st_aux.px,st_aux.py);

        /*if (l==0)
        {
            printf("Star:\n");
            long int fx,fy;
            for(fx=0;fx<D[i];fx++)
            {
                for(fy=0;fy<D[i];fy++)
                {
                    printf("%lf ",st_aux.star[(fx*D[i])+fy]);
                }
                printf("\n");
            }

            printf("SKY:\n");
            for(fx=0;fx<D2[i];fx++)
            {
                for(fy=0;fy<D2[i];fy++)
                {
                    printf("%lf ",sky_aux.data[(fx*D2[i])+fy]);
                }
                printf("\n");
            }
        }*/

        for(j=0 ; j<nx[l] ; j++)
        {
            st[i].xy[data + j]    = st_aux.xy[j];
            st[i].px[data + j] = st_aux.px[j];
            st[i].py[data + j]    = st_aux.py[j];
            st[i].xyabs[data + j]    = st_aux.xyabs[j];
            for(k=0; k<(D[i]*D[i]); k++)
                st[i].star[(D[i]*D[i]*(data+j))+k] = st_aux.star[(D[i]*D[i]*j)+k];
            for(k=0; k<(D2[i]*D2[i]); k++)
                sky[i].data[(D2[i]*D2[i]*(data+j))+k] = sky_aux.data[(D2[i]*D2[i]*j)+k];
        }

        free_MASTER_DATA(&sky_aux,1);
        free(st_aux.xy);
        free(st_aux.px);
        free(st_aux.py);
        free(st_aux.xyabs);
        free(st_aux.star);
    }
    //printf("Terminando la Fotometría \n");
    data += nx[l];
    free_MASTER_DATA(night,1);
    for (i=0 ; i<comp.cant ; i++)
    {
        comp.x[i] = st[i].px[data -1];
        comp.y[i] = st[i].py[data -1];
    }
  }
  /*****************************************************************************/


  /***************************************************************************************/
  for(i=0 ; i<comp.cant ; i++)
  {
   st[i].err=(float*) calloc(st[i].naxes[2],sizeof(float));
   st[i].errfluxes=(float*) calloc(st[i].naxes[2],sizeof(float));
   st[i].curve=get_mag(st[i].star,sky[i].data,st[i].naxes,sky[i].naxes,comp.D3[i],arguments.tinteger,st[i].err,st[i].errfluxes);
   st[i].fluxes=get_fluxes(st[i].star,sky[i].data,st[i].naxes,sky[i].naxes,comp.D3[i]);
   st[i].z= size_snr(st[i].naxes[2],lz);
   st[i].C_SNR = sign_noise(st[i].star,sky[i].data,st[i].naxes,D2[i],st[i].z,lz,comp.D3[i]);
   st[i].SNR[0]= sign_curve(st[i].curve,st[i].naxes[2]);
   for(j=1 ; j<lz ; j++)
   {
        st[i].SNR[j]= sign_curve(st[i].C_SNR[j-1],st[i].naxes[2]/st[i].z[j-1]);
   }
  }
  printf(" APRENDER saved the data who:\n");
  for(i=0 ; i<naper ; i++)
  {
     save_ALCDEF_complete(st,comp.cant,JD,header,nobj,i,name);
  }

  if(comp.cant == 1)
  {
    st[0].diff = get_diff_curve(st[0].curve,st[0].curve,st[0].naxes);
    save_hdf5(st[0],lz,JD);
    save_ALCDEF(st[0],lz,JD,header);
    free_STAR_DATA(st[0],lz);
  }
  else
  {
    st[0].diff = get_diff_curve(st[0].curve,st[1].curve,st[0].naxes);
    save_hdf5(st[0],lz,JD);
    save_ALCDEF(st[0],lz,JD,header);
    free_STAR_DATA(st[0],lz);
    for(i=1 ; i<comp.cant ; i++)
    {
        st[i].diff = get_diff_curve(st[i].curve,st[i-1].curve,st[i].naxes);
        save_hdf5(st[i],lz,JD);
        save_ALCDEF(st[i],lz,JD,header);
        free_STAR_DATA(st[i],lz);
    }
  }

  /***************************************************************************************/

  /***************************************************************************************/
  printf("\n|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|\n");
  printf(" APRENDER to finished\n");
  stop = clock();
  double segundos = (double)(stop - start) / CLOCKS_PER_SEC;
  printf(" APRENDER worked for %f Seg. \n",segundos);
  printf("|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|\n");
  /***************************************************************************************/


  /************************************/
  free(filemef_science[m]);
  free_MASTER_DATA(&bias,1);
  free_MASTER_DATA(&darks,1);
  free_MASTER_DATA(&flats,1);
  /************************************/

  /********* Free data Stars *******/
  free(st);
  free(sky);
  free(D);
  free(D2);
  free(nx);
  free(comp.D2);
  free(comp.D3);
  free_STAR_BIN(comp);
  free(JD);
  free(EXPTIME);
  free(xtmp);
  free(ytmp);
 /**********************************/

 }
 free(filemef_science);
  return 0;
}
