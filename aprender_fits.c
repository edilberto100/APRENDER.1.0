/**************************************************************
 *                           APPHI                            *
 * AN AUTOMATED PHOTOMETRY PIPELINE FOR HIGH CADENCE,         *
 *                     LARGE VOLUME DATA                      *
 * Edilberto Sanchez Moreno                                   *
 **************************************************************/
 #include "aprender.h"

int get_mef_size(char *filetxt)
{
    fitsfile *fptr;
    int status = 0, single = 0, hdupos,nkeys, mef=0;
    if (!fits_open_file(&fptr, filetxt, READONLY, &status))
    {
    fits_get_hdu_num(fptr, &hdupos);
    if (hdupos != 1 || strchr(filetxt, '[')) single = 1;
    for (; !status; hdupos++)  /* Main loop through each extension */
    {
        fits_get_hdrspace(fptr, &nkeys, NULL, &status);      /* get # of keywords */
        if (single) break;  /* quit if only listing a single header */
        fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
    }
    fits_close_file(fptr, &status);
    mef = hdupos -1;
    }
return mef;
}

LIST_TXT ** get_file_mef(char *filetxt, int mef)
{
    LIST_TXT ** files;
    files = (LIST_TXT**) calloc(mef,sizeof(LIST_TXT*));
    int i,j;
    for (i=0 ; i<mef ; i++)
    {
        int nfiles = 0;
        files[i] = read_list(filetxt,&nfiles);
        for (j=0 ; j<nfiles ; j++)
        {
            if(mef>1)
            name_mef(files[i][j].name,i+1);
            else
            name_mef(files[i][j].name,i);
        }
    }
return files;
}

void name_mef(char name[], int i)
{
  //char name2[] = "[";
  //char name3[] = " ";
  //char name4[] = "]";
  //sprintf(name3,"%d", i);
  //strcat(name, name2);
  //strcat(name, name3);
  //strcat(name, name4);
  //printf("%s\n",name);
  return;
}

 void save_fits_ld(FILE_FITS imagen, char name[])
{
  //** Crea el Archivo *.fits **//
  fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
  int status;
  status = 0;         /* initialize status before calling fitsio routines */
  fits_create_file(&fptr, name, &status);   /* create new file */

  /* Create the primary array image (16-bit short integer pixels */
  imagen.naxes[2]=imagen.naxes[2];
  fits_create_img(fptr, LONG_IMG, imagen.dim, imagen.naxes, &status);

  fits_write_img(fptr, TLONG, 1, imagen.naxes[0]*imagen.naxes[1]*(imagen.naxes[2]), imagen.data, &status);

  fits_close_file(fptr, &status);            /* close the file */
  fits_report_error(stderr, status);  /* print out any error messages */
}

 void save_fits_fl(FILE_FITS imagen, char name[])
{
  //** Crea el Archivo *.fits **//
  fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
  int status;
  status = 0;         /* initialize status before calling fitsio routines */
  fits_create_file(&fptr, name, &status);   /* create new file */

  /* Create the primary array image (16-bit short integer pixels */
  imagen.naxes[2]=imagen.naxes[2];
  fits_create_img(fptr, FLOAT_IMG, imagen.dim, imagen.naxes, &status);

  fits_write_img(fptr, FLOAT_IMG, 1, imagen.naxes[0]*imagen.naxes[1]*(imagen.naxes[2]), imagen.data, &status);

  fits_close_file(fptr, &status);            /* close the file */
  fits_report_error(stderr, status);  /* print out any error messages */
}

 LIST_TXT *read_list(char *filetxt, int *nlines)
{
	FILE *file;
	file = fopen(filetxt,"r");
	int i=0, large;
	nlines[0]= get_nlines(filetxt) -1;

	LIST_TXT *list;
	list= (LIST_TXT*) malloc(nlines[0]*sizeof(LIST_TXT));
    file = fopen(filetxt,"r");
	if (file == NULL){
		printf("\nError de apertura del archivo. \n\n");
        }
    else{
	    while (i < nlines[0])
	      {
	      fgets(list[i].name,500,file);
	      large = strlen(list[i].name);
	      if (list[i].name[large-1]== 's')
	      {
	        list[i].name[large-1]= 's';
	      }
	      else
	      {
          list[i].name[large-1]= '\0';
          }
          //printf("File:%s.\n",list[i].name);
          i++;
          }
        }
        fclose(file);
return list;
}

long int * get_diameter(char * data, int * dim)
{
    int i=0;
    int c1=0, count=0;
    int large = strlen(data);
    for(i=0; i<large; i++)
    {
        //printf("%c\n",data[i]);
        if(data[i] == ',' || data[i]==']')
        {
            count ++;
        }
    }
    //printf("datos = %d\n",count);
    dim[0]=count;
    long int * D1=(long int*) calloc(count,sizeof(long int));
    int l=1;
    int c=0;
    char dat[5];
    while(l<large)
    {
        dat[c1]=data[l];
        c1++;
        if(data[l] == ',' || data[l]==']')
        {
            dat[c1-1]='\0';
            D1[c]=atoi(dat);
            //printf("%ld\n",D1[c]);
            c1=0;
            c++;
        }
        l++;
    }
 return D1;
}

 STAR_BIN read_lst(char *filetxt, int *nlines)
{
	FILE *file;
    STAR_BIN st;
	file = fopen(filetxt,"r");
	int i=0,j, large;
	nlines[0]= get_nlines(filetxt) -1;

	LIST_TXT *list;
	list= (LIST_TXT*) malloc(nlines[0]*sizeof(LIST_TXT));
    st.x   = (long int*) calloc(nlines[0],sizeof(long int));
    st.y   = (long int*) calloc(nlines[0],sizeof(long int));
    st.xe  = (long int*) calloc(nlines[0],sizeof(long int));
    st.ye  = (long int*) calloc(nlines[0],sizeof(long int));
    st.x0  = (long int*) calloc(nlines[0],sizeof(long int));
    st.y0  = (long int*) calloc(nlines[0],sizeof(long int));
    st.bin = (long int*) calloc(nlines[0],sizeof(long int));
    st.dt = (double*) calloc(nlines[0],sizeof(double));
	st.cant = nlines[0];
    file = fopen(filetxt,"r");
	if (file == NULL)
        {
		printf("\nError de apertura del archivo. \n\n");
        }
        else
        {
	    while (i < nlines[0])
	    {
            int l=0;
            char A[10];
            char B[10];
            char C[10];
            char D[10];
            int a,b;
            int k0=0, k1=0, k2=0, k3=0;
            fgets(list[i].name,500,file);
            large = strlen(list[i].name);
            int cou=0;
            for(j=0;j<large;j++)
            {
                if(list[i].name[j]== ',')
                cou++;
            }
            if(cou == 1)
            {
                for(j=0;j<large;j++)
                {
                    if(list[i].name[j]== ',')
                        l++;

                    if(l==0 && list[i].name[j] != ',')
                    {
                        A[k0]=0;
                        A[k0]=list[i].name[j];
                        A[k0+1]='\0';
                        k0++;
                    }

                    if(l==1 && list[i].name[j] != ',' && list[i].name[j] != '\0')
                    {
                        B[k1]=0;
                        B[k1]=list[i].name[j];
                        B[k1+1]='\0';
                        k1++;
                    }
                }
            st.x[i]=atoi(A);
            st.y[i]=atoi(B);
            st.xe[i]=atoi(A);
            st.ye[i]=atoi(B);

            }

            if (cou == 3)
            {
                for(j=0;j<large;j++)
                {
                    if(list[i].name[j]== ',')
                        l++;

                    if(l==0 && list[i].name[j] != ',')
                    {
                        A[k0]=0;
                        A[k0]=list[i].name[j];
                        A[k0+1]='\0';
                        k0++;
                    }

                    if(l==1 && list[i].name[j] != ',' && list[i].name[j] != '\0')
                    {
                        B[k1]=0;
                        B[k1]=list[i].name[j];
                        B[k1+1]='\0';
                        k1++;
                    }

                    if(l==2 && list[i].name[j] != ',' && list[i].name[j] != '\0')
                    {
                        C[k2]=0;
                        C[k2]=list[i].name[j];
                        C[k2+1]='\0';
                        k2++;
                    }

                    if(l==3 && list[i].name[j] != ',' && list[i].name[j] != '\0')
                    {
                        D[k3]=0;
                        D[k3]=list[i].name[j];
                        D[k3+1]='\0';
                        k3++;
                    }
                }
            st.x[i]=atoi(A);
            st.y[i]=atoi(B);
            st.xe[i]=atoi(C);
            st.ye[i]=atoi(D);
            st.dt[i]=sqrt((st.xe[0]-st.x0[0])*(st.xe[0]-st.x0[0])+(st.ye[0]-st.y0[0])*(st.ye[0]-st.y0[0]));
            printf("dt=%lf\n",st.dt[i]);
            st.x0[i]=atoi(A);
            st.y0[i]=atoi(B);
            }


	      st.bin[i]=0;
          i++;
          }
        }
        for(i=0;i<st.cant;i++)
            printf("Star %d = [%ld,%ld] [%ld,%ld]\n",i,st.x[i],st.y[i],st.xe[i],st.ye[i]);
        fclose(file);
	free(list);
return st;
}

FILE_FITS *read_files(char *filetxt, int *nlines)
{
	FILE *file;
	file = fopen(filetxt,"r");
	int i=0, large;
	nlines[0]= get_nlines(filetxt) -1;

	LIST_TXT *list;
	list= (LIST_TXT*) malloc(nlines[0]*sizeof(LIST_TXT));
    file = fopen(filetxt,"r");
	if (file == NULL){
		printf("\nError de apertura del archivo. \n\n");
        }
    else{
	    while (i < nlines[0])
	      {
	      fgets(list[i].name,500,file);
	      large = strlen(list[i].name);
	      if (list[i].name[large-1]== 's')
	      {
	        list[i].name[large-1]= 's';
	      }
	      else
	      {
          list[i].name[large-1]= '\0';
          }
          i++;
          }
        }
        fclose(file);
     //printf("Entró\n");
     FILE_FITS * data = get_data_fits(list, nlines);
     //printf("Salió\n");
     free(list);


return data;
}

int get_nlines(char *filetxt)
{
FILE *file;
file = fopen(filetxt,"r");
int nline = 0;
char chars[500];
    if (file == NULL){
		printf("\nError de apertura del archivo. \n\n");
        }
    else{
	    while (feof(file) == 0)
	      {
	      fgets(chars,500,file);
          nline ++;
          }
        }
        fclose(file);
return nline;
}

FILE_FITS * get_data_fits(LIST_TXT *filesname, int * nfiles)
{
FILE_FITS * list= (FILE_FITS*) calloc(nfiles[0],sizeof(FILE_FITS));
long long int i,j,k;
int *V = (int*) calloc(nfiles[0],sizeof(int));
long long int P = 0;

for (i=0 ; i< nfiles[0] ; i++)
{
    fitsfile *fptr;
    int status=0;
    short nulval = -1;
    list[i].naxes[0]=0;
    list[i].naxes[1]=0;
    list[i].naxes[2]=0;
    list[i].dim = 0 ;
    fits_open_file(&fptr, filesname[i].name, READONLY, &status);
    fits_get_img_dim(fptr, &list[i].dim ,&status);
    fits_get_img_size(fptr,list[i].dim,list[i].naxes,&status);
    //printf("%d , %s , dim=%d , X=%ld, Y=%ld, Z=%ld\n",i, filesname[i].name, list[i].dim, list[i].naxes[0], list[i].naxes[1], list[i].naxes[2]);
    if (list[i].dim > 0)
    {
    V[i] = 1;
    P++;
    }
    else
    {
    V[i] = 0;
    }
    fits_close_file(fptr, &status);
}
free_FILE_FITS(list,nfiles[0]);
list= (FILE_FITS*) calloc(P,sizeof(FILE_FITS));
P = 0;
for (i=0 ; i< nfiles[0] ; i++)
{
//printf("V[%lld] = %d\n", i,V[i]);
    if(V[i] == 1)
    {
    fitsfile *fptr;
    int status=0;
    short nulval = -1;
    list[P].naxes[0]=0;
    list[P].naxes[1]=0;
    list[P].naxes[2]=0;
    list[P].dim = 0 ;
    fits_open_file(&fptr, filesname[i].name, READONLY, &status);
    fits_get_img_dim(fptr, &list[P].dim ,&status);
    fits_get_img_size(fptr,list[P].dim,list[P].naxes,&status);
    //printf("%s , dim=%d , X=%ld , Y=%ld, Z=%ld \n",filesname[i].name, list[P].dim, list[P].naxes[0], list[P].naxes[1], list[P].naxes[2]);
        if(list[P].dim == 2)
        {
        list[P].naxes[2]=1;
        }
        if(list[P].naxes[2]==0)
        {
            list[P].naxes[2]=1;
        }

    list[P].data = (long int*) calloc(list[P].naxes[0]*list[P].naxes[1]*list[P].naxes[2],sizeof(long int));
    fits_read_img(fptr, TULONG, 1, list[P].naxes[0]*list[P].naxes[1]*list[P].naxes[2], &nulval, list[P].data, 0 , &status);
    fits_close_file(fptr, &status);
    P++;
    }
}

nfiles[0] = P;
free(V);
return list;
}

void free_FILE_FITS(FILE_FITS * data, int N)
{
int i;
for(i=0; i<N; i++)
{
    free(data[i].data);
}
free(data);
return;
}

void file_name(char name[], char star[])
{
int i=0;
  while (name[i] > 46)
  {
  star[i]=name[i];
  i++;
  }
  star[i]='_';
  star[i+1]='\0';
return ;
}

double get_JD(char FileFits[], long int m)
{
fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
char jdtxt[256];
char jddat[256];
double JD=-1.0;
int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
int single = 0, hdupos, nkeys, ii;
double pos = -1.0;
//printf("%s \n",card);

if (!fits_open_file(&fptr, FileFits, READONLY, &status))
    {
    fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */
      /* List only a single header if a specific extension was given */
    if (hdupos != 1 || strchr(FileFits, '[')) single = 1;

    for (; !status; hdupos++)  /* Main loop through each extension */
    {
        fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* get # of keywords */
        //printf("Header listing for HDU #%d:\n", hdupos);
        for (ii = 1; ii <= nkeys; ii++) { /* Read and print each keywords */
           if (fits_read_record(fptr, ii, card, &status))break;
           if (card[0] == 74 && card[1] == 68  && card[2] == 32)  // 74=J, 68=D, 32=espacio
           {
               sprintf(jdtxt,"%s", card );
               pos = 0.0 ;
           }
        }
        //printf("END\n\n");  /* terminate listing with END */
        if (single) break;  /* quit if only listing a single header */
        fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
      }
      if (status == END_OF_FILE)  status = 0; /* Reset after normal error */
      fits_close_file(fptr, &status);
    }

    if (status) fits_report_error(stderr, status); /* print any error message */

    int large , c;

    if(pos==0.0){
        large = strlen(jdtxt);
        c=0;
        for (ii= 0 ; ii<large ; ii++)
        {
            if (jdtxt[ii] > 47 && jdtxt[ii]< 58 || jdtxt[ii] =='.')
            {
                jddat[c]=jdtxt[ii];
                c++;
            }
        }
        jddat[c]='\0';
        JD= atof(jddat);
        // Obteniendo el tiempo de integración
        double EXPTIME = get_EXPTIME(FileFits,m);
        double Y = ((EXPTIME/2.0) / (24*3600.0));
        JD = JD + Y;

        //printf("ds9 %s, JD=%lf\n",FileFits,JD);
    }
    else{
        JD=(double)m;
    }



	      //if (list[i].name[large-1]== 's')
return(JD);
}

double get_EXPTIME(char FileFits[], long int m)
{
fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
char jdtxt[256];
char jddat[256];
double JD=-1.0;
int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
int single = 0, hdupos, nkeys, ii;
double pos = -1.0;
//printf("%s \n",card);

if (!fits_open_file(&fptr, FileFits, READONLY, &status))
    {
    fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */
      /* List only a single header if a specific extension was given */
    if (hdupos != 1 || strchr(FileFits, '[')) single = 1;

    for (; !status; hdupos++)  /* Main loop through each extension */
    {
        fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* get # of keywords */
        //printf("Header listing for HDU #%d:\n", hdupos);
        for (ii = 1; ii <= nkeys; ii++) { /* Read and print each keywords */
           if (fits_read_record(fptr, ii, card, &status))break;
           if (card[0] == 'E' && card[1] == 'X'  && card[2] == 'P' && card[3] == 'T' && card[4] == 'I'  && card[5] == 'M')  // 74=J, 68=D, 32=espacio
           {
               sprintf(jdtxt,"%s", card );
               pos = 0.0 ;
           }
        }
        //printf("END\n\n");  /* terminate listing with END */
        if (single) break;  /* quit if only listing a single header */
        fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
      }
      if (status == END_OF_FILE)  status = 0; /* Reset after normal error */
      fits_close_file(fptr, &status);
    }

    if (status) fits_report_error(stderr, status); /* print any error message */

    int large , c;

    if(pos==0.0){
        large = strlen(jdtxt);
        c=0;
        for (ii= 0 ; ii<large ; ii++)
        {
            if (jdtxt[ii] > 47 && jdtxt[ii]< 58 || jdtxt[ii] =='.')
            {
                jddat[c]=jdtxt[ii];
                c++;
            }
        }
        jddat[c]='\0';
        JD= atof(jddat);
    }
    else{
        JD=(double)m;
    }


	      //if (list[i].name[large-1]== 's')
return(JD);
}

HEADER get_Header(char FileFits[], long int m)
{
fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
char jdtxt[1024];
char jddat[1024];
int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
int single = 0, hdupos, nkeys, ii;
double pos = -1.0;
HEADER header;
char parameter[1024];
//printf("%s \n",card);

if (!fits_open_file(&fptr, FileFits, READONLY, &status))
    {
    fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */
      /* List only a single header if a specific extension was given */
    if (hdupos != 1 || strchr(FileFits, '[')) single = 1;

    for (; !status; hdupos++)  /* Main loop through each extension */
    {
        fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* get # of keywords */
        printf("Header listing for HDU #%d:\n", hdupos);
        for (ii = 1; ii <= nkeys; ii++) { /* Read and print each keywords */
           if (fits_read_record(fptr, ii, card, &status))break;
              if (card[0] == 'O' && card[1] == 'B'  && card[2] == 'J' && card[3] == 'E' && card[4] == 'C'  && card[5] == 'T')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.object,"OBJECTNAME=");
                supr_text(header.object,parameter);
              }
              else if(card[0] == 'D' && card[1] == 'A'  && card[2] == 'T' && card[3] == 'E' && card[4] == '-'&& card[5] == 'O')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.DateObs,"SESSIONDATE=");
                supr_text(header.DateObs,parameter);
              }
              else if(card[0] == 'U' && card[1] == 'T'  && card[2] == ' ' && card[3] == ' ' && card[4] == ' ' && card[5] == ' ')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.UT,"SESSIONTIME=");
                supr_text(header.UT,parameter);
              }
              else if(card[0] == 'O' && card[1] == 'B'  && card[2] == 'S' && card[3] == 'E' && card[4] == 'R' && card[5] == 'V' && card[6] == 'E')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.Observer,"CONTACTNAME=");
                supr_text(header.Observer,parameter);
              }
              else if(card[0] == 'O' && card[1] == 'R'  && card[2] == 'I' && card[3] == 'G' && card[4] == 'I' && card[5] == 'N')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.Origin,"FACILITY=");
                supr_text(header.Origin,parameter);
              }
              else if(card[0] == 'L' && card[1] == 'A'  && card[2] == 'T' && card[3] == 'I' && card[4] == 'T' && card[5] == 'U')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.Latitud,"OBSLATITUDE=");
                supr_text(header.Latitud,parameter);
              }
              else if(card[0] == 'L' && card[1] == 'O'  && card[2] == 'N' && card[3] == 'G' && card[4] == 'I' && card[5] == 'T')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.Longitud,"OBSLONGITUDE=");
                supr_text(header.Longitud,parameter);
              }
              else if(card[0] == 'T' && card[1] == 'E'  && card[2] == 'L' && card[3] == 'E' && card[4] == 'S' && card[5] == 'C')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.Telescop,"TELESCOPE=");
                supr_text(header.Telescop,parameter);
              }
              else if(card[0] == 'C' && card[1] == 'C'  && card[2] == 'D' && card[3] == 'T' && card[4] == 'Y' && card[5] == 'P')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.CCDType,"DETECTOR=");
                supr_text(header.CCDType,parameter);
              }
              else if(card[0] == 'R' && card[1] == 'A'  && card[2] == ' ' && card[3] == ' ' && card[4] == ' ' && card[5] == ' ')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.RA,"OBJECTRA=");
                supr_text(header.RA,parameter);
              }
              else if(card[0] == 'D' && card[1] == 'E'  && card[2] == 'C' && card[3] == ' ' && card[4] == ' ' && card[5] == ' ')  //
              {
                //char parameter[256];
                sprintf(parameter,"%s", card );
                //strcat(header.Dec,"OBJECTDEC=");
                supr_text(header.Dec,parameter);
              }
              else if(card[0] == 'F' && card[1] == 'I'  && card[2] == 'L' && card[3] == 'T' && card[4] == 'E' && card[5] == 'R')  //
              {
                //char parameter[1024];
                sprintf(parameter,"%s", card );
                //strcat(header.Filter,"FILTER=");
                supr_text(header.Filter,parameter);
              }
              else if(card[0] == 'E' && card[1] == 'Q'  && card[2] == 'U' && card[3] == 'I' && card[4] == 'N' && card[5] == 'O')  //
              {
                sprintf(jdtxt,"%s", card );
                get_value(header.Equinox,jdtxt);
              }
              else if ((card[0] == 'E' && card[1] == 'X'  && card[2] == 'P' && card[3] == 'T' && card[4] == 'I' && card[5] == 'M') || (card[0] == 'E' && card[1] == 'X'  && card[2] == 'P' && card[3] == 'O' && card[4] == 'S' && card[5] == 'U') ) //
              {
                sprintf(jdtxt,"%s", card );
                get_value(header.ExpTime,jdtxt);
                //sprintf(header.ExpTime,"%s", card );
              }

        }

        printf("%s \n",header.ExpTime);
        printf("%s \n",header.Equinox);

        printf("END\n\n");  /* terminate listing with END */
        if (single) break;  /* quit if only listing a single header */
        fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
      }
      if (status == END_OF_FILE)  status = 0; /* Reset after normal error */
      fits_close_file(fptr, &status);
    }
    if (status) fits_report_error(stderr, status); /* print any error message */

return(header);
}

void supr_text(char id[], char parameter[])
{
int i, large , c=0;
int pos=0;
large = strlen(parameter);
char newtext[1024];
memset(newtext, ' ', 1024);

for (i= 0 ; i<large ; i++)
{
    if(parameter[i]=='\''){
        pos++;
        i++;
    }

    if(pos>0 && pos<2){
        newtext[c]=parameter[i];
        c++;
    }
}
newtext[c]='\0';
//sprintf(id,"%s",newtext);
strcat(id, newtext);
//printf("%s\n",id);
return ;
}

void get_value(char id[], char parameter[])
{
int i, large , c=0;
int pos=0;
large = strlen(parameter);
char newtext[1024];
memset(newtext, ' ', 1024);

for (i= 0 ; i<large ; i++)
{
    if (parameter[i] > 47 && parameter[i]< 58 || parameter[i] =='.')
            {
                newtext[c]=parameter[i];
                c++;
            }
}
newtext[c]='\0';
if (newtext[c-1] < 48 || newtext[c-1] > 57)
    newtext[c-1]='\0';
//sprintf(id,"%s",newtext);
strcat(id, newtext);
//printf("%s\n",id);
return ;
}
