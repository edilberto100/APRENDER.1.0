/**************************************************************
 *                           APPHI                            *
 * AN AUTOMATED PHOTOMETRY PIPELINE FOR HIGH CADENCE,         *
 *                     LARGE VOLUME DATA                      *
 * Edilberto Sanchez Moreno                                   *
 **************************************************************/
 #include "aprender.h"

void save_image_fl(hid_t file_id, float *pix, int NOBS, int WSIZE_x, int WSIZE_y, int NTEL, char name[])
{
int NDIM = 4;
hid_t space_id, set_id;
herr_t hval;
hsize_t dim[NDIM] ;
dim[0] = NOBS;
dim[1] = NTEL;
dim[2] = WSIZE_x;
dim[3] = WSIZE_y;

/*create the dataspace's image */
space_id = H5Screate_simple (NDIM, dim, NULL);
if (space_id < 0)
 {
 fprintf (stderr, "H5Screate_simple() failed\n");
 exit (1);
 }
/* create the dataset's image */
set_id = H5Dcreate (file_id, name, H5T_IEEE_F32LE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
if (set_id < 0)
 {
 fprintf (stderr, "H5Screate() failed\n");
 exit (1);
 }
/* write image to data file */
hval = H5Dwrite (set_id, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pix);
if (hval < 0)
 {
 fprintf (stderr, "H5Dwrite() failed\n");
 exit (1);
 }
/* free resources */
hval = H5Sclose (space_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Sclose() failed\n");
 exit (1);
 }
hval = H5Dclose (set_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Dclose() failed\n");
 exit (1);
 }
return;
}

void save_image_do(hid_t file_id, double *pix, int NOBS, int WSIZE_x, int WSIZE_y, int NTEL, char name[])
{
int NDIM = 4;
hid_t space_id, set_id;
herr_t hval;
hsize_t dim[NDIM] ;
dim[0] = NOBS;
dim[1] = NTEL;
dim[2] = WSIZE_x;
dim[3] = WSIZE_y;

/*create the dataspace's image */
space_id = H5Screate_simple (NDIM, dim, NULL);
if (space_id < 0)
 {
 fprintf (stderr, "H5Screate_simple() failed\n");
 exit (1);
 }
/* create the dataset's image */
set_id = H5Dcreate (file_id, name, H5T_IEEE_F64LE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
if (set_id < 0)
 {
 fprintf (stderr, "H5Screate() failed\n");
 exit (1);
 }
/* write image to data file */
hval = H5Dwrite (set_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pix);
if (hval < 0)
 {
 fprintf (stderr, "H5Dwrite() failed\n");
 exit (1);
 }
/* free resources */
hval = H5Sclose (space_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Sclose() failed\n");
 exit (1);
 }
hval = H5Dclose (set_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Dclose() failed\n");
 exit (1);
 }
return;
}

void save_vec_f(hid_t file_id, hid_t grp_id, char name[], float *lightcurve, int NOBS, int NTEL)
{
hid_t space_id, set_id;
herr_t hval;

hsize_t pdim[2] ;
pdim[0] = NTEL;
pdim[1] = NOBS;


/* create the dataspace's group */
space_id = H5Screate_simple (2, pdim, NULL);
if (space_id < 0)
 {
 fprintf (stderr, "H5Screate_simple() failed\n");
 exit (1);
 }
 /* lightcurve */
set_id = H5Dcreate (grp_id, name, H5T_IEEE_F32LE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
if (set_id < 0)
 {
 fprintf (stderr, "H5Dcreate() failed\n");
 exit (1);
 }
hval = H5Dwrite (set_id, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lightcurve);
if (hval < 0)
 {
  fprintf (stderr, "H5Dwrite() failed\n");
  exit (1);
 }
hval = H5Dclose (set_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Dclose() failed\n");
 exit (1);
 }
hval = H5Sclose (space_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Sclose() failed\n");
 exit (1);
 }
}

void save_vec_lf(hid_t file_id, hid_t grp_id, char name[], double *lightcurve, int NOBS, int NTEL)
{
hid_t space_id, set_id;
herr_t hval;

hsize_t pdim[2] ;
pdim[0] = NTEL;
pdim[1] = NOBS;


/* create the dataspace's group */
space_id = H5Screate_simple (2, pdim, NULL);
if (space_id < 0)
 {
 fprintf (stderr, "H5Screate_simple() failed\n");
 exit (1);
 }
 /* lightcurve */
set_id = H5Dcreate (grp_id, name, H5T_IEEE_F64LE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
if (set_id < 0)
 {
 fprintf (stderr, "H5Dcreate() failed\n");
 exit (1);
 }
hval = H5Dwrite (set_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lightcurve);
if (hval < 0)
 {
  fprintf (stderr, "H5Dwrite() failed\n");
  exit (1);
 }
hval = H5Dclose (set_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Dclose() failed\n");
 exit (1);
 }
hval = H5Sclose (space_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Sclose() failed\n");
 exit (1);
 }
}

void save_vec_i(hid_t file_id, hid_t grp_id, char name[], int *lightcurve, int NOBS, int NTEL)
{
hid_t space_id, set_id;
herr_t hval;

hsize_t pdim[2] ;
pdim[0] = NTEL;
pdim[1] = NOBS;


/* create the dataspace's group */
space_id = H5Screate_simple (2, pdim, NULL);
if (space_id < 0)
 {
 fprintf (stderr, "H5Screate_simple() failed\n");
 exit (1);
 }
 /* lightcurve */
set_id = H5Dcreate (grp_id, name, H5T_STD_I32LE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
if (set_id < 0)
 {
 fprintf (stderr, "H5Dcreate() failed\n");
 exit (1);
 }
hval = H5Dwrite (set_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, lightcurve);
if (hval < 0)
 {
  fprintf (stderr, "H5Dwrite() failed\n");
  exit (1);
 }
hval = H5Dclose (set_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Dclose() failed\n");
 exit (1);
 }
hval = H5Sclose (space_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Sclose() failed\n");
 exit (1);
 }
}

void save_image_i(hid_t file_id, int *pix, int NOBS, int WSIZE_x, int WSIZE_y, int NTEL, char name[])
{
int NDIM = 4;
hid_t space_id, set_id;
herr_t hval;
hsize_t dim[NDIM] ;
dim[0] = NOBS;
dim[1] = NTEL;
dim[2] = WSIZE_x;
dim[3] = WSIZE_y;

/*create the dataspace's image */
space_id = H5Screate_simple (NDIM, dim, NULL);
if (space_id < 0)
 {
 fprintf (stderr, "H5Screate_simple() failed\n");
 exit (1);
 }
/* create the dataset's image */
set_id = H5Dcreate (file_id, name, H5T_STD_I32LE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
if (set_id < 0)
 {
 fprintf (stderr, "H5Screate() failed\n");
 exit (1);
 }
/* write image to data file */
hval = H5Dwrite (set_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pix);
if (hval < 0)
 {
 fprintf (stderr, "H5Dwrite() failed\n");
 exit (1);
 }
/* free resources */
hval = H5Sclose (space_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Sclose() failed\n");
 exit (1);
 }
hval = H5Dclose (set_id);
if (hval < 0)
 {
 fprintf (stderr, "H5Dclose() failed\n");
 exit (1);
 }
return;
}

void get_vec_f(char * file)
{
    hid_t    file_id;
    hid_t    grp;
    herr_t   status;
    /*
     *  Example: open a file, open the root, scan the whole file.
     */
    file_id = H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT);
    grp = H5Gopen2(file_id, "/", H5P_DEFAULT);
    scan_group(grp);

    status = H5Fclose(file_id);
    status = H5Gclose (grp);
}

void scan_group(hid_t group)
{
    hid_t   grpid, dsid;
    herr_t  err;
    hsize_t  nobj;
    int MAX_NAME = 1024, i;
    int otype;
    char group_name[MAX_NAME];
    char memb_name[MAX_NAME];
    ssize_t  len;

    len = H5Iget_name(group, group_name, MAX_NAME);
    err = H5Gget_num_objs(group, &nobj);
    printf("Group Name: %s\n",group_name);

    for (i = 0; i < nobj; i++)
    {
		fflush(stdout);
		len = H5Gget_objname_by_idx(group, (hsize_t)i,memb_name, (size_t)MAX_NAME );
		fflush(stdout);
		fflush(stdout);
		otype =  H5Gget_objtype_by_idx(group, (size_t)i );

		switch(otype) {
			case H5G_GROUP:
				grpid = H5Gopen2(group,memb_name,H5P_DEFAULT);
				scan_group(grpid);
				H5Gclose(grpid);
				break;
			case H5G_DATASET:
				dsid = H5Dopen2(group,memb_name,H5P_DEFAULT);
				do_dset(dsid);
				H5Dclose(dsid);
				break;
			default:
				printf(" unknown?\n");
				break;
    }
   }

}

void do_dset(hid_t did)
{
    int MAX_NAME = 1024;
	hid_t tid;
	hid_t sid;
	hid_t dataset;
    hid_t memspace;
	hsize_t size;
	char ds_name[MAX_NAME];
    herr_t   status;
    int rank,status_n , i;
    HDF5_DATA_FL image;

	H5Iget_name(did, ds_name, MAX_NAME  );

	dataset = H5Dopen2(did, ds_name, H5P_DEFAULT);

	sid = H5Dget_space(dataset);       /* the dimensions of the dataset (not shown) */
    image.Rank      = H5Sget_simple_extent_ndims(sid);
    hsize_t     dims_out[image.Rank];
    status_n  = H5Sget_simple_extent_dims(sid, dims_out, NULL);
    long int tam = dims_out[0];

    printf("rank %d, dimensions %ld ", image.Rank,(long)(dims_out[0]));
    for (i=1 ; i<image.Rank ; i++)
    {
    tam = tam*dims_out[i];
    printf("x %ld ", (long)(dims_out[i]));
    }
    printf("\n");

    tid = H5Dget_type(did);
    H5T_class_t t_class;
	t_class = H5Tget_class(tid);

    memspace = H5Screate_simple(image.Rank,dims_out,NULL);
	size = H5Dget_storage_size(did);

    float  * data = (float*)calloc(tam,sizeof(float));
    status = H5Dread(dataset, tid, memspace, sid,H5P_DEFAULT, data);

    free(data);
    H5Sclose(memspace);
	H5Dclose(dataset);
    H5Sclose(sid);
    H5Tclose(tid);
}

float * get_img_f(char * file , char * section)
{
    hid_t file_id, grp_id1;
    file_id = H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT);
    grp_id1 = H5Gopen2(file_id, section, H5P_DEFAULT);
    hid_t    dataset, dataspace;
    herr_t   status;

    hid_t       datatype;
    H5T_class_t t_class;
    H5T_order_t order;

    size_t      size;
    int rank,status_n , i;
    hid_t       memspace;

    dataset = H5Dopen2(grp_id1, "STAR", H5P_DEFAULT);
    if( dataset < 0)
    printf(" Dataset \"STAR\" is not found. \n");
    printf("\"STAR\" dataset is open \n");

    datatype  = H5Dget_type(dataset);
    t_class     = H5Tget_class(datatype);  //---- Tipo de variable
    size  = H5Tget_size(datatype);
    printf(" Data size is %d \n", (int)size);

    dataspace = H5Dget_space(dataset);    /* dataspace handle */
    rank      = H5Sget_simple_extent_ndims(dataspace);
    hsize_t     dims_out[rank];
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    long int tam = dims_out[0];

    printf("rank %d, dimensions %ld ", rank,(long)(dims_out[0]));
    for (i=1 ; i<rank ; i++)
    {
    tam = tam*dims_out[i];
    printf("x %ld ", (long)(dims_out[i]));
    }
    printf("\n");


    memspace = H5Screate_simple(rank,dims_out,NULL);
    float * data = (float *) calloc(tam,sizeof(float));
    status = H5Dread(dataset, H5T_IEEE_F32LE, memspace, dataspace,H5P_DEFAULT, data);
    for(i=0 ; i<10 ; i++)
    printf("%f\n",data[i]);

    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Gclose (grp_id1);
    H5Fclose (file_id);

    return data;
}

void save_hdf5(STAR_DATA star, long int lz, double *jd)
{
  hid_t file_id, grp_id1, grp_id2, grp_id3, grp_id4;
  printf(" %s\n",star.name);
  long int i, p;
  file_id = H5Fcreate (star.name, H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  grp_id1 = H5Gcreate2 (file_id, "CURVE_CENTER", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  grp_id2 = H5Gcreate2 (file_id, "STAR_CENTER", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  grp_id3 = H5Gcreate2 (file_id, "SNR_DATA", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  grp_id4 = H5Gcreate2 (file_id, "DATA", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  double *fs= JD_hdf5(star.fluxes,jd,star.naxes[2]);
  double *ms= JDmag_hdf5(star.curve,jd,star.err,star.naxes[2]);
  double *dif= JD_hdf5(star.diff,jd,star.naxes[2]);

  save_txt(star.txt,jd,star.curve,star.err,star.naxes[2]);

   save_image_fl(grp_id2,star.star,star.naxes[2],star.naxes[0],star.naxes[1],1,"STAR");
   save_image_do(grp_id1,fs,1,star.naxes[2],2,1,"FLUXES");
   save_image_do(grp_id1,ms,1,star.naxes[2],3,1,"MAGNITUDE");
   save_image_do(grp_id1,dif,1,star.naxes[2],2,1,"DIF");
   //save_vec_f(file_id,grp_id3,"CURVE_DIFF",star.diff,1,star.naxes[2]);
   //save_vec_f(file_id,grp_id3,"SNR",star.SNR,1,lz);
   save_vec_f(file_id,grp_id4,"Rel_Mov",star.xy,1,star.naxes[2]);
   save_vec_f(file_id,grp_id4,"Abs_Mov",star.xyabs,1,star.naxes[2]);
   save_vec_f(file_id,grp_id4,"POS_X",star.px,1,star.naxes[2]);
   save_vec_f(file_id,grp_id4,"POS_Y",star.py,1,star.naxes[2]);

   save_vec_f(file_id,grp_id3,"SNR",star.errfluxes,1,star.naxes[2]);

   /*for (i=0; i<(lz-1) ; i++)
   {
   p = star.naxes[2]/star.z[i];
   char es[]= "c_" ;
   name_section(es,(i+1));
   save_vec_f(file_id,grp_id3,es,star.C_SNR[i],1,p);
   }*/

  free(fs);
  free(ms);
  H5Gclose (grp_id1);
  H5Gclose (grp_id2);
  H5Gclose (grp_id3);
  H5Gclose (grp_id4);
  H5Fclose (file_id);
}

void star_name(char name[], int x , int y, int B, int mef, char star[])
{
  sprintf(star,"%s", name);
  char xy[4];
  char ext[]  ="].hdf5";
  sprintf(xy,"%d", x);
  strcat(star, xy);
  strcat(star, "_");
  sprintf(xy,"%d", y);
  strcat(star, xy);
  strcat(star, "_");
  sprintf(xy,"%d", B);
  strcat(star,xy);
  strcat(star, "_[");
  sprintf(xy,"%d", mef);
  strcat(star,xy);
  strcat(star,ext);
  return;
}

void txt_name(char name[], int x , int y, int B, int mef, char star[])
{
  sprintf(star,"%s", name);
  char xy[4];
  char ext[]  ="].dat";
  sprintf(xy,"%d", x);
  strcat(star, xy);
  strcat(star, "_");
  sprintf(xy,"%d", y);
  strcat(star, xy);
  strcat(star, "_");
  sprintf(xy,"%d", B);
  strcat(star,xy);
  strcat(star, "_[");
  sprintf(xy,"%d", mef);
  strcat(star,xy);
  strcat(star,ext);
  return;
}

void name_section(char name[], long int i)
{
  char name2[] = " ";
  sprintf(name2,"%ld", i);
  strcat(name, name2);
  return;
}

double * JD_hdf5(float * fluxes, double * JD, long int dim)
{
double *data = (double *) calloc(dim*2,sizeof(double));
long int i;

for (i=0 ; i<dim ; i++)
{
   if(JD[0] == JD[1])
        data[i*2] = i;
   else
        data[i*2] = JD[i];

   data[(i*2)+1] = fluxes[i];
}

return data;
}

double * JDmag_hdf5(float * fluxes, double * JD, float *error, long int dim)
{
double *data = (double *) calloc(dim*3,sizeof(double));
long int i;

for (i=0 ; i<dim ; i++)
{
   if(JD[0] == JD[1])
        data[i*3] = i;
   else
        data[i*3] = JD[i];

   data[(i*3)+1] = fluxes[i];
   data[(i*3)+2] = error[i];
}

return data;
}

void save_txt(char name[], double *C1, float *C2, float *C3, long int dim)
{
    long int i;
    FILE *file = fopen(name, "wt");
    fprintf(file, "#TJD            Mag        Err\n");
    for(i=0 ; i< dim; i++)
    {
    fprintf(file, "%lf  %f  %f\n", C1[i], C2[i], C3[i]);
    }
    fclose(file);

return;
}

void save_ALCDEF(STAR_DATA star, long int lz, double *jd, HEADER header)
{
//save_txt(star.txt,jd,star.curve,star.err,star.naxes[2]);
    long int i;
    FILE *file = fopen(star.txt, "wt");

    fprintf(file, "STARTMETADATA\n");
    fprintf(file, "SUBMITPDS=TRUE\n");
    fprintf(file, "OBJECTNUMBER=%s\n",header.object);
    fprintf(file, "OBJECTNAME=%s\n",header.object);
    fprintf(file, "SESSIONDATE=%s\n",header.DateObs);
    fprintf(file, "SESSIONTIME=%s\n",header.UT);
    fprintf(file, "CONTACTNAME=%s\n",header.Observer);
    fprintf(file, "CONTACTINFO=NONE\n");
    fprintf(file, "OBSERVERS=%s\n",header.Observer);
    fprintf(file, "FACILITY=%s\n",header.Origin);
    fprintf(file, "MPCCODE=NONE\n");
    fprintf(file, "OBSLATITUDE=%s\n",header.Latitud);
    fprintf(file, "OBSLONGITUDE=%s\n",header.Longitud);
    fprintf(file, "TELESCOPE=%s\n",header.Telescop);
    fprintf(file, "DETECTOR=%s\n",header.CCDType);
    fprintf(file, "EXPOSURE=%s\n",header.ExpTime);
    fprintf(file, "EXPJD=NONE\n");
    fprintf(file, "OBJECTRA=%s\n",header.RA);
    fprintf(file, "OBJECTDEC=%s\n",header.Dec);
    fprintf(file, "FILTER=%s\n",header.Filter);
    fprintf(file, "MAGBAND=%s\n",header.Filter);
    fprintf(file, "DIFFERMAGS=FALSE\n");
    fprintf(file, "STANDARD=INTERNAL\n");
    fprintf(file, "LTCAPP=NONE\n");
    fprintf(file, "LTCTYPE=NONE\n");
    fprintf(file, "REDUCEDMAGS=NONE\n");
    fprintf(file, "DELIMITER=PIPE\n");
    fprintf(file, "ENDMETADATA\n");
    for(i=0 ; i< star.naxes[2]; i++)
    {
    fprintf(file, "DATA=%lf|%f|%f\n", jd[i], star.curve[i], star.err[i]);
    }
    fprintf(file, "ENDDATA\n");
    fclose(file);

return;
}

void save_ALCDEF_complete(STAR_DATA *star, long int cant, double *jd, HEADER header, int nobj, int aper,char name[])
{
//save_txt(star.txt,jd,star.curve,star.err,star.naxes[2]);
    long int i,j;
    char fileout[500];
    aper_name(name,star[aper*nobj].naxes[0],fileout);

    FILE *file = fopen(fileout, "wt");
    star[aper*nobj].naxes[0];


    fprintf(file, "STARTMETADATA\n");
    fprintf(file, "SUBMITPDS=TRUE\n");
    fprintf(file, "OBJECTNUMBER=%s\n",header.object);
    fprintf(file, "OBJECTNAME=%s\n",header.object);
    fprintf(file, "SESSIONDATE=%s\n",header.DateObs);
    fprintf(file, "SESSIONTIME=%s\n",header.UT);
    fprintf(file, "CONTACTNAME=%s\n",header.Observer);
    fprintf(file, "CONTACTINFO=NONE\n");
    fprintf(file, "OBSERVERS=%s\n",header.Observer);
    fprintf(file, "FACILITY=%s\n",header.Origin);
    fprintf(file, "MPCCODE=NONE\n");
    fprintf(file, "OBSLATITUDE=%s\n",header.Latitud);
    fprintf(file, "OBSLONGITUDE=%s\n",header.Longitud);
    fprintf(file, "TELESCOPE=%s\n",header.Telescop);
    fprintf(file, "DETECTOR=%s\n",header.CCDType);
    fprintf(file, "EXPOSURE=%s\n",header.ExpTime);
    fprintf(file, "EXPJD=NONE\n");
    fprintf(file, "OBJECTRA=%s\n",header.RA);
    fprintf(file, "OBJECTDEC=%s\n",header.Dec);
    fprintf(file, "FILTER=%s\n",header.Filter);
    fprintf(file, "MAGBAND=%s\n",header.Filter);
    fprintf(file, "DIFFERMAGS=FALSE\n");
    fprintf(file, "STANDARD=INTERNAL\n");
    fprintf(file, "LTCAPP=NONE\n");
    fprintf(file, "LTCTYPE=NONE\n");
    fprintf(file, "REDUCEDMAGS=NONE\n");
    fprintf(file, "DELIMITER=PIPE\n");
    fprintf(file, "ENDMETADATA\n");

    for(j=(aper*nobj) ; j<((aper*nobj)+nobj) ; j++)
    {
        for(i=0 ; i< star[j].naxes[2]; i++)
        {
            fprintf(file, "DATA=%lf|%f|%f\n", jd[i], star[j].curve[i], star[j].err[i]);
            //printf("DATA=%lf|%f|%f\n", jd[i], star[j].curve[i], star[j].err[i]);
        }
    }

    fprintf(file, "ENDDATA\n");
    fclose(file);

return;
}

void aper_name(char name[], long int x, char star[])
{
  sprintf(star,"%s", name);
  char xy[4];
  char ext[]  =".dat";
  sprintf(xy,"%ld", x);
  strcat(star, xy);
  strcat(star,ext);
  return;
}
