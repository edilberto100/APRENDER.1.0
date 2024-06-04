/**************************************************************
 *                           APPHI                            *
 * AN AUTOMATED PHOTOMETRY PIPELINE FOR HIGH CADENCE,         *
 *                     LARGE VOLUME DATA                      *
 * Edilberto Sanchez Moreno                                   *
 **************************************************************/

#ifndef APRENDER_H_
#define APRENDER_H_

#include <string.h>
#include <stdio.h>
#include <fitsio.h>
#include <hdf5.h>
#include <math.h>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>

struct arguments
{
  char *bias;
  char *dark;
  char *flat;
  char *science;
  char *ist;
  char *d1;
  int d2;
  int d3;
  double um;
  int np;
  int tinteger;
};

typedef struct
{
char object[256];
char DateObs[256];
char UT[256];
char Observer[256];
char Origin[256];
char Latitud[256];
char Longitud[256];
char Telescop[256];
char CCDType[256];
char ExpTime[256];
char RA[256];
char Dec[256];
char Equinox[256];
char Filter[256];
}HEADER;

typedef struct
{
float x;
float y;
}PROFIL;

typedef struct
{
long int x;
long int y;
}PROFILE;

typedef struct
{
float *x;
float *y;
}CORD;

typedef struct
{
char name[500];
}LIST_TXT;

typedef struct
{
int dim;
long int naxes[3];
long int *data;
}FILE_FITS;

typedef struct
{
long int naxes[3];
float * data;
}MASTER_DATA;

typedef struct
{
char name[500];
char txt[500];
long int naxes[3];
int x;
int y;
float *star;
float *curve;
float *fluxes;
float *errfluxes;
float *err;
float *diff;
float *SNR;
float **C_SNR;
long int *z;
float *xy;
float *px;
float *py;
float *xyabs;
}STAR_DATA;

typedef struct
{
long int *x;
long int *y;
long int *xe; //Xfin de trayectoria
long int *ye; //Yfin de trayectoria
long int *x0; //Xfin de trayectoria
long int *y0; //Yfin de trayectoria
double *dt; //delta pixel
long int *bin;
long int cant;
long int *D2;
long int *D3;
}STAR_BIN;

typedef struct
{
int Rank;
hsize_t *naxes;
float *data;
}HDF5_DATA_FL;

//This functions are in apphi_field.c
int apphi_field(struct arguments arguments);

//This functions are in Aphhi_fits.c
 long int * get_diameter(char * data, int * dim);
 LIST_TXT ** get_file_mef(char *filetxt, int mef);
 void name_mef(char name[], int i);
 int get_mef_size(char *filetxt);
 STAR_BIN read_lst(char *filetxt, int *nlines);
 FILE_FITS *read_files(char *filetxt, int *nlines);
 LIST_TXT *read_list(char *filetxt, int *nlines);
 int get_nlines(char *filetxt);
 FILE_FITS * get_data_fits(LIST_TXT *filesname, int * nfiles);
 void free_FILE_FITS(FILE_FITS * data, int N);
 void file_name(char name[], char star[]);
 void save_fits_ld(FILE_FITS imagen, char name[]);
 double get_JD(char FileFits[], long int m);
 double get_EXPTIME(char FileFits[], long int m);
  HEADER get_Header(char FileFits[], long int m);
 void supr_text(char *id, char *parameter);
 void get_value(char id[], char parameter[]);

 //This functions are in Aphhi_reduction.c
 MASTER_DATA master_bias(FILE_FITS * bias, int nfiles);
 MASTER_DATA master_darks(MASTER_DATA  bias, FILE_FITS * darks, int nfiles);
 MASTER_DATA master_flats(MASTER_DATA  bias, MASTER_DATA darks, FILE_FITS * flats, int nfiles);
 MASTER_DATA * image_reduction(MASTER_DATA bias, MASTER_DATA darks, MASTER_DATA flats, FILE_FITS * night, int nfiles);
 void free_MASTER_DATA(MASTER_DATA * data, int N);

 //This functions are in Aphhi_hdf5.c
 void save_image_fl(hid_t file_id, float *pix, int NOBS, int WSIZE_x, int WSIZE_y, int NTEL, char name[]);
 void save_image_do(hid_t file_id, double *pix, int NOBS, int WSIZE_x, int WSIZE_y, int NTEL, char name[]);
 void save_image_i(hid_t file_id, int *pix, int NOBS, int WSIZE_x, int WSIZE_y, int NTEL, char name[]);
 void save_vec_li(hid_t file_id, hid_t grp_id, char name[], long int *lightcurve, int NOBS, int NTEL);
 void save_vec_f(hid_t file_id, hid_t grp_id, char name[], float *lightcurve, int NOBS, int NTEL);
 void save_vec_lf(hid_t file_id, hid_t grp_id, char name[], double *lightcurve, int NOBS, int NTEL);
 void save_vec_i(hid_t file_id, hid_t grp_id, char name[], int *lightcurve, int NOBS, int NTEL);
 void get_vec_f(char * file);
 float * get_img_f(char * file , char * section);
 void save_hdf5(STAR_DATA star, long int lz, double *jd);
 void star_name(char name[], int x , int y, int B, int mef, char star[]);
 void txt_name(char name[], int x , int y, int B, int mef, char star[]);
 void name_section(char name[], long int i);
 double * JD_hdf5(float * fluxes, double * JD, long int dim);
 double * JDmag_hdf5(float * fluxes, double * JD, float *error, long int dim);
 void save_txt(char name[], double *C1, float *C2, float *C3, long int dim);
 void save_ALCDEF(STAR_DATA star, long int lz, double *jd, HEADER header);


 void scan_group(hid_t gid);
 void scan_attrs(hid_t oid);
 void do_link(hid_t gid, char *name) ;
 void do_dset(hid_t did);
 void do_dtype(hid_t tid);
 void do_attr(hid_t aid) ;
 void do_plist(hid_t pid);

//This functions are in Aphhi_search.c
STAR_BIN read_position(struct arguments arguments);
STAR_BIN conv_img(long int *data, float umbrall, long int *naxes, int D , float *comp);
float * create_gaussean(int d1);
float * create_mask_circle(int d1);
int * create_mask_ring(int d1, int d2);
void free_STAR_BIN(STAR_BIN datos);
STAR_BIN select_stars(STAR_BIN comp, long int * naxes , float *night);
float get_max(long int *bin,long int *px,long int *py, float *night, long int * naxes);


// Agregadas por Elias
CORD get_perfil(float* frame,long int* naxes);
float get_maximum(float* P,long int N);
float get_minimum(float* P,long int N);
float get_avg(float *V,long int N);
void normal_vec(float* P, long int N);
void filter_hp(float* P,long int N, float Umb);
long int* filter_positions(float* P,long int N, float Umb, long int *Cant, long int A);
STAR_BIN cut_star(STAR_BIN pst, PROFILE tam, float* FrameTotal, long int* naxes);
STAR_BIN obtain_stars(STAR_BIN stars, float* FrameTotal, long int* naxes, long int box);
void free_CORD(CORD v);
float get_moda(float *V , long int N);

//This functions are in aprender_photemtry.c
 //void get_centroids_list(int x, int y, MASTER_DATA * night, int N, int D, long int *naxes, float * xy, float * PX, float * PY, float * xyabs);
void get_centroids_list(int x, int y, MASTER_DATA * night, int N, int D, long int *naxes, float * xy, float * PX, float * PY, float * xyabs, double exptime);
void get_naxes(MASTER_DATA * night, int N, long int *naxes);
float * get_st(MASTER_DATA * night, int N, int D, long int *naxes, float * PX, float * PY);
float * get_star_xy(int x, int y, MASTER_DATA * night, int N, int D, long int *naxes, float * xy, float * PX, float * PY, float * xyabs);
float * get_star(int x, int y, MASTER_DATA * night, int N, int D, long int *naxes );
float * get_curve(float * star, float * sky , long int *naxes_star, long int *naxes_sky, long int D3);
float * get_mag(float * star, float * sky , long int *naxes_star, long int *naxes_sky, long int D3,double ti, float * err, float * errfluxes);
//float * get_mag(float * star, float * sky , long int *naxes_star, long int *naxes_sky, long int D3, int ti);
float * get_fluxes(float * star, float * sky , long int *naxes_star, long int *naxes_sky, long int D3);
void get_centroids_bin(float *bin, int B, long int *dx, long int *dy);
void supr_outfl(float *curva, long int naxes);
void free_STAR_DATA(STAR_DATA star,long int lz);
float * get_diff_curve(float * star, float * comp, long int *naxes);
long int * size_snr(long int z, long int lz);
float ** sign_noise(float *star, float *sky, long int *naxes, long int D, long int *z, long int lz, long int D3);
float * copy_vec(float * vec, long int *naxes);
float sign_curve(float *datos, long int z);
void save_fits_fl(FILE_FITS imagen, char name[]);

#endif
