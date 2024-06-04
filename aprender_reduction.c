/**************************************************************
 *                           APPHI                            *
 * AN AUTOMATED PHOTOMETRY PIPELINE FOR HIGH CADENCE,         *
 *                     LARGE VOLUME DATA                      *
 * Edilberto Sanchez Moreno                                   *
 **************************************************************/
 #include "aprender.h"

MASTER_DATA master_bias(FILE_FITS * bias, int nfiles)
{
MASTER_DATA master;

master.naxes[0]=bias[0].naxes[0];
master.naxes[1]=bias[0].naxes[1];
master.data = (float *) calloc(master.naxes[0]*master.naxes[1],sizeof(float));

long int i,j,k;
long int total = 0;

for (k=0 ; k< nfiles ; k++)
{
total += bias[k].naxes[2];
}
float *vector =  (float *) calloc(total,sizeof(float));
long int med = (total/2)-1;

for (i=0 ; i < master.naxes[0]*master.naxes[1]  ; i++)
{
long int plus = 0;
float temp;

    for (k=0 ; k< nfiles ; k++)
    {
        for (j=0 ; j< bias[k].naxes[2] ; j++)
        {
            vector [plus] = bias[k].data[(j*master.naxes[0]*master.naxes[1])+i];
            //if(i==1)
                //printf("%f \n",vector[plus]);
            plus++;
        }
    }

    for (k=0 ; k<total ; k++)
    {
        for (j=0 ; j<total-1 ; j++)
        {
            if (vector[j] > vector[j+1])
            {
                temp = vector[j];
                vector[j] = vector[j+1];
                vector[j+1] =temp;
            }
        }
    }

    //for (k=0 ; k<total ; k++)
    //    if(i==1){
            //printf("%f \n",vector[k]);
    //        }
    master.data[i] = vector[med];

}
free(vector);

/****************************
  hid_t file_id, grp_id;
  file_id = H5Fcreate ("Bias.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  save_image_fl(file_id,master.data,1,master.naxes[0],master.naxes[1],1,"star");
  H5Fclose (file_id);
//****************************/
return master;
}

MASTER_DATA master_darks(MASTER_DATA  bias, FILE_FITS * darks, int nfiles)
{
MASTER_DATA master;
master.naxes[0]=darks[0].naxes[0];
master.naxes[1]=darks[0].naxes[1];

master.data = (float *) calloc(master.naxes[0]*master.naxes[1],sizeof(float));
long int i,j,k;
long int total = 0;
for (k=0 ; k< nfiles ; k++)
{
total += darks[k].naxes[2];
}
float *vector =  (float *) calloc(total,sizeof(float));
long int med = (total/2)-1;

for (i=0 ; i < master.naxes[0]*master.naxes[1]  ; i++)
{
long int plus = 0;
float temp;

    for (k=0 ; k< nfiles ; k++)
    {
        for (j=0 ; j< darks[k].naxes[2] ; j++)
        {
            vector [plus] = darks[k].data[(j*master.naxes[0]*master.naxes[1])+i] - bias.data[i];
            plus++;
        }
    }

    for (k=0 ; k<total ; k++)
    {
        for (j=0 ; j<total-1 ; j++)
        {
            if (vector[j] > vector[j+1])
            {
                temp = vector[j];
                vector[j] = vector[j+1];
                vector[j+1] = temp;
            }
        }
    }
    master.data[i] = vector[med];//- bias.data[i];
}
free(vector);

/****************************
  hid_t file_id, grp_id;
  file_id = H5Fcreate ("Darks.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  save_image_fl(file_id,master.data,1,master.naxes[0],master.naxes[1],1,"star");
  H5Fclose (file_id);
//****************************/
return master;
}

MASTER_DATA master_flats(MASTER_DATA  bias, MASTER_DATA darks, FILE_FITS * flats, int nfiles)
{
MASTER_DATA master;

master.naxes[0]=flats[0].naxes[0];
master.naxes[1]=flats[0].naxes[1];

master.data = (float *) calloc(master.naxes[0]*master.naxes[1],sizeof(float));
long int i,j,k;
long int total = 0;

for (k=0 ; k< nfiles ; k++)
{
total += flats[k].naxes[2];
}

for (i=0 ; i < master.naxes[0]*master.naxes[1]  ; i++)
{
    for (k=0 ; k< nfiles ; k++)
    {
        for (j=0 ; j< flats[k].naxes[2] ; j++)
        {
            flats[0].data[i] += flats[k].data[(j*master.naxes[0]*master.naxes[1])+i] - bias.data[i] - darks.data[i];
        }
    }
    master.data[i] = flats[0].data[i]/(total+1);
}
long int p=0;
for (k=0 ; k<master.naxes[0]*master.naxes[1] ; k++)
{
    p += master.data[k];
}
p = p/(master.naxes[0]*master.naxes[1]);

for (i=0 ; i < master.naxes[0]*master.naxes[1]  ; i++)
{
    for (k=0 ; k< nfiles ; k++)
    {
        for (j=0 ; j< flats[k].naxes[2] ; j++)
        {
            flats[0].data[i] += (flats[k].data[(j*master.naxes[0]*master.naxes[1])+i] - bias.data[i] - darks.data[i])/p;
        }
    }
    master.data[i] = flats[0].data[i]/(total+1);
}

for (k=0 ; k<master.naxes[0]*master.naxes[1] ; k++)
{
    master.data[k] = master.data[k]/p;
}

/****************************
  hid_t file_id, grp_id;
  file_id = H5Fcreate ("Flats.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  save_image_fl(file_id,master.data,1,master.naxes[0],master.naxes[1],1,"star");
  H5Fclose (file_id);
//****************************/
return master;
}

MASTER_DATA * image_reduction(MASTER_DATA bias, MASTER_DATA darks, MASTER_DATA flats, FILE_FITS * night, int nfiles)
{
MASTER_DATA * data = (MASTER_DATA *) calloc(nfiles,sizeof(MASTER_DATA));
long long int i,j,k,l;
long long int total = 0;
for (k=0 ; k< nfiles ; k++)
{
    data[k].data = (float *) calloc(night[k].naxes[0]*night[k].naxes[1]*night[k].naxes[2],sizeof(float));
    data[k].naxes[0] = night[k].naxes[0];
    data[k].naxes[1] = night[k].naxes[1];
    data[k].naxes[2] = night[k].naxes[2];
    //printf("data_%lld: x=%ld, y=%ld, z=%ld\n",k,data[k].naxes[0],data[k].naxes[1],data[k].naxes[2]);
    if (data[k].naxes[0] > 0 && data[k].naxes[1] > 0 && data[k].naxes[2] > 0 )
    {
    for (j=0 ; j< data[k].naxes[2] ; j++)
    {
        for (i=0 ; i < data[k].naxes[1] ; i++)
        {
            for (l=0 ; l < data[k].naxes[0] ; l++)
            {
            data[k].data[(j*data[k].naxes[0]*data[k].naxes[1]) + (i*data[k].naxes[0]) + l] = (night[k].data[(j*data[k].naxes[0]*data[k].naxes[1]) + (i*data[k].naxes[0]) + l] - darks.data[(i*data[k].naxes[0]) + l] - bias.data[(i*data[k].naxes[0]) + l]) / flats.data[(i*data[k].naxes[0]) + l];
            }
        }
    total ++;
    }
    free(night[k].data);
    }
}
free(night);
//***************************
  hid_t file_id, grp_id;
  file_id = H5Fcreate ("Reduction.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  save_image_fl(file_id,data[0].data,data[0].naxes[2],data[0].naxes[0],data[0].naxes[1],1,"image_01");
  //save_image_fl(file_id,data[1].data,data[1].naxes[2],data[1].naxes[0],data[1].naxes[1],1,"image_02");
  //save_image_fl(file_id,data[2].data,data[2].naxes[2],data[2].naxes[0],data[2].naxes[1],1,"image_03");
  //save_image_fl(file_id,data[3].data,data[3].naxes[2],data[3].naxes[0],data[3].naxes[1],1,"image4");
  //save_image_fl(file_id,data[4].data,data[4].naxes[2],data[4].naxes[0],data[4].naxes[1],1,"image5");
  save_image_fl(file_id,bias.data,1,bias.naxes[0],bias.naxes[1],1,"bias");
  save_image_fl(file_id,darks.data,1,darks.naxes[0],darks.naxes[1],1,"dark");
  save_image_fl(file_id,flats.data,1,flats.naxes[0],flats.naxes[1],1,"flat");
  H5Fclose (file_id);
  /*for(i=0 ; i < data[0].naxes[1] ; i++){
    for(l=0 ; l < data[0].naxes[0] ; l++){

    printf("%f ",data[0].data[(i*data[0].naxes[0]) + l]);
    }
    printf("\n");
  }*/


//****************************/
return data;
}

void free_MASTER_DATA(MASTER_DATA * data, int N)
{
long int i;
for (i=0 ; i<N ; i++)
{
free(data[i].data);
}

return;
}

