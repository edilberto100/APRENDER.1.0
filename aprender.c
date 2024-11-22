#include <stdlib.h>
#include <argp.h>
#include "aprender.h"

const char *argp_program_version =
  "APPHi_field 1.0";
const char *argp_program_bug_address =
  "<taos2_apphi@astrosen.unam.mx>";

/* The options we understand. */
static struct argp_option options[] = {
  {"biaslist",  'b', "STRING",      0,
   "File list that include the Bias Files" },
  {"darklist",  'd', "STRING",      0,
   "File list that include the Dark Files" },
  {"flatlist",  'f', "STRING",      0,
   "File list that include the Flat Files" },
  {"sciencelist",   's', "STRING", 0,
  "File list that include the Scientfic Data" },
  {"starlist",  'l', "STRING",      0,
   "File list that include the star position x,y,dim" },
  {"list D1",  'x', "STRING",      0,
   "Diameter of a star [d1,d2,d3, ... , dn]." },
  {"D2",  'y', "int",      0,
   "Diametro menor del anillo" },
  {"D3",  'z', "int",      0,
   "Diametro mayor del anillo" },
  {"Integration Time",  't', "int",      0,
   "Integration time in seconds" },
  //{"umbral",   'u', "FLOAT", 0,
  //   "Binarization Umbral to detect stars (0.0 to 1.0)" },
  {"npixels",   'n', "INTEGER", 0,
     "Number of stack to calculate the SNR curve " },
  { 0 }
};


static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'b':
      arguments->bias = arg;
      break;
    case 'd':
      arguments->dark = arg;
      break;
    case 'f':
      arguments->flat = arg;
      break;
    case 's':
      arguments->science = arg;
      break;
    //case 'u':
    //  arguments->um = atof(arg);
      break;
    case 'n':
      arguments->np = atoi(arg);
      break;
    case 'l':
      arguments->ist =arg;
      break;
    case 'x':
      arguments->d1 =arg;
      break;
    case 'y':
      arguments->d2 =atoi(arg);
      break;
    case 'z':
      arguments->d3 =atoi(arg);
      break;
    case 't':
      arguments->tinteger =atoi(arg);
      break;
    case ARGP_KEY_ARG:
      /* if (state->arg_num >= 1)
       *{
       *  argp_usage(state);
       }*/
      break;
    case ARGP_KEY_END:
      /*if (state->arg_num < 2)
       *{
       *  argp_usage (state);
       }*/
      break;
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static char doc[] =
  "APPHi_field -- Automated Photometry Pipeline for High Cadence and  Large Volume Data -\va pipeline for  TAOSII Proyect http://taos2.astrosen.unam.mx";

static char args_doc[] = " ";

static struct argp argp = { options, parse_opt, args_doc, doc };

int main (int argc, char **argv)
{
  struct arguments arguments;

  arguments.bias="";
  arguments.flat="";
  arguments.dark="";
  arguments.science="";
  arguments.ist="";
  arguments.d1="";

  //argp_usage (state);

  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  printf("\n You arguments are the following:");
  printf("\n Bias list File: %s", arguments.bias);
  printf("\n Flat list File: %s", arguments.flat);
  printf("\n Dark list File: %s", arguments.dark);
  printf("\n Science Data list File: %s", arguments.science);
  printf("\n Stars list File: %s", arguments.ist);
  printf("\n Integration time: %d seg", arguments.tinteger);
  printf("\n D1 list File: %s", arguments.d1);
  printf("\n D2 : %d", arguments.d2);
  printf("\n D3: %d", arguments.d3);
  //printf("\n Umbral: %f", arguments.um);
  //printf("\n Number of pixel to generate LC: %d\n", arguments.np);

  /*STAR_BIN comp, data;
  int nist, ndiam,i,j;
  comp = read_lst(arguments.ist,&nist);
  long int * diameter= get_diameter(arguments.d1,&ndiam);
  for(i=0; i<ndiam;i++)
    printf("%ld \n",diameter[i]);*/

  /*data.cant=ndiam*comp.cant;
  data.bin=(long int*) calloc(data.cant,sizeof(long int));
  data.x=(long int*) calloc(data.cant,sizeof(long int));
  data.y=(long int*) calloc(data.cant,sizeof(long int));
  data.D2=(long int*) calloc(data.cant,sizeof(long int));
  data.D3=(long int*) calloc(data.cant,sizeof(long int));

  for(i=0; i<comp.cant;i++)
  {
      for(j=0; j<ndiam;j++)
        {
            data.x[(j*comp.cant)+i]=comp.x[i];
            data.y[(j*comp.cant)+i]=comp.y[i];
            data.bin[(j*comp.cant)+i]=diameter[j];
            data.D2[(j*comp.cant)+i]=arguments.d3;
            data.D3[(j*comp.cant)+i]=arguments.d2;
        }
  //  printf("%ld \n",diameter[i]);
  }*/
  /*for(i=0; i<data.cant;i++)
  {
    printf("Dato[%d] x=%ld, y=%ld, D1=%ld, D2=%ld, D3=%ld \n",i,data.x[i],data.y[i],data.bin[i],data.D2[i],data.D3[i]);
  }*/

  //STAR_BIN data = read_position(arguments);

  apphi_field(arguments);

  //free(data.x);
  //free(data.y);
  //free(data.bin);

  //free(data.D2);
  //free(data.D3);

  exit(0);

}
