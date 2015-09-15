#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "structs.h"

int close_streams(t_files filepkg)
{
  int i;

  for(i=0;i<filepkg->nrfiles;i++) {
    if(filepkg->bOpen[i]){
      if(filepkg->ftype[i]=='w')
        fprintf(filepkg->fil[i],"\n");
      fclose(filepkg->fil[i]);
    }
  }

  return 0;
}

int dealloc_data(t_frame pframes, t_clusters pclusters, int bQM)
{
  int i;

  free(pframes->r);
  free(pframes->v);
  free(pframes->f);
  free(pframes->resnr);
  free(pframes->restyp);
  free(pframes->partyp);
  free(pframes->Q);
  free(pframes->M);
  free(pframes->C6);
  free(pframes->C12);
  free(pframes->partyp_key);

  if(bQM){
    for(i=0;i<pframes->nr_parts;i++){
      free(pclusters->qmcluster[i]->qmndx);
      free(pclusters->qmcluster[i]->clndx);
      free(pclusters->qmcluster[i]);
    }
    free( pclusters->qmcluster );
    free( pclusters->dist_mat );
    free( pclusters->bool_mat );
    free( pclusters->indx_mat );
    free( pclusters->rQM );
  }

  return 0;
}
