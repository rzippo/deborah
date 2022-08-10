/* Computes average and variance for multiple replics
 * of LUDB-non-nested evaluation study
 * Output files are named "nnfull-replic_id-x-y.out"
 */


#include <stdio.h>
#include <string.h>
#include <math.h>


typedef struct info {
  int nodes;
  int perc;
  int flows;
  int cuts;
  double tcut;
  double tsimpl;
  double ludb;
  int nsimpl;
} info;

info simdata[20][512];

int main(int argc, char *argv[])
{
  char fname[260];
  int rep_start=1000, nreps=11;
  int nperc = 3;
  FILE *fin, *fout;
  int lines = 0;
  char suffix[64]="";

  if(argc<5)
  {
    printf("Usage:  %s <base filename> <replic_id_seed_start> <nr_replics> <nr_horizontal_categories (F%%)> [fname_suffix]\n",argv[0]);
    printf("   ex:  %s nnfull 1005 10 3 -0  ==>  nnfull-{1005..1014}-0.csv\n", argv[0]);
    return -1;
  }
  
  sscanf_s(argv[2],"%i",&rep_start);
  sscanf_s(argv[3],"%i",&nreps);
  sscanf_s(argv[4],"%i",&nperc);
  if(argc >= 6) sscanf_s(argv[5],"%s",suffix);
  
  printf("Processing replic IDs: %i to %i;  nr of horizontal categories = %i\n",rep_start,rep_start+nreps-1,nperc);

  for(int r=0; r<nreps; r++)
  {
    sprintf(fname, "%s-%i%s.csv", argv[1], r+rep_start, suffix);
    fin = fopen(fname, "rt");
    if(fin == nullptr) 
    {
      printf("ERROR: can't open %s\n",fname);
      break;
    }

    lines = 0;
    while(true)
    {
      int nodes = 0;
      int perc = 0;
      int flows = 0;
      int cuts = 0;
      double tcut = 0.0;
      double tsimpl = 0.0;
      double ludb = 0.0;
      int nsimpl = 0;
      
      int err = fscanf(fin, "%i,%i,%i,%i,%lf,%lf,%lf,%i", &nodes, &perc, &flows, &cuts, &tcut, &tsimpl, &ludb, &nsimpl);
      if(err <= 0) break;
      
      simdata[r][lines].nodes = nodes;
      simdata[r][lines].perc = perc;
      simdata[r][lines].flows = flows;
      simdata[r][lines].cuts = cuts;
      simdata[r][lines].tcut = tcut;
      simdata[r][lines].tsimpl = tsimpl;
      simdata[r][lines].ludb = ludb;
      simdata[r][lines].nsimpl = nsimpl;
      
      //printf("R=%i (%02i):  %i %i %i %i %lf %lf %lf %i\n", r, err, nodes, perc, flows, cuts, tcut, tsimpl, ludb, nsimpl);
      lines++; 
    }
    
    fclose(fin);
  }
  
  info sum[250];
  for(int l=0; l<lines; l++)
  {
    sum[l].nodes = simdata[0][l].nodes;
    sum[l].perc = simdata[0][l].perc;
    sum[l].flows=0;
    sum[l].cuts=0;
    sum[l].tcut=0.0;
    sum[l].tsimpl=0.0;
    sum[l].ludb=0.0;
    sum[l].nsimpl=0;
  }
  for(int r=0; r<nreps; r++)
  {
    for(int l=0; l<lines; l++)
    {
      sum[l].flows += simdata[r][l].flows;
      sum[l].cuts += simdata[r][l].cuts;
      sum[l].tcut += simdata[r][l].tcut;
      sum[l].tsimpl += simdata[r][l].tsimpl;
      sum[l].ludb += simdata[r][l].ludb;
      sum[l].nsimpl += simdata[r][l].nsimpl;
    }
  }
  for(int l=0; l<lines; l++)
  {
    sum[l].flows/=nreps;
    sum[l].cuts/=nreps;
    sum[l].tcut/=nreps;
    sum[l].tsimpl/=nreps;
    sum[l].ludb/=nreps;
    sum[l].nsimpl/=nreps;    
    printf("%i,%i,%i,%i,%lf,%lf,%lf,%i\n", sum[l].nodes, sum[l].perc, sum[l].flows, sum[l].cuts, sum[l].tcut, sum[l].tsimpl, sum[l].ludb, sum[l].nsimpl);
  }

  
  int nnodes = lines / nperc;

  printf("matrix sheet: t_cut      nodes=%i x cat=%i\n", nnodes, nperc);

  for(int n=0; n<nnodes; n++)
  {
    for(int p=0; p<nperc; p++)
    {
      int l = n*nperc+p;
      printf("%lf",sum[l].tcut);
      if(p<nperc-1) printf(",");
    }
    std::cout << std::endl;
  }

  for(int l=0; l<lines; l++)
  {
    double var1 = 0.0;
    double var2 = 0.0;
    for(int r=0; r<nreps; r++)
    {
      var1 += pow(simdata[r][l].tcut - sum[l].tcut, 2.0);
      var2 += pow(simdata[r][l].tsimpl - sum[l].tsimpl, 2.0);
    }
    var1 /= nreps;
    var2 /= nreps;
    double stddev1 = sqrt(var1);
    double stddev2 = sqrt(var2);
    printf("variance l=%02i:  t_cut=%lf (%lf)   t_simpl=%lf (%lf)\n", l, var1, stddev1/sum[l].tcut, var2, stddev2/sum[l].tsimpl);
  }

  
  printf("matrix sheet: t_simpl   nodes=%i x cat=%i\n", nnodes, nperc);
  
  for(int n=0; n<nnodes; n++)
  {
    for(int p=0; p<nperc; p++)
    {
      int l = n*nperc+p;
      printf("%lf",sum[l].tsimpl);
      if(p<nperc-1) printf(",");
    }
    std::cout << std::endl;
  }


  printf("matrix sheet: ludb_numeric   nodes=%i x cat=%i\n", nnodes, nperc);
  
  for(int n=0; n<nnodes; n++)
  {
    for(int p=0; p<nperc; p++)
    {
      int l = n*nperc+p;
      printf("%lf",sum[l].ludb);
      if(p<nperc-1) printf(",");
    }
    std::cout << std::endl;
  }

  
  return 0;
}
