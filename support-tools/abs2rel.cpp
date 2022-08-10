/* Abs2Rel - Converts per-cut LUDB delays to normalized values
 * against the minimum LUDB found.
 * Requires DEBORAH compilation without TANDEM_LUDB_TIMINGS
 */

#include <stdio.h>
#include <string.h>


typedef struct {
 int n;
 double v;
 int c;
} info;

#define NMAX 8192

int main(int argc, char*argv[])
{
char txt[512];
info data[NMAX];

if(argc < 2) {printf("Usage: %s <inputfile>\n",argv[0]); return -1;}
FILE *fin = fopen(argv[1], "rt");
if(fin == nullptr) {printf("Error opening file\n"); return -2;}

int num = 0;
double vmin = 9999999999.0;
int nmin = 0;
while(true)
{
 txt[0] = 0;
 int p1,p3;
 double p2;
 if(fgets(txt, 512, fin) == nullptr) break;
 int err = sscanf_s(txt, "%u %lf %u", &p1, &p2, &p3);
 if(err < 3) break;
 data[num].n = p1; data[num].v = p2; data[num].c = p3;
 if(data[num].v < vmin) {vmin = data[num].v; nmin=num; }
 num++;
 if(num >= NMAX) {printf("Fatal error: too many records (max %i)\n",NMAX); break;}
}
fclose(fin);

for(int i=0; i<num; i++) printf("%03i %.4lf %i\n", data[i].n,100.0*data[i].v/vmin,data[i].c);

return 0;
}
