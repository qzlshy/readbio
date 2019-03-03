# readbio
readfile.h <br/>
a c++ library to read pdb psf dcd file <br/>
fit.h <br/>
fit a protein structure to another <br/>

An example: <br/>
```
#include "readfile.h"
#include "fit.h"

int main(int argc,char **argv) {
  string res[25]= {"ALA","ARG","ASN","ASP","CYS", \
                        "GLN","GLU","GLY","HSD","HSE", \
                        "HSP","ILE","LEU","LYS","MET", \
                        "PHE","PRO","SER","THR","TRP", \
                        "TYR","VAL","HIS","HIE","MSE"};
                        
  Read_psf psf;
  psd.loadfile("./myfile.psf");
  
  int *cidx,cn;
  psf.choose("resname",res,25);  //choose protein
  psf.choose("segid","C");       //choose segment
  psf.choose("name","CA");       //choose CA atom
  cn=psf.choose_getn();
  cidx=new int[cn];
  psf.choose_get(cidx);
  psf.choose_reset();
  
  int frame;
  float *x,*y,*z;
  float *x_t,*y_t,*z_t;
  Read_dcd dcd;
  frame=dcd.loadfile("./myfile.dcd");
  x=new float[cn];
  y=new float[cn];
  z=new float[cn];
  x_t=new float[cn];
  y_t=new float[cn];
  z_t=new float[cn];
  dcd.readstep();
  dcd.getstep(x,y,z,cidx,cn);
  dcd.closefile();
  
  dcd.loadfile("./myfile.dcd");
  for(int i=0;i<frame;i++) {
    dcd.readstep();
    dcd.getstep(x_t,y_t,z_t,cidx,cn);
    fit(x,y,z,x_t,y_t,z_t,x_t,y_t,z_t,cn,cn);
  }
}

```
