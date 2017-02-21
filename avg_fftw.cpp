/* Suruchi 05/10/2015
This code calculates fourier transform of a matrix h(x,y)-h_avg(x,y)
h(x,y) was obtained from trajectory file and stored at each time step

This code uses C++ library fftw and saves kx^2+ky^2 in Angstrom square

Also saves <|h(kx,ky)|^2> in second column in units of Angstrom square, 

I AM USING sqrt(N) as the normalization constant
*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <ctime>
#define PI  3.1415926  /* pi */
#include <fftw3.h> 
using namespace std;
int main(int argc, char **argv) {

  int frame_count= 0 ;
  int tfirst = 501 ;        //first frame to be read
  int tlast = 7000 ;        //last frame to be read
  double Lx = 120.0 ;       //Angstrom
  double Ly = 120.0 ;
  double dL = 1.0 ;         //grid size in Angstrom used in calculating hxy
  int Nx =int(Lx/dL) ; 
  int Ny =int(Ly/dL) ; 
  string filepath = "/Users/suruchi/Desktop/Research/Capillary_waves/Analysis/L12/data/" ;
  string filename = "h_";
  string end = ".dat";

/* Read average h_xy */
  double h_avg[Nx][Ny] ;
  
  double x_temp[Nx*Ny] ;
  double y_temp[Nx*Ny] ;
  double h_temp[Nx*Ny] ;
  string current_file2=filepath+"hxy_avg"+end;
  int i=0 ;
  ifstream infile (current_file2.c_str());
  if (infile.is_open())
  {
 	while (! infile.eof())
        {	infile.ignore(10000,'\n');  //ignore first line of the file
                infile >>x_temp[i]>>y_temp[i]>>h_temp[i] ;
                i ++;
        }
       infile.close();
       infile.clear();
    }
   else cout<< "cant open file" << endl;
   
   for(int i=0;i<Nx;i++) {
	for(int j=0;j<Ny;j++) {
        h_avg[i][j] = h_temp[j+i*Ny] ;
        }
   }


/* Fourier transform variables */
    double h_Re_ij ;
    double h_Im_ij ;
    double kxy_sq[Nx][Ny] ;             //stores kx^2+ky^2
    double abs_hk_sq[Nx][Ny] ;          //stores instantaneous |h(kx,ky)|^2
    double Avg_hk_sq[Nx][Ny] ;          //stores <|h(kx,ky)|^2>
    double area_var ;                   //stores A-Ao=LxLy/(NxNy)^2 *sum_kx sum_ky <|h(kx,ky)|^2> kx^2+ky^2
      
    for(int i = 0; i < Nx; i++)
    {
        for(int j = 0; j < Ny; j++)
        {
            kxy_sq[i][j] = 2*PI*(i+1)/Lx * 2*PI*(i+1)/Lx + 2*PI*(j+1)/Ly *2*PI*(j+1)/Ly ;
            Avg_hk_sq[i][j] = 0.0 ;
        }
    }
    /* instantaneous h_xy calculations  */
    
    string current_file;
    for(int t=tfirst ; t<tlast;t++) {
        frame_count++ ;
        stringstream num;
        num<<t;
        current_file =filepath+ filename + num.str() + end;
        int p=0 ;
        ifstream infile (current_file.c_str());
        if (infile.is_open())
        {
        	while (! infile.eof())
            {
                //infile.ignore(10000,'\n');
                infile >>x_temp[p]>>y_temp[p]>>h_temp[p] ;
                p ++;
            }
       		infile.close();
       		infile.clear();
        }
        else cout<< "cant open file" << endl;
        
        
        for(int i=0;i<Nx;i++) {
            for(int j=0;j<Ny;j++) {
                h_temp[j+i*Ny] = 10.0*(h_temp[j+i*Ny] - h_avg[i][j]) ; //convert into Angstrom
            	abs_hk_sq[i][j] = 0.0 ;
                //h[i][j] = h_temp[j+i*Ny] ;
            }
        }
        
        fftw_complex *in, *out;
        fftw_plan pf;
        
        // Declare one-dimensional contiguous arrays of dimension Nx*Ny
        in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
        out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
        
        //get input data in row major 1-D form , order is important
        
        for (int i=0; i<Nx*Ny; i++) {
            in[i][0] = h_temp[i] ;          //htemp was already in 1D format
            in[i][1] = 0.0 ;                //complex part is zero
        }
        
        // create plan, an object containing necessary parameters
        pf = fftw_plan_dft_2d(Nx, Ny, in,out,FFTW_FORWARD,FFTW_ESTIMATE);
        
      
        // Perform FFT
        fftw_execute(pf);
        
        // Get output data into 2-d form

        for(int i = 0; i < Nx; i++)
        {
            for(int j = 0; j < Ny; j++)
            {
                h_Re_ij = out[j+i*Ny][0];
                h_Im_ij = out[j+i*Ny][1];
                abs_hk_sq[Nx-1-i][Ny-1-j]  =  h_Re_ij* h_Re_ij/Nx/Ny + h_Im_ij*h_Im_ij/Nx/Ny ;
            }
        }
 	
       area_var=0 ;
       for(int m=0;m<Nx;m++) {
        for(int n=0;n<Ny;n++) {
         	Avg_hk_sq[m][n]+= abs_hk_sq[m][n] ;
		area_var += (abs_hk_sq[m][n])*kxy_sq[m][n] ;
	}
	}	       
        cout<<t<<"\t"<<"\t"<<area_var*Lx*Ly<<endl ;
        fftw_destroy_plan(pf);
        fftw_free(in);
        fftw_free(out);
       
        
    } //close frames
    
    cout<<"#frame_count  ="<<frame_count<<endl ;
    
	/*Dividing abs_sq_hk by frames counted and save to file*/
    
	string currentfile ;
	currentfile = "Hk_avg_fftw.dat";
	FILE * hfile;
    hfile = fopen (currentfile.c_str(), "w");
    fprintf (hfile, "#kx_sq+ky_sq     hk\n");
	
    for(int m=0;m<Nx;m++) {
        for(int n=0;n<Ny;n++) {
            Avg_hk_sq[m][n]/= frame_count ; //divide by total frames
            fprintf(hfile, "%8.4f\t%8.4f\n",kxy_sq[m][n],Avg_hk_sq[m][n] ) ;
        }
	}

	fclose (hfile);
	return 0 ;
}
