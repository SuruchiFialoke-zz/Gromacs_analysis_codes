/* This Code analyzes trajectories from MD Simulations Using Gromacs 
 * Trajectory data is saved in XTC files 
 * Outpus radial distribution functions
 * by Suruchi Fialoke 
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <ctime>
#include <malloc.h>
#include </opt/gromacs/4.5.5/include/xdrfile/xdrfile_xtc.h> // xdr include file locations

        using namespace std;

        int main(int argc, char **argv) {

        /* XTC variables */
        XDRFILE *xd;                  // the xdr file holder
        int natoms=0;                 // number of atoms
        int step;                     // the step counter for the simulation
        float time;                   // simulation time
        matrix box;                   // box coordinates in a 3x3 matrix
        rvec *xtc_coor;               // atom coordinates in a 2-D matrix
        float prec;                   // precision of the xtc file

        cout<<"******* 3D g(r) of water(NVT, 20ns)*********\n";

        /* Program variables */

        int num_water= 4000 ;           //number of water in the simulation box
        float L = 5.0 ;                 //box length

        int firstsnap,lastsnap ;        //first and last snapshots to loop
        cout<<"Enter first time snap\n";
        cin>>firstsnap ;
        cout<<"Enter last time snap\n";
        cin>>lastsnap ;
        int snap_num = 0 ;              //accounts for the frame that is being looped
        int snap_count = 0 ;            // accounts from the frames that are counted

        int oxy_in_frame = 0 ;          //keeps a track of # of oxygen counted
        int n_bins = 1000 ;             // number of bins in the histogram
        float dr = 0.5*L/n_bins ;       //delta r for increments in r
        double g_r[n_bins] ;            // g(r) stores # in each bin 
        float rho_bulk = num_water/(L*L*L) ;    // bulk density

        double rij,rijsq,dx,dy,dz ;      //variables for instantaneous rij calculations

        /* End of Program variables */

        for(int k=0;k<n_bins;k++)
                {
                        g_r[k] = 0;
                }

        /* Memory allocation to read coordinates */
        read_xtc_natoms(argv[1], &natoms);
        xtc_coor = (rvec *) malloc(natoms*sizeof(rvec));

        if(xtc_coor==0)
                {
                cout<<"Insufficient memory to load .xtc file.\n";
                return 0;
                 }
        /* Open the xtc file and loop through each frame. */
        xd=xdrfile_open(argv[1],"r");
        while( ! read_xtc(xd, natoms, &step, &time, box, xtc_coor, &prec) )
        {
                snap_num++ ;
                if((snap_num>firstsnap)&&(snap_num<lastsnap))
                {
                        snap_count++ ;
                        oxy_in_frame = 0 ;

                        for (int i=0; i<num_water -1; i++)
                        {
                                for (int j=i+1; j<num_water; j++)
                                 {      /*Correcting for Periodic Boindary Conditions to calculate distances*/
                                        dx = (xtc_coor[3*i][0]-xtc_coor[3*j][0]) ;
                                        dx = dx - L*round(dx/L);                             
                                        dy = (xtc_coor[3*i][1]-xtc_coor[3*j][1]) ;
                                        dy = dy - L*round(dy/L);
                                        dz = (xtc_coor[3*i][2]-xtc_coor[3*j][2]) ;
                                        dz = dz - L*round(dz/L);

                                        rijsq = dx*dx+dy*dy+dz*dz ;

                                        rij=sqrt(rijsq) ;
                                        if(rij<0.5*L){
                                                g_r[(int)(rij/dr)] += 2.0;
                                                oxy_in_frame++ ;
                                        }
                                   }
                             }
                }

         }
         if(step==0) {
                cout<<step<<"\t"<<natoms<<endl;
                cout<<box[0][0]<<"\t"<<box[1][1]<<"\t"<<box[2][2]<<endl;
         }

        for(int k=0;k<n_bins;k++)
        {
                g_r[k] /= snap_count*num_water;
                g_r[k] /= 4.0*M_PI*k*k*pow(dr,3)*rho_bulk;
                cout<<dr*(k+0.5)<<"\t"<<g_r[k]<<endl;
        }
        return 0 ;
}
                                      
