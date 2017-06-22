/* This instantaneous interface code is for individual trajectories
 * generates one interface xtc file per xtc trajectory
 * Uses Marching cubes algorithm to search for interface
 * Works if pillar atoms are frozen in all direcion, calculates pillar density in first step
 */

// Suruchi 09/05/2013

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <ctime>
#include <malloc.h>
#include <xdrfile_xtc.h> // xdr include file

using namespace std;
#include "miscfun.h"
#include "Marching_Cube.h"

int main(int argc, char **argv) {
    
    
    /* XTC variables */
    XDRFILE *xdin;                  // the xdr file holder
    XDRFILE *xdout;
    int natoms=0;                 // number of atoms
    int step;                     // the step counter for the simulation
    float time;                   // simulation time
    matrix box;                   // box coordinates in a 3x3 matrix
    rvec *xtc_coor;               // atom coordinates in a 2-D matrix
    float prec;                   // precision of the xtc file
    
    /* variables from grofiles*/
    
    int frame = 0;
    int skip_frame = 0;
    int frame_count = 0;
    int npil = 768; // VARIABLE
    int nwater =1589; // VARIABLE
    
    // variables for look up table, function calphi defined in miscfun.h
    
    double l;
    double sigw = 0.24;
    double sigp = 0.16;  //slightly less than water
    double cutoff = 0.7;
    double dl = 0.01;
    int npoints = int(cutoff/dl) + 1;
    double LUT_phiw[npoints];
    double LUT_phip[npoints];
    for (int i=0; i<npoints; i++)
    {
        l = i*dl;
        LUT_phiw[i] = calphi(l, sigw, cutoff);
        LUT_phip[i] = calphi(l, sigp, cutoff);
    }
    
    // set up grids for coarsed grained density
    double dL = 0.1;
    double dgrid[3]; //contains delta x, delta y and delta z
    int ngrids[3];  //number of grids in each directions
    
    // variables for coarse grained density calculations
    double*** rhow;
    double*** rhop;
    double*** rho;
    int index[3];
    int iwater;
    int Ninc = int(cutoff/dL);
    double phix, phiy, phiz;
    double rx, ry, rz;
    int nx, ny, nz;
    int nrx, nry, nrz;
    
    // Marching cubes variables
    
    double rhow0 = 33.4; // #/nm^3
    double rhop0 = 40.0; // WARNING: need to test diff rhop0!
    double rhoc = 0.5; //global cut-off
    int NII = 30000; // VARIABLE 50000 for large system , number of marching cubes coordinates
    int II = 0;
    rvec* ii_coor;
    int i1, j1, k1;
    int cubefactor;
    int vertexflag;
    int cubecount;
    double vertexdist;
    double avg_nii = 0.0;
    double avg_area = 0.0;
    double avg_vol = 0.0;
    double area, volume;
    double a2, b2, c2;
    
    
    /* variables that go into marching cube function*/
    
    int ntri;
    CUBEINFO grid;
    TRIANGLE tri[5]; //contains co-ordinates of triangles that from the interface
    
    
    
    /* Memory allocation to read coordinates */
    read_xtc_natoms(argv[1], &natoms);
    xtc_coor = (rvec *) malloc(natoms*sizeof(rvec));
    
    if(xtc_coor==0){
        cout<<"Insufficient memory to load .xtc file.\n";
        return 0;
    }
    
    
    /* Open the xtc file and loop through each frame. */
    
    xdin=xdrfile_open(argv[1],"r"); //read_file
    xdout=xdrfile_open("PWII_NNN.xtc","w"); //write to file
    
    cout << "# Time\tNii\tavg_nii\tarea\tavg_area\tvolume\tavg_vol" << endl;
    
    while( ! read_xtc(xdin, natoms, &step, &time, box, xtc_coor, &prec) )
    {
        frame ++;
        if (frame > skip_frame)
        {
            frame_count ++;
            
            // find the box size and set up grids
            for (int i=0; i<3; i++)
            {
                ngrids[i] = int(box[i][i]/dL);
                dgrid[i] = box[i][i]/ngrids[i];
            }
            rhow = allocate3D(ngrids[0], ngrids[1], ngrids[2]);
            //rhop = allocate3D(ngrids[0], ngrids[1], ngrids[2]);
            rho = allocate3D(ngrids[0], ngrids[1], ngrids[2]);
            
            // compute coarse grain density from pillar atoms
            // pillar atoms are fixed so calculate only for the first-step
            if(frame_count==1) {
                rhop = allocate3D(ngrids[0], ngrids[1], ngrids[2]);
                for (int i=0; i<npil; i++)
                {
                    for (int dim=0; dim<3; dim++)
                    {
                        index[dim] = int( xtc_coor[i][dim]/dgrid[dim] );
                    }
                    
                    for (int xg=0; xg<2*Ninc; xg++)
                    {
                        nx = index[0] - Ninc + 1 + xg;
                        if (nx < 0) nx += ngrids[0];
                        else if (nx >= ngrids[0]) nx -= ngrids[0];
                        
                        rx = abs( (Ninc-1-xg)*dgrid[0] + (xtc_coor[i][0] - index[0]*dgrid[0]) );
                        nrx = int(rx/dl);
                        phix = LUT_phip[nrx] + (LUT_phip[nrx+1] - LUT_phip[nrx]) * (rx - nrx*dl) / dl;
                        
                        for (int yg=0; yg<2*Ninc; yg++)
                        {
                            ny = index[1] - Ninc + 1 + yg;
                            if (ny < 0) ny += ngrids[1];
                            else if (ny >= ngrids[1]) ny -= ngrids[1];
                            
                            ry = abs( (Ninc-1-yg)*dgrid[1] + (xtc_coor[i][1] - index[1]*dgrid[1]) );
                            nry = int(ry/dl);
                            phiy = LUT_phip[nry] + (LUT_phip[nry+1] - LUT_phip[nry]) * (ry - nry*dl) / dl;
                            
                            for (int zg=0; zg<2*Ninc; zg++)
                            {
                                nz = index[2] - Ninc + 1 + zg;
                                if (nz < 0)
                                nz += ngrids[2];
                                else if (nz >= ngrids[2])
                                nz -= ngrids[2];
                                
                                rz = abs( (Ninc-1-zg)*dgrid[2] + (xtc_coor[i][2] - index[2]*dgrid[2]) );
                                nrz = int(rz/dl);
                                phiz = LUT_phip[nrz] + (LUT_phip[nrz+1] - LUT_phip[nrz]) * (rz - nrz*dl) / dl;
                                
                                rhop[nx][ny][nz] += phix*phiy*phiz;
                            }
                        }
                    }
                }  //close npill
            } //close pillar calculations in first_frame
            //cout<<"testing" ;
            
            // compute coarse grain density from OW at each time-step
            // the first OW is at npil
            for (int i=0; i<nwater; i ++)
            {
                iwater = npil + i*3;
                for (int dim=0; dim<3; dim++)
                {
                    index[dim] = int( xtc_coor[iwater][dim]/dgrid[dim] );
                }
                
                for (int xg=0; xg<2*Ninc; xg++)
                {
                    nx = index[0] - Ninc + 1 + xg;
                    if (nx < 0) nx += ngrids[0];
                    else if (nx >= ngrids[0]) nx -= ngrids[0];
                    rx = abs( (Ninc-1-xg)*dgrid[0] + (xtc_coor[iwater][0] - index[0]*dgrid[0]) );
                    nrx = int(rx/dl);
                    phix = LUT_phiw[nrx] + (LUT_phiw[nrx+1] - LUT_phiw[nrx]) * (rx - nrx*dl) / dl;
                    
                    for (int yg=0; yg<2*Ninc; yg++)
                    {
                        ny = index[1] - Ninc + 1 + yg;
                        if (ny < 0) ny += ngrids[1];
                        else if (ny >= ngrids[1]) ny -= ngrids[1];
                        ry = abs( (Ninc-1-yg)*dgrid[1] + (xtc_coor[iwater][1] - index[1]*dgrid[1]) );
                        nry = int(ry/dl);
                        phiy = LUT_phiw[nry] + (LUT_phiw[nry+1] - LUT_phiw[nry]) * (ry - nry*dl) / dl;
                        
                        for (int zg=0; zg<2*Ninc; zg++)
                        {
                            nz = index[2] - Ninc + 1 + zg;
                            if (nz < 0) nz += ngrids[2];
                            else if (nz >= ngrids[2]) nz -= ngrids[2];
                            rz = abs( (Ninc-1-zg)*dgrid[2] + (xtc_coor[iwater][2] - index[2]*dgrid[2]) );
                            nrz = int(rz/dl);
                            phiz = LUT_phiw[nrz] + (LUT_phiw[nrz+1] - LUT_phiw[nrz]) * (rz - nrz*dl) / dl;
                            
                            rhow[nx][ny][nz] += phix*phiy*phiz;
                        }
                    }
                }
            } //close nwater
            
            // normalize and combine the density distribution
            for (int i=0; i<ngrids[0]; i++)
            {
                for (int j=0; j<ngrids[1]; j++)
                {
                    for (int k=0; k<ngrids[2]; k++)
                    {
                        rho[i][j][k] = rhop[i][j][k]/rhop0 + rhow[i][j][k]/rhow0;
                    }
                }
            }
            
            // search for ii using marching cube
            ii_coor = (rvec *) calloc(NII, sizeof(rvec));
            II = 0;
            area = 0.0;
            cubecount = 0;
            
            for (int i=0; i<ngrids[0]; i++)
            {
                for (int j=0; j<ngrids[1]; j++)
                {
                    for (int k=0; k<ngrids[2]; k++)
                    {
                        i1 = i + 1;
                        j1 = j + 1;
                        k1 = k + 1;
                        if (i1 >= ngrids[0]) i1 -= ngrids[0];
                        if (j1 >= ngrids[1]) j1 -= ngrids[1];
                        if (k1 >= ngrids[2]) k1 -= ngrids[2];
                        
                        grid.v[0] = rho[i][j][k];
                        grid.v[1] = rho[i][j1][k];
                        grid.v[2] = rho[i1][j1][k];
                        grid.v[3] = rho[i1][j][k];
                        grid.v[4] = rho[i][j][k1];
                        grid.v[5] = rho[i][j1][k1];
                        grid.v[6] = rho[i1][j1][k1];
                        grid.v[7] = rho[i1][j][k1];
                        
                        cubefactor = 0;
                        for (int v=0; v<8; v++)
                        {
                            if (grid.v[v] <= rhoc)
                            {
                                cubefactor ++;
                            }
                        }
                        if (cubefactor >=4)
                        cubecount ++;
                        
                        grid.p[0].x = i*dgrid[0]; grid.p[0].y = j*dgrid[1];grid.p[0].z = k*dgrid[2];
                        grid.p[1].x = i*dgrid[0]; grid.p[1].y = (j+1)*dgrid[1];grid.p[1].z = k*dgrid[2];
                        grid.p[2].x = (i+1)*dgrid[0];grid.p[2].y = (j+1)*dgrid[1]; grid.p[2].z = k*dgrid[2];
                        grid.p[3].x = (i+1)*dgrid[0]; grid.p[3].y = j*dgrid[1];  grid.p[3].z = k*dgrid[2];
                        grid.p[4].x = i*dgrid[0];grid.p[4].y = j*dgrid[1];grid.p[4].z = (k+1)*dgrid[2];
                        grid.p[5].x = i*dgrid[0];grid.p[5].y = (j+1)*dgrid[1];grid.p[5].z = (k+1)*dgrid[2];
                        grid.p[6].x = (i+1)*dgrid[0]; grid.p[6].y = (j+1)*dgrid[1]; grid.p[6].z = (k+1)*dgrid[2];
                        grid.p[7].x = (i+1)*dgrid[0];grid.p[7].y = j*dgrid[1];grid.p[7].z = (k+1)*dgrid[2];
                        
                        ntri = marching_cube(grid, rhoc, tri);
                        for (int t=0; t<ntri; t++)
                        {
                            a2 = (tri[t].p[0].x-tri[t].p[1].x)*(tri[t].p[0].x-tri[t].p[1].x) + (tri[t].p[0].y-tri[t].p[1].y)*(tri[t].p[0].y-tri[t].p[1].y) + (tri[t].p[0].z-tri[t].p[1].z)*(tri[t].p[0].z-tri[t].p[1].z);
                            b2 = (tri[t].p[1].x-tri[t].p[2].x)*(tri[t].p[1].x-tri[t].p[2].x) + (tri[t].p[1].y-tri[t].p[2].y)*(tri[t].p[1].y-tri[t].p[2].y) + (tri[t].p[1].z-tri[t].p[2].z)*(tri[t].p[1].z-tri[t].p[2].z);
                            c2 = (tri[t].p[2].x-tri[t].p[0].x)*(tri[t].p[2].x-tri[t].p[0].x) + (tri[t].p[2].y-tri[t].p[0].y)*(tri[t].p[2].y-tri[t].p[0].y) + (tri[t].p[2].z-tri[t].p[0].z)*(tri[t].p[2].z-tri[t].p[0].z);
                            area += 1.0/4.0 * sqrt( 2.0*(a2*b2 + b2*c2 + c2*a2) - a2*a2 - b2*b2 - c2*c2 );
                        }
                        
                        for (int l=0; l<ntri; l++)
                        {
                            for (int m=0; m<3; m++)
                            {
                                vertexflag = 0;
                                for (int n=0; n<II; n++)
                                {
                                    vertexdist = 0.0;
                                    vertexdist += (tri[l].p[m].x - ii_coor[n][0]) * (tri[l].p[m].x - ii_coor[n][0]);
                                    vertexdist += (tri[l].p[m].y - ii_coor[n][1]) * (tri[l].p[m].y - ii_coor[n][1]);
                                    vertexdist += (tri[l].p[m].z - ii_coor[n][2]) * (tri[l].p[m].z - ii_coor[n][2]);
                                    if (vertexdist < 1e-6)
                                    {
                                        vertexflag = 1;
                                        break;
                                    }
                                }
                                
                                if (vertexflag == 0)
                                {
                                    ii_coor[II][0] = tri[l].p[m].x;
                                    ii_coor[II][1] = tri[l].p[m].y;
                                    ii_coor[II][2] = tri[l].p[m].z;
                                    II ++;
                                }
                            }
                        } 
                    }// loop over k 
                }// loop over j
            }//loop over i
            
            volume = cubecount*dgrid[0]*dgrid[1]*dgrid[2];
            avg_nii = ( (frame_count-1) * avg_nii + II ) / frame_count;
            avg_area = ( (frame_count-1) * avg_area + area ) / frame_count;
            avg_vol = ( (frame_count-1) * avg_vol + volume ) / frame_count;
            
            cout<<time<<"\t"<<II<<"\t"<<avg_nii<<"\t"<<area<<"\t"<<avg_area<<"\t"<<volume<<"\t"<<avg_vol<<endl ;	
            //printf ("%8.1f%8d%10.3f%10.3f%10.3f%10.3f%10.3f\n", time, II, avg_nii, area, avg_area, volume, avg_vol);
            
            if (II > NII)
            cout << "Need assign more space to store II ! at time " << time << endl;
            else
            write_xtc(xdout, NII, step, time, box, ii_coor, prec);
            
            if (frame_count == 1)
            {
                FILE * outfile;
                string filename;
                string head = "PWII_NNN_time";
                string ext = ".gro";
                stringstream ss;
                ss << time;
                filename = head + ss.str() + ext;
                outfile = fopen (filename.c_str(), "w");
                fprintf (outfile, "Instantaneous interfaces by Marching Cube\n");
                fprintf (outfile, "%-10d\n", NII);
                
                for (int i=0; i<NII; i++)
                fprintf (outfile, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", i+1, "TRI", "II", i+1, ii_coor[i][0], ii_coor[i][1], ii_coor[i][2]);
                
                fprintf (outfile, "%10.5f%10.5f%10.5f\n", box[0][0], box[1][1], box[2][2]);
                fclose (outfile);
            }
            
            // deallocate the memory of rhow and rho at each frame
            deallocate3D(rhow, ngrids[0], ngrids[1], ngrids[2]);
            deallocate3D(rho, ngrids[0], ngrids[1], ngrids[2]);
            free (ii_coor);
            ii_coor = NULL;
            
        }// if cut frame 
    }
    
    deallocate3D(rhop, ngrids[0], ngrids[1], ngrids[2]);
    xdrfile_close(xdin);
    xdrfile_close(xdout);
    
    return 0;
}
