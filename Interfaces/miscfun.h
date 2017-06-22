#ifndef MISCFUN_H
#define MISCFUN_H

double calphi(double x, double sig, double cutoff)
{
    //double sig = 0.24;
    //double cutoff = 0.7;
    double phic, C;
    double phix;
    
    phic = exp(-cutoff*cutoff/(2.0*sig*sig));
    C = 1.0 / ( pow(2.0*M_PI,0.5) * sig * erf(cutoff / (pow(2.0,0.5) * sig)) - 2.0*cutoff*phic );
    if (abs(x) <= cutoff)
    {
        phix = C * ( exp(-x*x/(2.0*sig*sig)) - phic );
    }
    else
    {
        phix = 0.0;
    }
    return phix;
};

double*** allocate3D(int x, int y, int z)
{
    double*** the_array = new double**[x];
    for ( int i=0; i<x; i++ )
    {
        the_array[i] = new double*[y];
        for ( int j=0; j<y; j++ )
        {
            the_array[i][j] = new double[z];
            for ( int k=0; k<z; k++ )
            {
                the_array[i][j][k] = 0.0;
            }
        }
    }
    return the_array;
};

void deallocate3D(double*** the_array, int x, int y, int z)
{
    for ( int i=0; i<x; i++ )
    {
        for ( int j=0; j<y; j++ )
        {
            delete [] the_array[i][j];
        }
        delete [] the_array[i];
    }
    delete [] the_array;
};

#endif



