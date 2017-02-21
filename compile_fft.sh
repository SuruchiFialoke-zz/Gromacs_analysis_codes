fileroot=$1

fftwlibdir="/opt/local/lib"
fftwincldir="/opt/local/include"

g++ -DCPLUSPLUS -I$fftwincldir -L $fftwlibdir -lfftw3 -o $fileroot.out $fileroot.cpp -Wall

#g++ -o FFTW_2d -I<location of include> -L<location of library> -lfftw3 FFTW2d.cpp
