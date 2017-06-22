
module load gromacs-gcc

fileroot=$1

#xdrlibdir="/opt/gromacs/4.5.5/lib"
#xdrincldir="/opt/gromacs/4.5.5/include/xdrfile"

xdrlibdir="/usr/lib64"
xdrincldir="/usr/include/xdrfile"

g++ -DCPLUSPLUS -I $xdrincldir -L $xdrlibdir -lxdrfile -lm -o $fileroot.out $fileroot.cpp -Wall

