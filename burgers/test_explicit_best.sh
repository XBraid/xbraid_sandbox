# For TODO look at test_implicit.sh 
#
# $$ egrep "iterations         |Problem size:|Setup:" test.out1 
#
# $$ ex-burgers -nu 1 -nx 64 -nt 64 -eps 0 -prob 0 -xL 2 -xR 1 -st 0 -mi 50 -sc 1 -ml 25 -fmg 1 -asc 1 -a 0.24
#
# You can actually run with a=0.48, so long as you don't exactly solve divergent cfl=1.0 coarse-grid. 


export problem1=" -nx 32  -nt 128 "
export problem2=" -nx 64  -nt 256 "
export problem3=" -nx 128 -nt 512 "
export problem4=" -nx 256 -nt 1024 "
export problem5=" -nx 512 -nt 2048 "
export problem6=" -nx 1024 -nt 4096 "
export problem7=" -nx 2048 -nt 8192 "

export PDEParams=" -eps 0 -prob 1 -a 0.24 -xL 2 -xR 1 -st 0 -skip 0 -tfinal 2.2"
export BraidParams=" -nu 2  -mi 50 -ml 45 -fmg 1 -tol 1e-10 "

echo "Parameters are  "
echo $PDEParams
echo $BraidParams
echo " "
echo " "

echo "---------------------------------------------------------------------- "
echo "---------------------------------------------------------------------- "
echo " "
echo "Setup:  Spatial coarsen every OTHER level, First in time and then space, add more relaxation and more FMG, append  -sc 1 -asc 1 -nu 2 -fmg 2 "
echo " "; echo " "; echo "Problem size: $problem1"
ex-burgers $problem1 $PDEParams $BraidParams -ml 10 -sc 1 -asc 1 -nu 2 -fmg 2
echo " "; echo " "; echo "Problem size: $problem2"
ex-burgers $problem2 $PDEParams $BraidParams -ml 12 -sc 1 -asc 1 -nu 2 -fmg 2
echo " "; echo " "; echo "Problem size: $problem3"
ex-burgers $problem3 $PDEParams $BraidParams -ml 14 -sc 1 -asc 1 -nu 2 -fmg 2
echo " "; echo " "; echo "Problem size: $problem4"
ex-burgers $problem4 $PDEParams $BraidParams -ml 16 -sc 1 -asc 1 -nu 2 -fmg 2
echo " "; echo " "; echo "Problem size: $problem5"
ex-burgers $problem5 $PDEParams $BraidParams -ml 18 -sc 1 -asc 1 -nu 2 -fmg 2
echo " "; echo " "; echo "Problem size: $problem6"
ex-burgers $problem6 $PDEParams $BraidParams -ml 20 -sc 1 -asc 1 -nu 2 -fmg 2
echo " "; echo " "; echo "Problem size: $problem7"
ex-burgers $problem7 $PDEParams $BraidParams -ml 22 -sc 1 -asc 1 -nu 2 -fmg 2


