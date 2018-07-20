# In the context of a scaling study...
#
# Possible tests:
#  Other dt/dx, possibly smaller ones
#
#
# Possible implementations:
#  Implement Jacobi, and compare Jacobi and FCF
#  Implement a "Galerkin" coarse-grid
#  Implement a higher-order time-stepper (but not in space)
#
# $$ egrep "iterations         |Problem size:|Setup:" test.out1 
#
# $$ ex-burgers -nu 1 -nx 128 -nt 128  -eps 0 -prob 0 -a 1 -xL 2 -xR 1 -st 1 -mi 50 -sc 1 -ml 25 -fmg 1 -asc 1

export problem1=" -nx 128  -nt 128 "
export problem2=" -nx 256  -nt 256 "
export problem3=" -nx 512  -nt 512 "
export problem4=" -nx 1024 -nt 1024 "
export problem5=" -nx 2048 -nt 2048 "

export problem6="-asdfasfd " 
#export problem6=" -nx 4096 -nt 4096 "

# a = 0.1136 matches his 0.05 (a dt/dx)
export PDEParams=" -eps 0 -prob 0 -a 0.1136 -xL 2 -xR 2 -st 0 -skip 0 -tfinal 2.2"
export BraidParams=" -nu 1  -mi 50 -ml 45 -fmg 1 -tol 1e-10 "

echo "Parameters are  "
echo $PDEParams
echo $BraidParams
echo " "
echo " "

echo "---------------------------------------------------------------------- "
echo "---------------------------------------------------------------------- "
echo "Setup:  Two-level, No spatial coarsening tests"
echo " "; echo " "; echo "Problem size: $problem1"
ex-burgers $problem1 $PDEParams $BraidParams -ml 2
echo " "; echo " "; echo "Problem size: $problem2"
ex-burgers $problem2 $PDEParams $BraidParams -ml 2
echo " "; echo " "; echo "Problem size: $problem3"
ex-burgers $problem3 $PDEParams $BraidParams -ml 2
echo " "; echo " "; echo "Problem size: $problem4"
ex-burgers $problem4 $PDEParams $BraidParams -ml 2
echo " "; echo " "; echo "Problem size: $problem5"
ex-burgers $problem5 $PDEParams $BraidParams -ml 2
echo " "; echo " "; echo "Problem size: $problem6"
ex-burgers $problem6 $PDEParams $BraidParams -ml 2

echo "---------------------------------------------------------------------- "
echo "---------------------------------------------------------------------- "
echo " "
echo "Setup:  Spatial coarsen on every level, append -sc 1"
echo " "; echo " "; echo "Problem size: $problem1"
ex-burgers $problem1 $PDEParams $BraidParams -sc 1 
echo " "; echo " "; echo "Problem size: $problem2"
ex-burgers $problem2 $PDEParams $BraidParams -sc 1 
echo " "; echo " "; echo "Problem size: $problem3"
ex-burgers $problem3 $PDEParams $BraidParams -sc 1 
echo " "; echo " "; echo "Problem size: $problem4"
ex-burgers $problem4 $PDEParams $BraidParams -sc 1 
echo " "; echo " "; echo "Problem size: $problem5"
ex-burgers $problem5 $PDEParams $BraidParams -sc 1 
echo " "; echo " "; echo "Problem size: $problem6"
ex-burgers $problem6 $PDEParams $BraidParams -sc 1 


echo "---------------------------------------------------------------------- "
echo "---------------------------------------------------------------------- "
echo " "
echo "Setup:  Spatial coarsen every OTHER level, First in time and then space, append -sc 1, -asc 1"
echo " "; echo " "; echo "Problem size: $problem1"
ex-burgers $problem1 $PDEParams $BraidParams -sc 1 -asc 1
echo " "; echo " "; echo "Problem size: $problem2"
ex-burgers $problem2 $PDEParams $BraidParams -sc 1 -asc 1
echo " "; echo " "; echo "Problem size: $problem3"
ex-burgers $problem3 $PDEParams $BraidParams -sc 1 -asc 1
echo " "; echo " "; echo "Problem size: $problem4"
ex-burgers $problem4 $PDEParams $BraidParams -sc 1 -asc 1
echo " "; echo " "; echo "Problem size: $problem5"
ex-burgers $problem5 $PDEParams $BraidParams -sc 1 -asc 1
echo " "; echo " "; echo "Problem size: $problem6"
ex-burgers $problem6 $PDEParams $BraidParams -sc 1 -asc 1


echo "---------------------------------------------------------------------- "
echo "---------------------------------------------------------------------- "
echo " "
echo "Setup:  Spatial coarsen every OTHER level, First in time and then space, add more relaxation on coarse levels, append  -sc 1 -asc 1 -nu 2 -nu0 1 "
echo " "; echo " "; echo "Problem size: $problem1"
ex-burgers $problem1 $PDEParams $BraidParams -sc 1 -asc 1 -nu 2 -nu0 1
echo " "; echo " "; echo "Problem size: $problem2"
ex-burgers $problem2 $PDEParams $BraidParams -sc 1 -asc 1 -nu 2 -nu0 1
echo " "; echo " "; echo "Problem size: $problem3"
ex-burgers $problem3 $PDEParams $BraidParams -sc 1 -asc 1 -nu 2 -nu0 1
echo " "; echo " "; echo "Problem size: $problem4"
ex-burgers $problem4 $PDEParams $BraidParams -sc 1 -asc 1 -nu 2 -nu0 1
echo " "; echo " "; echo "Problem size: $problem5"
ex-burgers $problem5 $PDEParams $BraidParams -sc 1 -asc 1 -nu 2 -nu0 1
echo " "; echo " "; echo "Problem size: $problem6"
ex-burgers $problem6 $PDEParams $BraidParams -sc 1 -asc 1 -nu 2 -nu0 1


echo "---------------------------------------------------------------------- "
echo "---------------------------------------------------------------------- "
echo " "
echo "Setup:  Spatial coarsen every OTHER level, First in time and then space, add more powerful FMG, nfmg_Vcyc=2, append  -sc 1 -asc 1 -fmg 2 "
echo " "; echo " "; echo "Problem size: $problem1"
ex-burgers $problem1 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2
echo " "; echo " "; echo "Problem size: $problem2"
ex-burgers $problem2 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2
echo " "; echo " "; echo "Problem size: $problem3"
ex-burgers $problem3 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2
echo " "; echo " "; echo "Problem size: $problem4"
ex-burgers $problem4 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2
echo " "; echo " "; echo "Problem size: $problem5"
ex-burgers $problem5 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2
echo " "; echo " "; echo "Problem size: $problem6"
ex-burgers $problem6 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2

echo "---------------------------------------------------------------------- "
echo "---------------------------------------------------------------------- "
echo " "
echo "Setup:  Spatial coarsen every OTHER level, First in time and then space, add more powerful FMG, nfmg_Vcyc=2 and more relaxation on coarse levels, append  -sc 1 -asc 1 -fmg 2 -nu 2 -nu0 1 "
echo " "; echo " "; echo "Problem size: $problem1"
ex-burgers $problem1 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2 -nu 2 -nu0 1
echo " "; echo " "; echo "Problem size: $problem2"
ex-burgers $problem2 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2 -nu 2 -nu0 1
echo " "; echo " "; echo "Problem size: $problem3"
ex-burgers $problem3 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2 -nu 2 -nu0 1
echo " "; echo " "; echo "Problem size: $problem4"
ex-burgers $problem4 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2 -nu 2 -nu0 1
echo " "; echo " "; echo "Problem size: $problem5"
ex-burgers $problem5 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2 -nu 2 -nu0 1
echo " "; echo " "; echo "Problem size: $problem6"
ex-burgers $problem6 $PDEParams $BraidParams -sc 1 -asc 1 -fmg 2 -nu 2 -nu0 1


echo "---------------------------------------------------------------------- "
echo "---------------------------------------------------------------------- "
echo " "
echo "Setup:  Two level, Spatial coarsen every OTHER level, First in space and then time, append -sc 1, -asc 2"
echo " "; echo " "; echo "Problem size: $problem1"
ex-burgers $problem1 $PDEParams $BraidParams -sc 1 -asc 2 -ml 2
echo " "; echo " "; echo "Problem size: $problem2"
ex-burgers $problem2 $PDEParams $BraidParams -sc 1 -asc 2 -ml 2
echo " "; echo " "; echo "Problem size: $problem3"
ex-burgers $problem3 $PDEParams $BraidParams -sc 1 -asc 2 -ml 2
echo " "; echo " "; echo "Problem size: $problem4"
ex-burgers $problem4 $PDEParams $BraidParams -sc 1 -asc 2 -ml 2
echo " "; echo " "; echo "Problem size: $problem5"
ex-burgers $problem5 $PDEParams $BraidParams -sc 1 -asc 2 -ml 2
echo " "; echo " "; echo "Problem size: $problem6"
ex-burgers $problem6 $PDEParams $BraidParams -sc 1 -asc 2 -ml 2


echo "---------------------------------------------------------------------- "
echo "---------------------------------------------------------------------- "
echo " "
echo "Setup:  Spatial coarsen every OTHER level, First in space and then time, append -sc 1, -asc 2"
echo " "; echo " "; echo "Problem size: $problem1"
ex-burgers $problem1 $PDEParams $BraidParams -sc 1 -asc 2 -ml 12
echo " "; echo " "; echo "Problem size: $problem2"
ex-burgers $problem2 $PDEParams $BraidParams -sc 1 -asc 2 -ml 14
echo " "; echo " "; echo "Problem size: $problem3"
ex-burgers $problem3 $PDEParams $BraidParams -sc 1 -asc 2 -ml 16
echo " "; echo " "; echo "Problem size: $problem4"
ex-burgers $problem4 $PDEParams $BraidParams -sc 1 -asc 2 -ml 18
echo " "; echo " "; echo "Problem size: $problem5"
ex-burgers $problem5 $PDEParams $BraidParams -sc 1 -asc 2 -ml 20
echo " "; echo " "; echo "Problem size: $problem6"
ex-burgers $problem6 $PDEParams $BraidParams -sc 1 -asc 2 -ml 22

