
#export problem=" -nx 64 -nt 256 "
export problem=" -nx 256 -nt 1024 "
export PDEParams=" -eps 0 -prob 1 -a 0.24 -xL 2 -xR 1 -st 0 -skip 0 -tfinal 1.9"
export BraidParams=" -nu 2 -ml 12 -fmg 1 -tol 1e-10 -asc 1 -sc 1"

echo "Parameters are  "
echo $PDEParams
echo $BraidParams
echo "Problem size: $problem1"
echo " "
echo " "
echo "Setup:  Spatial coarsen every OTHER level, First in time and then space, add more relaxation and more FMG"

# Plots for 1 iteration
rm *out*
ex-burgers $problem $PDEParams $BraidParams -mi 1
python viz-burgers-for-powerpoint.py
mv burgers_step256.pdf  burgers_step256_iter1.pdf
mv burgers_step128.pdf burgers_step128_iter1.pdf
mv burgers_step0.pdf burgers_step0_iter1.pdf
mv burgers_3D.png burgers_3D_iter1.png
mv burgers_imshow.pdf burgers_imshow_iter1.pdf


# Plots for 2 iteration
rm *out*
ex-burgers $problem $PDEParams $BraidParams -mi 2
python viz-burgers-for-powerpoint.py
mv burgers_step256.pdf  burgers_step256_iter2.pdf
mv burgers_step128.pdf burgers_step128_iter2.pdf
mv burgers_step0.pdf burgers_step0_iter2.pdf
mv burgers_3D.png burgers_3D_iter2.png
mv burgers_imshow.pdf burgers_imshow_iter2.pdf


# Plots for 3 iteration
rm *out*
ex-burgers $problem $PDEParams $BraidParams -mi 3
python viz-burgers-for-powerpoint.py
mv burgers_step256.pdf  burgers_step256_iter3.pdf
mv burgers_step128.pdf burgers_step128_iter3.pdf
mv burgers_step0.pdf burgers_step0_iter3.pdf
mv burgers_3D.png burgers_3D_iter3.png
mv burgers_imshow.pdf burgers_imshow_iter3.pdf



# Plots for 4 iteration
rm *out*
ex-burgers $problem $PDEParams $BraidParams -mi 4
python viz-burgers-for-powerpoint.py
mv burgers_step256.pdf  burgers_step256_iter4.pdf
mv burgers_step128.pdf burgers_step128_iter4.pdf
mv burgers_step0.pdf burgers_step0_iter4.pdf
mv burgers_3D.png burgers_3D_iter4.png
mv burgers_imshow.pdf burgers_imshow_iter4.pdf





