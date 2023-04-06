#!/bin/sh
cd ~/researchpaper/paper6/fortrancode/
gfortran fortran_code1.f90 -L/usr/local/lib -llapack -lrefblas -o code1.out
./code1.out

gfortran fortran_code2.f90 -L/usr/local/lib -llapack -lrefblas -o code2.out
./code2.out

gfortran fortran_code3.f90 -L/usr/local/lib -llapack -lrefblas -o code3.out
./code3.out

gfortran fortran_code4.f90 -L/usr/local/lib -llapack -lrefblas -o code4.out
./code4.out

cd ~/researchpaper/paper6/fortrancode2/
gfortran fortrancode_solution21.f90 -L/usr/local/lib -llapack -lrefblas -o code1.out
./code1.out

gfortran fortrancode_solution22.f90 -L/usr/local/lib -llapack -lrefblas -o code2.out
./code2.out

gfortran fortrancode_solution23.f90 -L/usr/local/lib -llapack -lrefblas -o code3.out
./code3.out

gfortran fortrancode_solution24.f90 -L/usr/local/lib -llapack -lrefblas -o code4.out
./code4.out


cd ~/researchpaper/paper6/data/case_a/
gnuplot ./case_a_s11.gp
gnuplot ./case_a_s22.gp
gnuplot ./case_a_t12.gp
gnuplot ./case_a_u.gp
gnuplot ./case_a_v.gp
cd ~/researchpaper/paper6/texfig/case_a_s11/
epstopdf ./case_a_s11-inc.eps
pdflatex ./case_a_s11.tex
cp ./case_a_s11.pdf ~/researchpaper/paper6/arXiv/fig/fig2c.pdf
cd ~/researchpaper/paper6/texfig/case_a_s22/
epstopdf ./case_a_s22-inc.eps
pdflatex ./case_a_s22.tex
cp ./case_a_s22.pdf ~/researchpaper/paper6/arXiv/fig/fig2b.pdf
cd ~/researchpaper/paper6/texfig/case_a_t12/
epstopdf ./case_a_t12-inc.eps
pdflatex ./case_a_t12.tex
cp ./case_a_t12.pdf ~/researchpaper/paper6/arXiv/fig/fig2d.pdf
cd ~/researchpaper/paper6/texfig/case_a_u/
epstopdf ./case_a_u-inc.eps
pdflatex ./case_a_u.tex
cp ./case_a_u.pdf ~/researchpaper/paper6/arXiv/fig/fig2e.pdf
cd ~/researchpaper/paper6/texfig/case_a_v/
epstopdf ./case_a_v-inc.eps
pdflatex ./case_a_v.tex
cp ./case_a_v.pdf ~/researchpaper/paper6/arXiv/fig/fig2f.pdf

cd ~/researchpaper/paper6/data/case_b/
gnuplot ./case_b_s11.gp
gnuplot ./case_b_s22.gp
gnuplot ./case_b_t12.gp
gnuplot ./case_b_u.gp
gnuplot ./case_b_v.gp
cd ~/researchpaper/paper6/texfig/case_b_s11/
epstopdf ./case_b_s11-inc.eps
pdflatex ./case_b_s11.tex
cp ./case_b_s11.pdf ~/researchpaper/paper6/arXiv/fig/fig3c.pdf
cd ~/researchpaper/paper6/texfig/case_b_s22/
epstopdf ./case_b_s22-inc.eps
pdflatex ./case_b_s22.tex
cp ./case_b_s22.pdf ~/researchpaper/paper6/arXiv/fig/fig3b.pdf
cd ~/researchpaper/paper6/texfig/case_b_t12/
epstopdf ./case_b_t12-inc.eps
pdflatex ./case_b_t12.tex
cp ./case_b_t12.pdf ~/researchpaper/paper6/arXiv/fig/fig3d.pdf
cd ~/researchpaper/paper6/texfig/case_b_u/
epstopdf ./case_b_u-inc.eps
pdflatex ./case_b_u.tex
cp ./case_b_u.pdf ~/researchpaper/paper6/arXiv/fig/fig3e.pdf
cd ~/researchpaper/paper6/texfig/case_b_v/
epstopdf ./case_b_v-inc.eps
pdflatex ./case_b_v.tex
cp ./case_b_v.pdf ~/researchpaper/paper6/arXiv/fig/fig3f.pdf

cd ~/researchpaper/paper6/data/case_c/
gnuplot ./case_c_s11.gp
gnuplot ./case_c_s22.gp
gnuplot ./case_c_t12.gp
gnuplot ./case_c_u.gp
gnuplot ./case_c_v.gp
cd ~/researchpaper/paper6/texfig/case_c_s11/
epstopdf ./case_c_s11-inc.eps
pdflatex ./case_c_s11.tex
cp ./case_c_s11.pdf ~/researchpaper/paper6/arXiv/fig/fig4c.pdf
cd ~/researchpaper/paper6/texfig/case_c_s22/
epstopdf ./case_c_s22-inc.eps
pdflatex ./case_c_s22.tex
cp ./case_c_s22.pdf ~/researchpaper/paper6/arXiv/fig/fig4b.pdf
cd ~/researchpaper/paper6/texfig/case_c_t12/
epstopdf ./case_c_t12-inc.eps
pdflatex ./case_c_t12.tex
cp ./case_c_t12.pdf ~/researchpaper/paper6/arXiv/fig/fig4d.pdf
cd ~/researchpaper/paper6/texfig/case_c_u/
epstopdf ./case_c_u-inc.eps
pdflatex ./case_c_u.tex
cp ./case_c_u.pdf ~/researchpaper/paper6/arXiv/fig/fig4e.pdf
cd ~/researchpaper/paper6/texfig/case_c_v/
epstopdf ./case_c_v-inc.eps
pdflatex ./case_c_v.tex
cp ./case_c_v.pdf ~/researchpaper/paper6/arXiv/fig/fig4f.pdf

cd ~/researchpaper/paper6/data/case_d/
gnuplot ./case_d_s11.gp
gnuplot ./case_d_s22.gp
gnuplot ./case_d_t12.gp
gnuplot ./case_d_u.gp
gnuplot ./case_d_v.gp
cd ~/researchpaper/paper6/texfig/case_d_s11/
epstopdf ./case_d_s11-inc.eps
pdflatex ./case_d_s11.tex
cp ./case_d_s11.pdf ~/researchpaper/paper6/arXiv/fig/fig5c.pdf
cd ~/researchpaper/paper6/texfig/case_d_s22/
epstopdf ./case_d_s22-inc.eps
pdflatex ./case_d_s22.tex
cp ./case_d_s22.pdf ~/researchpaper/paper6/arXiv/fig/fig5b.pdf
cd ~/researchpaper/paper6/texfig/case_d_t12/
epstopdf ./case_d_t12-inc.eps
pdflatex ./case_d_t12.tex
cp ./case_d_t12.pdf ~/researchpaper/paper6/arXiv/fig/fig5d.pdf
cd ~/researchpaper/paper6/texfig/case_d_u/
epstopdf ./case_d_u-inc.eps
pdflatex ./case_d_u.tex
cp ./case_d_u.pdf ~/researchpaper/paper6/arXiv/fig/fig5e.pdf
cd ~/researchpaper/paper6/texfig/case_d_v/
epstopdf ./case_d_v-inc.eps
pdflatex ./case_d_v.tex
cp ./case_d_v.pdf ~/researchpaper/paper6/arXiv/fig/fig5f.pdf


cd ~/researchpaper/paper6/data/ck12/
gnuplot ./ck12_case_a.gp
gnuplot ./ck12_case_b.gp
gnuplot ./ck12_case_c.gp
gnuplot ./ck12_case_d.gp

cd ~/researchpaper/paper6/texfig/ck12_case_a/
epstopdf ./ck12_case_a-inc.eps
pdflatex ./ck12_case_a.tex
cp ./ck12_case_a.pdf ~/researchpaper/paper6/arXiv/fig/fig6a.pdf

cd ~/researchpaper/paper6/texfig/ck12_case_b/
epstopdf ./ck12_case_b-inc.eps
pdflatex ./ck12_case_b.tex
cp ./ck12_case_b.pdf ~/researchpaper/paper6/arXiv/fig/fig6b.pdf

cd ~/researchpaper/paper6/texfig/ck12_case_c/
epstopdf ./ck12_case_c-inc.eps
pdflatex ./ck12_case_c.tex
cp ./ck12_case_c.pdf ~/researchpaper/paper6/arXiv/fig/fig6c.pdf

cd ~/researchpaper/paper6/texfig/ck12_case_d/
epstopdf ./ck12_case_d-inc.eps
pdflatex ./ck12_case_d.tex
cp ./ck12_case_d.pdf ~/researchpaper/paper6/arXiv/fig/fig6d.pdf


cd ~/researchpaper/paper6/data/dk/
gnuplot ./dk_case_a.gp
gnuplot ./dk_case_b.gp
gnuplot ./dk_case_c.gp
gnuplot ./dk_case_d.gp

cd ~/researchpaper/paper6/texfig/dk_case_a/
epstopdf ./dk_case_a-inc.eps
pdflatex ./dk_case_a.tex
cp ./dk_case_a.pdf ~/researchpaper/paper6/arXiv/fig/fig7a.pdf

cd ~/researchpaper/paper6/texfig/dk_case_b/
epstopdf ./dk_case_b-inc.eps
pdflatex ./dk_case_b.tex
cp ./dk_case_b.pdf ~/researchpaper/paper6/arXiv/fig/fig7b.pdf

cd ~/researchpaper/paper6/texfig/dk_case_c/
epstopdf ./dk_case_c-inc.eps
pdflatex ./dk_case_c.tex
cp ./dk_case_c.pdf ~/researchpaper/paper6/arXiv/fig/fig7c.pdf

cd ~/researchpaper/paper6/texfig/dk_case_d/
epstopdf ./dk_case_d-inc.eps
pdflatex ./dk_case_d.tex
cp ./dk_case_d.pdf ~/researchpaper/paper6/arXiv/fig/fig7d.pdf

cd ~/researchpaper/paper6/arXiv/
pdflatex ./arXiv_style.tex
