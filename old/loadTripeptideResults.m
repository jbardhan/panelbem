fftsvd_bem = struct('A',0,'B',0,'C',0);
fftsvd_bem.A = load('../meshes/tripeptide/A_res1.txt');
fftsvd_bem.B = load('../meshes/tripeptide/B_res1.txt');
fftsvd_bem.C = load('../meshes/tripeptide/C_res1.txt');
fftsvd_bem.A = fftsvd_bem.A';
fftsvd_bem.B = fftsvd_bem.B';
fftsvd_bem.C = fftsvd_bem.C';
% the files are saved out a column at a time... which means a row
% at a time in the text file, thus the transpose!