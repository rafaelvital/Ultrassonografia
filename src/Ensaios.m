% Registro de ensaio

load('Func.mat');
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(96,3);
% 
% [kgrid.t_array,dt] = makeTime(kgrid, medium.sound_speed, 0.1, 0.4*10^-3);
% 
% cria os source
% fun = zeros(1,length(kgrid.t_array));
% fun(1: floor(10^-5/dt)) = 1;
% fun = 50 * sin(2*pi*kgrid.t_array*10^5) .* fun;
% fun = filterTimeSeries(kgrid, medium, fun,'ZeroPhase',true);

% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(96,1);
% [ensaio17] = gerarEnsaio(kgrid, medium,func,'Quadrado');
% ensaio17.medium = medium;
% ensaio17.mask = mask;
% ensaio17.c_indice = c_indice;
% ensaio17.func = func;

[kgrid, medium, mask, c_indice] = definirImpulsos(192);
[imp] = gerarEnsaio(kgrid, medium,func,'Quadrado');
imp.medium = medium;
imp.mask = mask;
imp.c_indice = c_indice;
imp.func = func;

save('imp2' ,'imp');

% Ensaio 1
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,3);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 20;
% I = 3
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96, CFL 0.3


% Ensaio 2
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,3);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 20;
% I = 3
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% funcZero = N96, CFL 0.1, faseZero

% Ensaio 3
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,3);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 20;
% I = 3
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% funcZero = N96
% Borda Quadrada

% Ensaio 4
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,3);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 60;
% I = 3
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% funcZero = N96
% Borda Quadrada

% Ensaio 5
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(96,3);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 96
% nSource = 20;
% I = 3
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 5;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Ensaio 6
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(96,3);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 96
% nSource = 60;
% I = 3
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 5;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Ensaio 7
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,3);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 96
% nSource = 60;
% I = 3
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 5;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Ensaio 8
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,3);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 60;
% I = 3
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Ensaio 9
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(96,3);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 96
% nSource = 60;
% I = 3
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada
% 
% Ensaio 10
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,6);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 60;
% I = 6
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Ensaio 11
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,15);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 60;
% I = 15
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Ensaio 12
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,6);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 60;
% I = 6
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda dilatada

% Ensaio 13
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,9);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 60;
% I = 9
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Ensaio 14
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,6);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 40;
% I = 6
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Ensaio 15
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,6);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 20;
% I = 6
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Ensaio 16
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,6);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 10;
% I = 6
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Ensaio 17
% [kgrid, medium, mask, c_indice] = poucoHeterogeoC(256,6);
% [ensaio1] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 30;
% I = 6
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada

% Imp 1
% [kgrid, medium, mask, c_indice] = definirInpulso(256);
% [inpt] = gerarEnsaio(kgrid, medium);
% N = 256
% nSource = 30;
% I = 1400
% CFL = 0.1;
% Tf = 0.4*10^-3;
% step = 20;
% frequenciaOnda = 10^5;
% intensidadeOnda = 50;
% func = N96
% Borda Quadrada
