function [source] = defineSouce(kgrid,bordaTorso,nsource,fun)
% Dada uma matriz bordaTorso com forma de uma circunferencia deformada
% seleciona quais pontos dessa borda teram uma fonte emissore que emite o
% sinal fun

%Entrada:
% kgrid: kgrid de entrada usado para saber Nx e Ny
% bordaTorso: uma matriz que define a forma de uma circunferencia deformada
% nsource: seleciona o numero de fontes emissores igualmente espaçadas na
% borda do torso
% fun: Função de excitaçao de cada fonte

% define a mascara de u da mesma forma que define os sensores
[sensor] = defineSensor(kgrid,bordaTorso,'nsensor',nsource);
source.p_mask = sensor.mask;
source.p = repmat(fun,nsource,1);

end