function [ensaio] = gerarEnsaio(kgrid, medium, func, tipoBorda, nSource)
% Gera os dados do ensaio (possui hardcoding)
% Um ensaio � composto pelos 'sinais reais' recebidos pelos sensores que
% circundam  o objeto.
% Os sensores e fontes ficam localizados sobre a uma linha fechada que 
% circunda o  objeto. 
% Os sensores s�o continuos nessa linha(Borda)e as nSource fontes emissores
%  est�o igalmente espa�adas.

% Entrada:
% kgrid = kgrid onde ser� feita os ensaios (definidos Nx e Ny)
% medium = medium do objeto
% func = fun��o de excita��o das fontes emissoras
% tipoBorda = Define de que forma os sensores circundan o objeto. O tipo
% pode ser 'Torso' ou 'Quadrado'. No caso 'Torso' os sensores ficam
% localizadas na borda do torso. No caso quadrado, os sensores formam um
% anel quadrado que circunda o objeto. Em ambos os casos os sensores formam
% uma 'circunferencia deformada' continua.
% nSource = Define o numero de fontes emissores igualmente espa�adas ao
% longo da borda.

% Saida:
% ensaio = Struct que armazena os dados recebidos e as condi��es de
% simu��es.
% A ideia � que ao fim do processo de setup a Struct ensaio tenha os
% campos. Sensor_data , kgrid, medium, source, sensorTorso, func 
% (acessiveis por essa fun��o)e os campos  mask, c_indice (n�o acessiveis 
% por essa fun��o)

% Defini��o de parametros temporais
CFL = 0.1;
Tf = 0.4 * 10^-3;

% Defini��o do pml
if(kgrid.Nx >= 256)
    PMLSize = [20,20];
elseif (kgrid.Nx >= 96)
    PMLSize = [20,10];
else
    PMLSize = [10,5];
end

% Define bordaTorso segundo o tipoBorda
if(strcmp(tipoBorda,'Torso'))
    [bordaTorso,~] = takeBoderTorso(kgrid);
elseif(strcmp(tipoBorda,'Quadrado'))
    p1 = [21/96,11/96];
    p2 = [21/96,86/96];
    p3 = [76/96,86/96];
    p4 = [76/96,11/96];
    Nsc = [kgrid.Nx,kgrid.Ny];
    bordaTorso = makeLine(kgrid.Nx, kgrid.Ny, ceil(Nsc.*p1),ceil(Nsc.*p2));
    bordaTorso = bordaTorso | makeLine(kgrid.Nx,kgrid.Ny, ceil(Nsc.*p2),ceil(Nsc.*p3));
    bordaTorso = bordaTorso | makeLine(kgrid.Nx,kgrid.Ny, ceil(Nsc.*p4),ceil(Nsc.*p3));
    bordaTorso = bordaTorso | makeLine(kgrid.Nx,kgrid.Ny, ceil(Nsc.*p4),ceil(Nsc.*p1));
else
    disp('tipoBorda n�o identificado');
end

%  Define os parametros temporais do tipoBorda
[kgrid.t_array,~] = makeTime(kgrid, medium.sound_speed, CFL, Tf);


% define os sensores na borda do torso
[sensorTorso] = defineSensor(kgrid,bordaTorso,'distSensor',1);
sensorTorso.record = {'p'};

%Interpola a fun��o de excita��o conforme o vetor de tempo de kgrid
F = griddedInterpolant(0:Tf/length(func):Tf-Tf/length(func),func, 'cubic');
fun = F(kgrid.t_array);

% A fun��o j� possui espectro adequado para ser simulado
% fun = filterTimeSeries(kgrid, medium, fun,'PlotSignals',true,'PlotSpectrums',true);

% define as fontes emissores
[source] = defineSouce(kgrid,bordaTorso,nSource,fun); 

% Realiza nSource simula��es da propaga��o direta da onda, armazenando a
% informa��o em sensor_data
for i = 1:nSource
    sourceAux.p_mask = source.p_mask;
    sourceAux.p = zeros(size(source.p));
    sourceAux.p(i,:) = source.p(i,:);    
    sensor_data(i) = kspaceFirstOrder2D(kgrid, medium, sourceAux, sensorTorso,'PMLSize',PMLSize,'DataCast', 'single','PlotSim' ,false);
end

%Atribui os campos viziveis pertirnente ao ensaio 
ensaio.sensor_data = sensor_data;
ensaio.kgrid = kgrid;
ensaio.medium = medium;
ensaio.source = source;
ensaio.sensorTorso = sensorTorso;
ensaio.func = func;

end