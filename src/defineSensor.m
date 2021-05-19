function [sensor] = defineSensor(kgrid, bordaTorso, tipo, arg)
% Dada uma matriz bordaTorso com forma de uma circunferencia deformada 
% seleciona quais pontos dessa borda ser� um sensor segundo o tipo e arg

% Entrada: 
% kgrid: kgrid de entrada usado para saber Nx e Ny
% bordaTorso: uma matriz que define a forma de uma circunferencia deformada
% Tipo: seleciona o tipo de algortimo que ser� usado para definir os
% sensores, tem dois tipos 'nsensor' e 'distSensor'. No algoritmo n sensor
% define-se nSensores igualmente espa�ados na circunferencia deformada. No
% algortimo distSensor define-se a distancia entre os sensores na
% circunferencia deformada.
% arg: argumento do tipo, isto �, escolhe-se ou o numero de sensores ou 
% distancia entre eles.

% Saida:
% sensor: estrutura onde sensor.mask � setado para 1 onde deve existir um
% sensor


sensor.mask = zeros(kgrid.Nx, kgrid.Ny);
tamanhoBorda = sum(sum(bordaTorso));
nsensor =  [];
distSensor = [];

if (strcmp(tipo,'nsensor'))    
    nsensor = arg;
    distSensor  = floor(tamanhoBorda/nsensor);
    rest = mod(tamanhoBorda,nsensor);    
elseif (strcmp(tipo,'distSensor'))
    distSensor = arg;
    rest = 0;
    nsensor = floor(tamanhoBorda/distSensor);
else
    error('Tipo')
end

% Numera a circunferencia deformada de 1 at� N, onde 1 � o primeiro pixel
% na horizontal e vertical e conta no sentido horarios,
i = 1;
while i <= kgrid.Nx
    j=1;
    while j <=kgrid.Ny
        if(bordaTorso(i,j)==1)
            bordaTorso(i,j) = 0;
            bordaContada = zeros(kgrid.Nx, kgrid.Ny);
            bordaContada(i,j) = 1;
            bordaContada = contaBorda(bordaContada,bordaTorso,i,j,1,tamanhoBorda);
            
            % for�a a saida do loop
            bordaTorso(i,j) = 1;
            j = kgrid.Ny+1;
            i = kgrid.Nx+1;
        end
        j = j+1;
    end
    i = i+1;
end

% Define quais numeros da borda deve ter um sensor
coord = zeros(nsensor,1);
coord(1:rest) = ((1:rest)-1)*(distSensor+1)+1;
coord((rest+1):nsensor) = (rest-1)*(distSensor+1)+1+...
    (1:nsensor-rest)*distSensor+1;

% Com a contagem da borda e os numeros que devem possuir um sensor, cria o
% sensor mask
for i=1:length(coord)
    sensor.mask(bordaContada==coord(i)) = 1;
end


end

function  bordaContada = contaBorda(bordaContada,bordaTorso,i,j,n,tamanhoBorda)
% Fun��o iterativa que numera os pixel de uma linha de 1 at� a linha
% terminar. A linha n�o pode ter buracos, nem cruzada e nao pode ser fechada,
%isto �, ela deve ter um inicio e fim bem definidos
% A contagem � feito no sentido horario, verificando o pixel ao norte,
% depois a nordeste e assim por diante.

% Entrada:
% bordaContada: Borda contada at� o momento
% bordaTorso: linha que ainda deve ser contada
% i,j: Ponto atual da leitura
% n: quantidade de pontos j� contados
% tamanhoBorda: quantidade total de pontos da linha que � igual ao numero de 
% pontos a serem contados

if(n==tamanhoBorda)
    %     termina fun��o iterativa
else
    if(bordaTorso(i+1,j) == 1)
        i = i+1;
    elseif(bordaTorso(i+1,j+1) == 1)
        i = i+1;
        j = j+1;
    elseif(bordaTorso(i,j+1) == 1)
        j = j+1;
    elseif(bordaTorso(i-1,j+1) == 1)
        i = i-1;
        j = j+1;
    elseif(bordaTorso(i-1,j) == 1)
        i = i-1;
    elseif(bordaTorso(i-1,j-1) == 1)
        i = i-1;
        j = j-1;
    elseif(bordaTorso(i,j-1) == 1)
        j = j-1;
    elseif(bordaTorso(i+1,j-1) == 1)
        i = i+1;
        j = j-1;
    end
    
    n = n+1;
    bordaContada(i,j) = n;
    bordaTorso(i,j) = 0;
    bordaContada = contaBorda(bordaContada,bordaTorso,i,j,n,tamanhoBorda);
end

end