function [sensorRealDataInterpolada,error] =  interpolateSensorData(sensorRealData,medium, kgrid, coordenadasVizinhos, numeroVizinho, coordenadaSource)
% Através dos dados dos sensorRealData, das informação de quais são os
% sensores efetivamente usados e da fonte emissora. Interpola-se os dados
% formando o sensorRealDataInterpolada

% sensorRealDataInterpolada tem a mesma dimensão de sensorRealData. Mas
% fingimos que conhecemos sensorRealData somente em alguns sensores. A
% partir disso, interpolamos os dados usando coordenadasVizinhos,
% numeroVizinho, velocidade do som no meio e a coordenada da fonte emissora,
% para definir o sensorRealDataInterpolada nos sensores que fingimos não
% conhecer e os que efetivamente são usados.

% Entrada:
% sensorRealData: Dados coletados nos sensores na borda do torso. Nesa
% função fingimos que temos acesso a apenas alguns dos sensores de
% sensorRealData.
% medium: medium no qual o dados serão interpolados, a velocidade do som é
% usada para calcular o tempo de percurso entre a fonte emissora e os
% sensores
% kgrid: kgrid da interpolação, usado para definir as coordenadas x e y das
% fronteiros de cada posição do medium
% coordenadasVizinhos(eixo X ou Y, numero do sensorTomo,Vizinho 1 ou 2 ou EleMesmo, cartesiana ou indice):
% Indica qual é a coordenada cartesiana ou indice dos sensores efetivamente
% usadmos que estão mais proximo ao n-ésimo sensorTomo.
% numeroVizinho(numero do vizinho 1 ou 2, numero do sensorTomo): Indica os
% numeros dos sensores efetivamente usadmos que estão mais proximo ao n-ésimo
% sensorTomo. (Quando zero, indica que o n-ésimo é um sensor que está sendo
% realmente utlizado)
% coordenadaSource(eixo X ou Y,1, cartesiana ou indice): Coordenada da X e
% Y da fonte emissora em coordenadacas cartesianas ou indice. A segunda
% dimensão de coordenadaSource é 1 para manter o padrão.

% Saida:
% sensorRealDataInterpolada: Saida dos dados interpolados
% error(nSensor): mse entre o dado interpolado e o coletado que fingimos
% não existir.

sizeTempo = size(sensorRealData.p,2);
sizeSensor = size(sensorRealData.p,1);
sensorRealDataInterpolada.p = zeros(sizeSensor,sizeTempo);
error = zeros(sizeSensor,1);

for iSensor = 1:sizeSensor
    if(numeroVizinho(1,iSensor) == 0) % O sensor é um sensor efetivamente usado
        sensorRealDataInterpolada.p(iSensor,:) = sensorRealData.p(iSensor,:);
    else
        % Calcula o tempo que o sinal partindo da fonte emissora demora
        % para chegar no sensor, no sensor mais proximo(vizinho 1) e no
        % segundo sensor mais proximo(vizinho 2)
        tempoVizinho1 = calculaTempo(coordenadasVizinhos(:,iSensor,1,1), coordenadaSource(:,1,1),...
            coordenadasVizinhos(:,iSensor,1,2), coordenadaSource(:,1,2),medium, kgrid);
        tempoVizinho2 = calculaTempo(coordenadasVizinhos(:,iSensor,2,1), coordenadaSource(:,1,1),...
            coordenadasVizinhos(:,iSensor,2,2), coordenadaSource(:,1,2),medium, kgrid);
        tempoSensor = calculaTempo(coordenadasVizinhos(:,iSensor,3,1), coordenadaSource(:,1,1),...
            coordenadasVizinhos(:,iSensor,3,2), coordenadaSource(:,1,2),medium, kgrid);
        
        % Diferença de tempo entre o tempo de chegado do sinal no vizinho e
        % no sensor. A diferença pode ser tanto positiva quanto negativa.
        % A diferença é em numero de timeStep
        difVizinho1 = round((tempoVizinho1-tempoSensor)/kgrid.dt);
        difVizinho2 = round((tempoVizinho2-tempoSensor)/kgrid.dt);
        
        % Calcula a distancia entre o sensor e o vizinho 1
        distanciaVizinho1 = sqrt((coordenadasVizinhos(1,iSensor,1,1)-coordenadasVizinhos(1,iSensor,3,1)) * ...
            (coordenadasVizinhos(1,iSensor,1,1)-coordenadasVizinhos(1,iSensor,3,1)) +  ...
            (coordenadasVizinhos(2,iSensor,1,1)-coordenadasVizinhos(2,iSensor,3,1)) ...
            * (coordenadasVizinhos(2,iSensor,1,1)-coordenadasVizinhos(2,iSensor,3,1)));
        
        % Calcula a distancia entre o sensor e o vizinho 1
        distanciaVizinho2 = sqrt((coordenadasVizinhos(1,iSensor,2,1)-coordenadasVizinhos(1,iSensor,3,1)) * ...
            (coordenadasVizinhos(1,iSensor,2,1)-coordenadasVizinhos(1,iSensor,3,1)) + ...
            (coordenadasVizinhos(2,iSensor,2,1)-coordenadasVizinhos(2,iSensor,3,1)) * ...
            (coordenadasVizinhos(2,iSensor,2,1)-coordenadasVizinhos(2,iSensor,3,1)));
        
        % Calcula o peso de cada vizinho de acordo com a distancia. O
        % vizinho mais proximo tem um peso maior.
        proporcaoVizinho1 = distanciaVizinho2 / (distanciaVizinho1 + distanciaVizinho2);
        proporcaoVizinho2 = 1-proporcaoVizinho1;
        
        % Define o range dos indices que sofreram alguma alteraçao, de modo que quando operado
        % não se acesse indices menor do que 1 ou maior do que sizeTempo.
        % Por exemplo, se difVizinho1= -2 e difVizinho1 = 3;
        % rangeVector = 3:(sizeTempo-3)
        rangeVector = (1-min([0 difVizinho1 difVizinho2])):(sizeTempo-max([0 difVizinho1 difVizinho2]));
        
        % Faz a interpolação linear dos dados aquisitados pelos vizinhos
        % pesados pela distancia do vizinho até o sensor, considerando a
        % diferença de tempo entre a chegada do sinal no sensor e os
        % vizinhos.
        sensorRealDataInterpolada.p(iSensor,rangeVector) = ...
            proporcaoVizinho1 * sensorRealData.p(numeroVizinho(1,iSensor), difVizinho1+rangeVector) + ...
            proporcaoVizinho2 * sensorRealData.p(numeroVizinho(2,iSensor), difVizinho2+rangeVector);
        
        % Calcula o erro de interpolação
        error(iSensor) = sqrt(sum( (sensorRealDataInterpolada.p(iSensor,:)-sensorRealData.p(iSensor,:)) .* ...
            (sensorRealDataInterpolada.p(iSensor,:)-sensorRealData.p(iSensor,:))) ./  ...
            sum((sensorRealData.p(iSensor,:) .* sensorRealData.p(iSensor,:))));
        
    end
end

end


function [tempo] = calculaTempo(catDest,catOri,indDest,indOri, medium, kgrid)
% Calcula o tempo necessario para que um onda em linha reta passando por um
% medium heterogenio demora para percorrer entre o origem e o destino.

% Entrado:
% catDest(X ou Y): Coordenada cartesiana do destino
% catOri(X ou Y): Coordenada cartesiana da origem
% indDest(X ou Y): Indice do destino
% indOri(X ou Y): Indice da origem
% medium: medium com a velocida do som em cada gridPoint
% kgrid: kgrid com as coordenadas dos gridPoint

% Saida:
% tempo: Tempo necessario para que um onda em linha reta


theta = abs(atan((catDest(1)-catOri(1))/(catDest(2)-catOri(2))));

if (catDest(1) < catOri(1))
    ix = -1;
elseif (catDest(1) > catOri(1))
    ix = 1;
else
    ix = 0;
end

if (catDest(2) < catOri(2))
    iy = -1;
elseif (catDest(2) > catOri(2))
    iy = 1;
else
    iy = 0;
end

tempo = 0;
if(ix*iy ~= 0)
    fatorX = 1/sin(theta);
    fatorY= 1/cos(theta);
    tgTheta = tan(theta);
    cotgTheta = 1/tgTheta;
    
    while(indDest(1) ~= indOri(1) || indDest(2) ~= indOri(2))
        distanciaParadeX = (kgrid.x_vec(indOri(1))-catOri(1))*(kgrid.x_vec(indOri(1))-catOri(1));
        distanciaParadeY = (kgrid.y_vec(indOri(2))-catOri(2))*(kgrid.y_vec(indOri(2))-catOri(2));
        
        diffDist = distanciaParadeX*fatorX - distanciaParadeY*fatorY;
        if(diffDist < 0)
            tempo = tempo + sqrt(distanciaParadeX*fatorX) / medium.sound_speed(indOri(1),indOri(2));
            indOri(1) = indOri(1)+ix;
            catOri(1) = catOri(1)+ix*distanciaParadeX;
            catOri(2) = catOri(2)+iy*distanciaParadeX*cotgTheta;
        elseif (diffDist > 0)
            tempo = tempo + sqrt(distanciaParadeY*fatorY) / medium.sound_speed(indOri(1),indOri(2));
            indOri(2) = indOri(2)+iy;
            catOri(1) = catOri(1)+ix*distanciaParadeY*tgTheta;
            catOri(2) = catOri(2)+iy*distanciaParadeY;
        else
            tempo = tempo + sqrt(distanciaParadeY*fatorY) / medium.sound_speed(indOri(1),indOri(2));
            indOri(1) = indOri(1)+ix;
            indOri(2) = indOri(2)+iy;
            catOri(1) = catOri(1)+ix*distanciaParadeX;
            catOri(2) = catOri(2)+iy*distanciaParadeY;
        end
    end
elseif(ix == 0)
    while(indDest(2) ~= indOri(2))
        distanciaParadeY = abs(kgrid.y_vec(indOri(2))-catOri(2));
        tempo = tempo + distanciaParadeY / medium.sound_speed(indOri(1),indOri(2));
        indOri(2) = indOri(2)+iy;
        catOri(2) = catOri(2)+iy*distanciaParadeY;
    end
else
    while(indDest(1) ~= indOri(1))
        distanciaParadeX = abs(kgrid.x_vec(indOri(1))-catOri(1));
        tempo = tempo + distanciaParadeX / medium.sound_speed(indOri(1),indOri(2));
        indOri(1) = indOri(1)+ix;
        catOri(1) = catOri(1)+ix*distanciaParadeX;
    end
end

distanciaParadeX = (catOri(1)-catDest(1))*(catOri(1)-catDest(1));
distanciaParadeY = (catOri(2)-catDest(2))*(catOri(2)-catDest(2));
tempo = tempo + sqrt(distanciaParadeX+distanciaParadeY) / medium.sound_speed(indDest(1),indDest(2));

end