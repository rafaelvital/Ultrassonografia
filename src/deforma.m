function forma = deforma(forma,prob,iteracao,deslocamento)
% Dado uma forma (matriz binaria) a deforma fazendo sucessivas dilatações e
% erosões na borda da forma de modo randomico.

% Entrada:
% foma = forma a ser deformada
% prob = Probabilidade de deformar os pontos da borda
% iteração = Numero de interaçõs, lembrando que uma iteração é composta de
% uma deformação e uma erosão
% deslocamento = Defino o intervalo maximo de deslocamento que a forma pode 
% sofrer no sentido vertical e horizontal. Lembrando que será sorteado um 
% numero entre [-deslocamento, deslocamento] para deslocar a forma

% Saida:
% forma = forma deformada

%% deforma a forma

se = strel('diamond',1);
for i=1:iteracao
    
    %identifica a borda interna da forma    
    formaErode = imerode(forma,se);
    border = forma-formaErode;
    
    %encontra o indice das bordas  e sorteia nPontosNaBorda numeros aleatorios     
    indice = find(border == 1);
    n = length(indice);
    rNumber = rand(n, 1);
    
    % Define quais pontos da borda vao mudar de 1 para 0;
    forma(border == 1) = (rNumber >  prob);
      
    %Na nova forma, identifica a borda  extena
    formaDil = imdilate(forma,se);
    border = formaDil - forma;
    
    %encontra o indice das bordas  e sorteia nPontosNaBorda numeros aleatorios     
    indice = find(border == 1);
    n = length(indice);
    rNumber = rand(n, 1);
    
    %encontra o indice da bordas externa que vao ser adicionas na forma,
    forma(border == 1) = (rNumber >  1-prob);
    
    % A cada 10 iterações faz um close para tapar os possiveis buracos
    if (mod(iteracao,10) == 0)
        forma = imclose(forma,se);
    end
    
end

%close para tapar os possiveis buracos
forma = imclose(forma,se);


%% Deslovamento

% Escolhe o deslocamento (em numero de pontos) que a forma sofrerá na
% vertical e horizontal
des = randi(2*deslocamento+1,1,2) - deslocamento - 1;

% Desloca a forma
forma2 = forma;
for i= (1+abs(des(1))): (size(forma,1)-abs(des(1)))
    for j=(1+abs(des(2))):(size(forma,2)-abs(des(2)))
        forma(i,j) = forma2(i+des(1),j+des(2));
    end
end


end