function [medium,mask,c_indice] = aproximacaoInicialPoucoHeterogeneoC(medium, mask, iterationForma,deslocamento,c_indice,cIntensity)
% Cria uma medium aproximad do objeto tomografico.
% Essa aproxima��o ocorre alterando-se a forma das regioes do objeto, isto
% �, alterando-se o dominio do pulm�o, cora��o, coluna e meio.
% Ou alterando-se o velocidade de propaga��o do som nessas regi�es

% Entrada:
% medium: medium do objeto tomografado que deseja-se criar uma aproxima��o
% mask: Marcara das regi�es (pulm�o, cora��o, coluna e meio) que forma o objeto
% iterationForma: quantas deforma��es de forma ser�o executadas
% deslocamento: Amplitude maxima do deslocamento da forma
% c_indice: Velocidade do som das regi�es que formam o objeto
% cIntensity: intensidade da varia��o do c_indice que o medium aproximado
% ter� com rela��o ao medium do objeto tomografado

% Saida:
% medium: Medium aproximado
% mask: Mascar� da regioes do objeto aproximado
% c_indice: Velocidade do som das regioes do objeto aproximado

%% Forma

% Deforma e desloca as mascaras da coluna, cora��o e pulmao
mask.colunaMask =  deforma(mask.colunaMask,0.5,iterationForma,deslocamento);
mask.coracaoMask =  deforma(mask.coracaoMask,0.5,iterationForma,deslocamento);

% separa o pulmao em duas partes e entao deforma cada pulm�o separadamene,
% depois junta cada um dos pulm�es
p = regionprops(mask.pulmoesMask,'PixelIdxList');
p1 = boolean(zeros(size(mask.pulmoesMask)));
p2 = p1;
p1(p(1).PixelIdxList) = 1;
p2(p(2).PixelIdxList) = 1;
p1 = deforma(p1,0.5,iterationForma,deslocamento);
p2 = deforma(p2,0.5,iterationForma,deslocamento);

mask.pulmoesMask = p1 | p2;

%% Intensidade

% Altera a velocidade de som no meio atr�ves de uma fun��o gaussiana
% c_meio permanece igual
c_indice.c_coracao = c_indice.c_coracao*(1+0.1*cIntensity*randn); 
c_indice.c_osso = c_indice.c_osso*(1+2*cIntensity*randn);
c_indice.c_pulmao = c_indice.c_pulmao*(1+cIntensity*randn);

%% Cria��o do medium

% Altera a velocidade do som do  medium de entrada de acordo com as novas
% mask e os novos indices
medium.sound_speed(mask.torsoMask) =  c_indice.c_meio;
medium.sound_speed(mask.colunaMask) =  c_indice.c_osso;
medium.sound_speed(mask.pulmoesMask) =  c_indice.c_pulmao;
medium.sound_speed(mask.coracaoMask) =  c_indice.c_coracao;

end