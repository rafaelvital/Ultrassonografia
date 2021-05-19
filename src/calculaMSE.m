function [MSE] = calculaMSE(ensaio,resposta)
% Dado o ensaio realizada a imagem tomografica(resposta) calcula o MSE de
% cada uma das partes: Pulmão, coração, coluna e meio.
% MSE é um vertor coluna onde
% MSE(1) = MSE Pulmão
% MSE(2) = MSE Coração
% MSE(3) = MSE Coluna
% MSE(4) = MSE Meio
% MSE(5) = NMSE Pulmão
% MSE(6) = NMSE Coração
% MSE(7) = NMSE Coluna
% MSE(8) = NMSE Meio


%  Pulmao
c_objeto = resposta(resize(ensaio.mask.pulmoesMask,size(resposta),'nearest'));
MSE(1,1) = sqrt(sum((c_objeto - ensaio.c_indice.c_pulmao).* (c_objeto - ensaio.c_indice.c_pulmao)) ./ length(c_objeto));
MSE(1,1+4) = 100*sqrt(sum((c_objeto - ensaio.c_indice.c_pulmao).* (c_objeto - ensaio.c_indice.c_pulmao)) ./ (length(c_objeto).* ensaio.c_indice.c_pulmao .* ensaio.c_indice.c_pulmao));

%  Coração
c_objeto = resposta(resize(ensaio.mask.coracaoMask,size(resposta),'nearest'));
MSE(1,2) = sqrt(sum((c_objeto - ensaio.c_indice.c_coracao).* (c_objeto - ensaio.c_indice.c_coracao)) ./ length(c_objeto));
MSE(1,2+4) = 100*sqrt(sum((c_objeto - ensaio.c_indice.c_coracao).* (c_objeto - ensaio.c_indice.c_coracao)) ./ (length(c_objeto).* ensaio.c_indice.c_coracao .* ensaio.c_indice.c_coracao));

%  Osso
c_objeto = resposta(resize(ensaio.mask.colunaMask,size(resposta),'nearest'));
MSE(1,3) = sqrt(sum((c_objeto - ensaio.c_indice.c_osso).* (c_objeto - ensaio.c_indice.c_osso)) ./ length(c_objeto));
MSE(1,3+4) = 100*sqrt(sum((c_objeto - ensaio.c_indice.c_osso).* (c_objeto - ensaio.c_indice.c_osso)) ./ (length(c_objeto).* ensaio.c_indice.c_osso .* ensaio.c_indice.c_osso));

%  Meio
mask = logical(ensaio.mask.torsoMask-ensaio.mask.pulmoesMask- ensaio.mask.coracaoMask-ensaio.mask.colunaMask);

c_objeto = resposta(resize(mask,size(resposta),'nearest'));
MSE(1,4) = sqrt(sum((c_objeto - ensaio.c_indice.c_meio).* (c_objeto - ensaio.c_indice.c_meio)) ./ length(c_objeto));
MSE(1,4+4) = 100*sqrt(sum((c_objeto - ensaio.c_indice.c_meio).* (c_objeto - ensaio.c_indice.c_meio)) ./ (length(c_objeto).* ensaio.c_indice.c_meio .* ensaio.c_indice.c_meio));

end