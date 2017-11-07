clear
close all
clc

pkg load statistics

function visualize(matrix, numVisualize, height, width, filename)
  matrix_visual = zeros(height*sqrt(numVisualize),width*sqrt(numVisualize));

  figure;
  ix = 1;
  ix1 = 1;
  for i = 1:sqrt(numVisualize)
      ix2 = 1;
      for j = 1:sqrt(numVisualize)
          curV = matrix(:,ix);
          curV = reshape(uint8(curV),[height, width]);
          matrix_visual(ix1:(ix1+height-1), ix2:(ix2+width-1)) = curV;
          ix2 = ix2 + width;
          ix = ix+1;
      end
      ix1 = ix1 + height;
  end

  imagesc(matrix_visual); axis ij
  caxis([0 255]);
  colormap(gray);
  imwrite(matrix_visual,colormap(gray),filename);
end

function [W, H] = em_algorithm_step (V, alpha_w, beta_w, alpha_h, beta_h, W_t, H_t, F, K, N)
  V_hat = W_t * H_t;
  W = ( W_t .* ( (V ./ V_hat) * H_t.' ) + (alpha_w - 1) * ones(F, K) ) ./ (ones(F, N) * H_t.' + beta_w *  ones(F, K));
  H = ( H_t .* ( W_t.' * (V ./ V_hat) ) + (alpha_h - 1) * ones(K, N) ) ./ (W_t.' * ones(F, N) + beta_h *  ones(K, N));
endfunction

height = 112;
width = 92;
numVisualize = 400; %must be a square
F = height * width
K = 25
N = 400

iterations = 15
alpha_h = 1
alpha_w = 1
beta_w = 1
beta_h = 1

filename = strcat("alph_", num2str(alpha_h), "_alpw_", num2str(alpha_w), 
"_beth_", num2str(beta_w), "_betw_", num2str(beta_h), ".png")

% This will load the matrix V
load('./attfaces.mat');

for beta_h = [0.1,0.2,1,5,10]
  for beta_w = [0.1,0.2,1,5,10]
    filename = strcat("alph_", num2str(alpha_h), "_alpw_", num2str(alpha_w),
    "_beth_", num2str(beta_w), "_betw_", num2str(beta_h), ".png")
    
    W_t = random('gam', alpha_w, beta_w, F, K);
    H_t = random('gam', alpha_h, beta_h, K, N);
    
    for c = 1:iterations  
      [W_t, H_t] = em_algorithm_step (V, alpha_w, beta_w, alpha_h, beta_h, W_t, H_t, F, K, N);
    end
    
    V_hat = W_t * H_t;
    
    visualize(V_hat, numVisualize, height, width, strcat("V_", filename))
    visualize(W_t, K, height, width, strcat("W_", filename))
  end
end



