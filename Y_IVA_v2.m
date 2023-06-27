function [Y_1,Y_2]=Y_IVA_v2(X_1,X_2,W_1,W_2,D,n)
    
           [~,b,D]=size(X_1);
           [~,e,~]=size(X_2);
           
%          Y_abfg_test=zeros(n*2,round(n_samples/2*(1-fator_treino)),D);
%          Y_left_test=zeros(n*2,round(n_samples/2*(1-fator_treino)),D);

           Y_1=zeros(n*2,b,D);
           Y_2=zeros(n*2,e,D);

      for it=1:D

          Y_1(1:n,:,it)=W_1(:,:,it)*X_1(:,:,it);
          Y_1(n+1:n*2,:,it)=W_2(:,:,it)*X_1(:,:,it);

          Y_2(1:n,:,it)=W_1(:,:,it)*X_2(:,:,it);
          Y_2(n+1:n*2,:,it)=W_2(:,:,it)*X_2(:,:,it);

      end


end