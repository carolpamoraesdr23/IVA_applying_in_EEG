function[V_class1,V_class2,Z_train_class1,Z_train_class2,X_train_class1,X_train_class2,X_test_class1,X_test_class2]=upload_dataset_whitening(ii,ff,trials,window_samples,s,fator_treino,n,D)

    %%
    %Upload Dataset   
    
    

    for ns = 1:length(s)

        load(strcat('BCICIV_calib_ds1',num2str(s(ns)),'.mat'))  
        
              
        
        p=1;
        t=1;
        q=1;
        for it = 1:trials
            
            if (mrk.y(it)==1)
            X_abfg_original(:,q:q+window_samples-1,ns)=cnt(mrk.pos(it):mrk.pos(it)+window_samples-1,ii:ff)';
            q=q+window_samples;
            elseif (mrk.y(it)==-1)            
            X_left_original(:,t:t+window_samples-1,ns)=cnt(mrk.pos(it):mrk.pos(it)+window_samples-1,ii:ff)'; 
            t=t+window_samples;
            end           
%             X(:,p:p+window_samples-1,ns)=cnt(mrk.pos(it):mrk.pos(it)+window_samples-1,1:n)';
%             labels_imagery(1,it,ix)=mrk.pos(it);

            X_no_control_original(:,p:p+window_samples-1,ns)=cnt(mrk.pos(it)+window_samples:mrk.pos(it)+2*window_samples-1,ii:ff)';
            p=p+window_samples;


            
            
        end

    end

        
    %%
    %Whitening Dataset
    %Dados Imagery Training
    X_class1_original=double(X_abfg_original);
    X_class2_original=double(X_left_original);
    
    [X_class1]=shuffle_upload(X_class1_original,D,window_samples);
    [X_class2]=shuffle_upload(X_class2_original,D,window_samples);

   
    [~,n_samples,~]=size(X_class1);
    
    n_samples_split=(n_samples)*fator_treino;     
     
    X_train_class1=X_class1(:,1:n_samples_split,:); 
    X_test_class1=X_class1(:,n_samples_split+1:end,:); 
    X_train_class2=X_class2(:,1:n_samples_split,:); 
    X_test_class2=X_class2(:,n_samples_split+1:end,:); 


    X_test_class_total=zeros(n,int32(2*(n_samples-n_samples_split)),D);
    %Garantir que a media nula e Normalização
    for imt=1:D

        X_test_class_total(:,:,imt)=[X_test_class1(:,:,imt) X_test_class2(:,:,imt)];
        
    end
    
    M_train_class1=zeros(n,n_samples_split,D);
    M_test_class_total=zeros(n,int32(2*(n_samples-n_samples_split)),D);
    M_train_class2=zeros(n,n_samples_split,D);

        for im=1:D
            
            %Training Data class 1
            %Média nula
            
            M_train_class1(:,:,im)=mean(mean(X_train_class1(:,:,im)));
            X_train_class1(:,:,im)=X_train_class1(:,:,im)-M_train_class1(:,:,im);            
            %Normalização
            X_train_class1(:,:,im)=X_train_class1(:,:,im)./norm(X_train_class1(:,:,im));
            
            
            %Training Data class 2
            %Média nula
            
            M_train_class2(:,:,im)=mean(mean(X_train_class2(:,:,im)));
            X_train_class2(:,:,im)=X_train_class2(:,:,im)-M_train_class2(:,:,im);            
            %Normalização
            X_train_class2(:,:,im)=X_train_class2(:,:,im)./norm(X_train_class2(:,:,im));
            
            %Test Data 
            %Média nula
            
            M_test_class_total(:,:,im)=mean(mean(X_test_class_total(:,:,im)));
            X_test_class_total(:,:,im)=X_test_class_total(:,:,im)-M_test_class_total(:,:,im);            
            %Normalização
            X_test_class_total(:,:,im)=X_test_class_total(:,:,im)./norm(X_test_class_total(:,:,im));

        end
  
	[~,ms1,~]=size(X_test_class_total);

        for imt2=1:D                       
                      
            X_test_class1(:,:,imt2)=X_test_class_total(:,1:int32(ms1/2),imt2);
            X_test_class2(:,:,imt2)=X_test_class_total(:,int32(ms1/2)+1:ms1,imt2);
        
        end


    %Branqueamento do Sinal
    %Calculo autovalores, autovetores e v
   
        for iv=1:D

            %Whitening class 1             
            [E1, D1] = eig(cov(X_train_class1(:,:,iv)'));
            V_class1(:,:,iv) =E1*D1^(-0.5)*E1';  
            
            %Whitening class 2            
            [E2, D2] = eig(cov(X_train_class2(:,:,iv)'));
            V_class2(:,:,iv) =E2*D2^(-0.5)*E2'; 

        end



    %Sinal Branqueado
    
        for iz=1:D

            %Signal whitened class 1            
            Z_train_class1(:,:,iz)=V_class1(:,:,iz)*X_train_class1(:,:,iz);
            
            %Signal whitened class 2
            
            Z_train_class2(:,:,iz)=V_class2(:,:,iz)*X_train_class2(:,:,iz);
     

        end    
     
         
   
    
end
