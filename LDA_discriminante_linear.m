function[labels_predict_test_total,Mdl,LDA_M_confu_test,LDA_M_confu_train,LDA_Acuracia_test,LDA_Acuracia_train]=LDA_discriminante_linear(s,input_X_train,input_y_labels_train,input_X_test,input_y_labels_test)
    
    [~,b,~]=size(input_X_train);    
    [~,e,~]=size(input_X_test);
    
    vector_position=1:200;

    
    %Os sujeitos c,d,e são dados articifiais
    ix=0;

    for ns = 1:length(s)
        
        ix=ix+1; 
        
        %Dados de Treinamento
        X_train=input_X_train(:,:,ix)';
        y_labels_train=input_y_labels_train(:,:,ix)';

        %Dados de Teste
        X_test=input_X_test(:,:,ix)';
        y_labels_test=input_y_labels_test(:,:,ix)';      
        
         
        Mdl = fitcdiscr(X_train,y_labels_train);
        
        labels_predict_test = predict(Mdl,X_test);

        labels_predict_test_total(:,:,ix)=labels_predict_test;
        
        %Matriz de Confusão
        LDA_M_confu_test(:,:,ix)=confusionmat(y_labels_test,labels_predict_test);
        
        
        %Acurácia de Teste
        acerto_test=0;
        cont_im_test=0;

       
        for it=1:e

            if y_labels_test(it,1)==labels_predict_test(it,1)

                acerto_test=acerto_test+1;         
        
            end


        end

           LDA_Acuracia_test(ix,:)=(acerto_test/e)*100;
                       
        %Acurácia de Treino
        
        labels_predict_train = predict(Mdl,X_train);
        
        %Matriz de Confusão
        LDA_M_confu_train(:,:,ix)=confusionmat(y_labels_train,labels_predict_train);
        
        acerto_train=0;
        cont_im_train=0;

        for it=1:b

            if y_labels_train(it,1)==labels_predict_train(it,1)

                acerto_train=acerto_train+1;

            end


        end

            
            LDA_Acuracia_train(ix,:)=(acerto_train/b)*100;
       

    end
    

end