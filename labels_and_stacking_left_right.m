function[tensor_all_coefs_a_3D,labels]=labels_and_stacking_left_right(coef,tensor_coefs_a_4D)


    ix=0;
    [~,b,c,d]=size(tensor_coefs_a_4D);
       
    %Os sujeitos c,d,e são dados articifiais

    for isuj = 1:d
               
        %Empilhamento dos coeficientes do vetor a        
        t=1;
        
        for icanal=1:c            
            
        tensor_all_coefs_a_3D(t:t+coef-1,:,isuj)=tensor_coefs_a_4D(:,:,icanal,isuj);
        t=t+coef;
        
        end    
        
        %Carregando os labels        
        lab_left=ones(1,int32(b/2));
        lab_right=2*ones(1,int32(b/2));
        
        labels(:,:,isuj)=[lab_left lab_right];
        
        
       

    end
    
 
end