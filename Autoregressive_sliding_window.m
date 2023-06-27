function[tensor_coefs_a_4D_test]=Autoregressive_sliding_window(sliding_window,Y_test,D,coef)

    [r,n_samples,D]=size(Y_test);
    
    for i_dataset=1:D

               
        for ic=1:r
            
           t=0;      
           for is=1:sliding_window:(n_samples-(sliding_window-1))

               y=Y_test(ic,is:is+sliding_window-1,i_dataset);   
               
               if y(1,end)==0
                break
               end 
               
               y=normalize(y);

               for id=1:coef

                   Y_desloc(id,:)=y(1,id:(sliding_window-coef+id)); 

               end

               [a,b]=size(Y_desloc);
               Ry_desloc=(Y_desloc*Y_desloc')/b;

               ry_desloc=Y_desloc(:,1:b-1)*(y(1,(coef+1):end))'./(b-1);

               t=t+1;
               iii=inv(Ry_desloc);
               if iii(1,1)~=Inf
                Matriz_a(:,t)=inv(Ry_desloc)*ry_desloc;              
               end           

           end 
 
           [g,h,m]=size(Matriz_a);
           Tensor_a(:,1:h,ic)=Matriz_a;

        end
        
        tensor_coefs_a_4D_test(:,1:h,:,i_dataset)= Tensor_a;

    end
   
end