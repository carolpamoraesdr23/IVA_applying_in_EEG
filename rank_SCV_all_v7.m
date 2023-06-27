function[SCVs1,rank_SCVs_ps1]=rank_SCV_all_v7(Sigma_N1)  
    [a,~,c]=size(Sigma_N1);
    vector1=zeros(10000,5);
    t1=1;   
       
    
    for kk=1:c
            
            for ii=1:a
                for jj=1:a
            
                    if (abs(Sigma_N1(ii,jj,kk))>0.01) && (ii~=jj)
                        vector1(t1,1)=ii;
                        vector1(t1,2)=jj;
                        vector1(t1,3)=abs(Sigma_N1(ii,jj,kk));
                        vector1(t1,5)=1;
                      
                        
                        if abs(Sigma_N1(ii,jj,kk))>0.1
                            vector1(t1,4)=1;
                        else 
                            vector1(t1,4)=0;
                        end
        
                        t1=t1+1;
                    end
            
                end
            end
     end
    
    

    
    SCVs1=vector1(1:t1-1,:);
    [value1,position1]=sort(SCVs1(:,3),'descend');
    
    for pp=1:t1-1
        rank_SCVs1(pp,:)=SCVs1(position1(pp,1),:);
    end
    
    
    for pp=1:a
        vv=1;
        for rr=1:t1-1
            if rank_SCVs1(rr,1)==pp
                rank_SCVs_ps1(vv,1,pp)=rank_SCVs1(rr,2);
                rank_SCVs_ps1(vv,2,pp)=rank_SCVs1(rr,3);
                vv=vv+1;
          
            end

        end
    end




end
