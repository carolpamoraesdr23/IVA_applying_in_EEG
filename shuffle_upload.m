function[X_shuffle]=shuffle_upload(X_class,D,window_samples)

 [~,a,~]=size(X_class);
 b=int32(a/window_samples);
 
 vector_position=1:window_samples:a;
 
 for is=1:D
     
     t=1;
     random_vector=randperm(b);
     
     for ii=1:b
         
         X_shuffle(:,t:t+window_samples-1,is)=X_class(:,vector_position(1,random_vector(1,ii)):vector_position(1,random_vector(1,ii))+window_samples-1,is); 
         t=t+window_samples;

     end
 end

end