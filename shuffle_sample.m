function[X_shuffle]=shuffle_sample(X_class)

[~,a,~]=size(X_class);

random_vector=randperm(a);

for ii=1:a

    X_shuffle(:,ii)=X_class(:,random_vector(1,ii)); 


end


end