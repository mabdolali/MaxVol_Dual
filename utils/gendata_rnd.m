function [ M, W, H ] = gendata_rnd(m,r,purity,Nif,Nii)
% m : dimension
% r : #endmembers
% Nif : # of points on each facet
% Nii : # of points inside facet

W = rand(m,r);
% generating points on facets
H1=zeros(r,sum(Nif));
facets=nchoosek(1:r,r-1);
concentration = 1;
if purity < 0.4
    concentration = 1000;
end
for i=1:length(facets)
    starti=sum(Nif(1:i-1))+1;
    endi=sum(Nif(1:i));
    H1(facets(i,:),starti:endi)=sample_dirichlet((concentration/(r-1))*ones(1,r-1),Nif(i))';
    for j = starti : endi
        while max( H1(facets(i,:),j) ) > purity
            H1(facets(i,:),j) = sample_dirichlet((concentration/(r-1))*ones(1,r-1),1)';
        end
    end
end

H1=H1./repmat(sum(H1,1),[r 1]);
M1=W*H1;

% generating points inside the facet
H2 = sample_dirichlet((concentration/r)*ones(1,r),Nii)';
for j = 1 : Nii
        while max( H2(:,j) ) > (purity) 
            H2(:,j) = sample_dirichlet((concentration/r)*ones(1,r),1)';
        end
end
H2=H2./repmat(sum(H2,1),[r 1]);
M2=W*H2;
H= [H1 H2];
M=[M1 M2];
end

