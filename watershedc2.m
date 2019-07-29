%   
%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%

function [cl]=watershedc(C, h)

N = size(C,1);
cl = zeros(N,1);

[r,c]=find(C);
v = C(find(C(:)));


if h>0
    maxv = max(v);
    mask = (v<(maxv-h));
    v = v.*mask+~mask*(maxv-h);
    
    ind = find(~mask);
    Cm = sparse(r(ind),c(ind),v(ind), size(C,1), size(C,1));
    
    cl = ff(Cm);
    
        
end


ind = r>c;
r = r(ind);
c = c(ind);
v = v(ind);



[v, ind] = sort(v, 'descend');
c = c(ind);
r = r(ind);



next_cluster_number=max(cl)+1;

for i=1:length(c)
    
    c1 = c(i);
    c2 = r(i);
    
    if (cl(c1)==0)&(cl(c2)==0)
        cl(c1) = next_cluster_number;
        cl(c2) = next_cluster_number;
        next_cluster_number = next_cluster_number+1;
    end
    
    if (cl(c1)==0)&(cl(c2)~=0)
        cl(c1) = cl(c2);
    end
    
    if (cl(c1)~=0)&(cl(c2)==0)
        cl(c2) = cl(c1);
    end
    
end









return
end

%
%   function cl = ff(C)
%
%	Forest-fire algorithm with the matrix of links
% 
%	v is the vector of values in matrix L
%	c is the vector of cluster numbers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cl = ff(C)

    C0 = (C>0);
    N = size(C,1);
	cl=zeros(N,1);                  % List of clusters
	cc=1;                           % Current cluster number
	n=0;
    
	for cn=1:N			%
            
        fraction_done = cn/N;
        %disp(char(repmat(8,1,50)))
            
            
        if cl(cn)==0				% Not classified yet
    							% Start new cluster

                v = sparse(cn, 1, 1, N, 1);
                
                
                fv = cn;
                dn=1;
                
                nstep=0;
                
                while(dn>0)
                    nstep=nstep+1;
                    nv = C0*v+v;
                    fnv = find(nv);
                    dn = length(fnv)-length(fv);
                    v=nv;
                    fv = fnv;
                end
                
                if nstep>1 % this is to exclude cluster with only one
                % element

                    ind = find(v);
                
                    cl(ind)=cc;
                
                    cc = cc+1;
                end
                fraction_done = cn/N;
        end
    end
    return
end



