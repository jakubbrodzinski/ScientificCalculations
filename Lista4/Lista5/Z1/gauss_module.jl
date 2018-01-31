module GaussModule
#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

using FileHandling

export findXiDom;
export gaussDomKdownMright;
export gaussDom;
export gauss;
export findXi;
export gaussKdownMright;

function gaussKdownMright(A::SparseMatrixCSC{Float64, Int64},b::Vector{Float64},row::Int64,k::Int64,m::Int64)
    for i in (row+1):(row+k)
        factor=Float64(A[i,row]/A[row,row]);
        for g in row:(row+m-1)
            A[i,g]=A[i,g]-factor*A[row,g];
        end
        b[i]=b[i]-factor*b[row];
    end
end

function findXi(A::SparseMatrixCSC{Float64, Int64},b::Vector{Float64},row::Int64,depth::Int64)
    for i in (row+1):(row+depth-1)
        b[row]=b[row]-A[row,i]*b[i];
        A[row,i]=0;
    end
    b[row]=b[row]/A[row,row];
    A[row,row]=1;
end

function gauss(A::SparseMatrixCSC{Float64, Int64},b::Vector{Float64},n::Int64,l::Int64)
    temp=l;
    for i in 1:(n-l)
        if temp>2
            gaussKdownMright(A,b,i,temp-1,l+1);
        else
            gaussKdownMright(A,b,i,temp+l-1,l+1);
        end

        temp=temp-1;
        if temp==0
            temp=l;
        end;
    end
    temp=l;
    for i in (n-l+1):n
        gaussKdownMright(A,b,i,temp-1,temp);
        temp=temp-1;
    end
    temp=1;
    for i in n:-1:(n-l)
        findXi(A,b,i,temp)
        temp=temp+1;
    end
    for i in (n-l-1):-1:1
        findXi(A,b,i,l+1)
    end
    return b;
end

function gaussDom(A::SparseMatrixCSC{Float64, Int64},b::Vector{Float64},n::Int64,l::Int64)
    p=Vector{Int64}();
    s=zeros(n);
    rows = rowvals(A)
    vals = nonzeros(A)
    for i in 1:n
        push!(p,i);
        for j in nzrange(A, i)
            row = rows[j]
            val = abs(vals[j])
            if s[row]<val
                s[row]=val;
            end
        end
    end

    temp=l;
    for i in 1:(n-l)
        if temp>2
            gaussDomKdownMright(A,b,i,temp-1,2*l,p,s);
        else
            gaussDomKdownMright(A,b,i,temp+l-1,2*l,p,s);
        end

        temp=temp-1;
        if temp==0
            temp=l;
        end;
    end
    temp=l;
    for i in (n-l+1):n
        gaussDomKdownMright(A,b,i,temp-1,temp,p,s);
        temp=temp-1;
    end
    temp=1;
    for i in n:-1:(n-l)
        findXiDom(A,b,i,temp,p,s)
        temp=temp+1;
    end
    for i in (n-l-1):-1:1
        findXiDom(A,b,i,2*l,p,s)
    end
    return convertBWithPermutation(b,p,n);
end

function gaussDomKdownMright(A::SparseMatrixCSC{Float64, Int64},b::Vector{Float64},row::Int64,k::Int64,m::Int64,p::Vector{Int64},s::Vector{Float64})
    r=p[row];
    max=abs(A[r,row])/s[r];
    swap=row;
    for i in (row+1):(row+k)
        r=p[i];
        max_local=abs(A[r,row])/s[r];
        if max_local>max
            max=max_local;
            swap=i;
        end
    end
    p[row],p[swap]=p[swap],p[row];
    s[row],s[swap]=s[swap],s[row];


    for i in (row+1):(row+k)
        factor=Float64(A[p[i],row]/A[p[row],row]);
        max=row+m-1;
        if max>A.n
            max=A.n
        end
        for g in row:max
            A[p[i],g]=A[p[i],g]-factor*A[p[row],g];
        end
        b[p[i]]=b[p[i]]-factor*b[p[row]];
    end
end

function findXiDom(A::SparseMatrixCSC{Float64, Int64},b::Vector{Float64},row::Int64,depth::Int64,p::Vector{Int64},s::Vector{Float64})
    max=row+depth-1;
    if max>A.n
        max=A.n;
    end
    for i in (row+1):max
        b[p[row]]=b[p[row]]-A[p[row],i]*b[p[i]];
        A[p[row],i]=0;
    end
    b[p[row]]=b[p[row]]/A[p[row],row];
    A[p[row],row]=1;
end

end #module
