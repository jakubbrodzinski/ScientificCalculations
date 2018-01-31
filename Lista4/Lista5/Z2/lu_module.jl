module LUGauss
#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

using FileHandling

export gaussKdownMrightLU;
export LUmatrixes;
export LUmatrixesDom;
export gaussDomKdownMrightLU;

function gaussKdownMrightLU(A::SparseMatrixCSC{Float64, Int64},row::Int64,k::Int64,m::Int64)
    r=Int64[];
    c=Int64[];
    v=Float64[];
    push!(r,row);
    push!(c,row);
    push!(v,1.0);
    for i in (row+1):(row+k)
        factor=Float64(A[i,row]/A[row,row]);
        A[i,row]=0.0;
        for g in (row+1):(row+m-1)
            A[i,g]=A[i,g]-factor*A[row,g];
        end
        push!(r,i);
        push!(c,row);
        push!(v,factor);
    end
    return (r,c,v);
end

function LUmatrixes(A::SparseMatrixCSC{Float64, Int64},n::Int64,l::Int64)
    row=Int64[];
    column=Int64[];
    value=Float64[];
    temp=l;
    for i in 1:(n-l)
        if temp>2
            (a,b,c)=gaussKdownMrightLU(A,i,temp-1,l+1);
            append!(row,a);
            append!(column,b);
            append!(value,c);
        else
            (a,b,c)=gaussKdownMrightLU(A,i,temp+l-1,l+1);
            append!(row,a);
            append!(column,b);
            append!(value,c);
        end

        temp=temp-1;
        if temp==0
            temp=l;
        end;
    end
    temp=l;
    for i in (n-l+1):n
        (a,b,c)=gaussKdownMrightLU(A,i,temp-1,temp);
        append!(row,a);
        append!(column,b);
        append!(value,c);

        temp=temp-1;
    end
    return (sparse(row,column,value),A);
end

function LUmatrixesDom(A::SparseMatrixCSC{Float64, Int64},vector_b::Vector{Float64},n::Int64,l::Int64)
    sparse_row=Int64[];
    sparse_column=Int64[];
    sparse_value=Float64[];

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
            (a,b,c)=gaussDomKdownMrightLU(A,i,temp-1,2*l,p,s);
            append!(sparse_row,a);
            append!(sparse_column,b);
            append!(sparse_value,c);
        else
            (a,b,c)=gaussDomKdownMrightLU(A,i,temp+l-1,2*l,p,s);
            append!(sparse_row,a);
            append!(sparse_column,b);
            append!(sparse_value,c);
        end

        temp=temp-1;
        if temp==0
            temp=l;
        end;
    end
    temp=l;
    for i in (n-l+1):n
        (a,b,c)=gaussDomKdownMrightLU(A,i,temp-1,temp,p,s);
        append!(sparse_row,a);
        append!(sparse_column,b);
        append!(sparse_value,c);

        temp=temp-1;
    end
    U=convertSparseMatrixWithPermutation(A,p,n);
    L=sparse(sparse_row,sparse_column,sparse_value)
    properb=convertBWithPermutation(vector_b,p,n);
    return (L,U,properb);
end

function gaussDomKdownMrightLU(A::SparseMatrixCSC{Float64, Int64},row::Int64,k::Int64,m::Int64,p::Vector{Int64},s::Vector{Float64})
    s_r=Int64[];
    s_c=Int64[];
    s_v=Float64[];
    push!(s_r,row);
    push!(s_c,row);
    push!(s_v,1.0);
#computing dominant element
    max=abs(A[p[row],row])/s[p[row]];
    swap=row;
    for i in (row+1):(row+k)
        max_local=abs(A[p[i],row])/s[p[i]];
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
        A[p[i],row]=0.0;
        for g in (row+1):max
            A[p[i],g]=A[p[i],g]-factor*A[p[row],g];
        end
        push!(s_r,i);
        push!(s_c,row);
        push!(s_v,factor);
    end
    return (s_r,s_c,s_v);
end

end #module
