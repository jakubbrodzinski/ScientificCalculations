module BlockSys
#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#
using FileHandling

export gaussDom;
export gauss;

export LUmatrixes;
export LUmatrixesDom;

export findXWithLU;

export computeError;

function computeError(x::Vector{Float64})
    exact=ones(Float64,length(x));
    return norm(exact - x)/norm(exact);
end
#=
    Function solves the equasion LUx=b, using the fact that we know L and U representation.
    Firstly it solves Lz=b and then Ux=z.
    Input:
        L,U - SparseMatrixCSC that L*U=A
        n - L's and U's sizes
        b - right side's vector
        l - size of the small L,U's small blocks
    Variables:
        offset - variable with which we control that loops won't look into empty cells.
    Output:
        x - vector where solution is stored
=#
function findXWithLU(L::SparseMatrixCSC{Float64, Int64},U::SparseMatrixCSC{Float64, Int64},b::Vector{Float64},n::Int64,l::Int64)
    offset=l-1;
    for i in 1:(l-2)
        for j in 1:offset
            b[i+j]=b[i+j]-L[i+j,i]*b[i];
        end
        offset=offset-1;
    end
    offset=l+1
    for i in (l-1):(n-2)
        for j in 1:offset
            b[i+j]=b[i+j]-L[i+j,i]*b[i];
        end
        offset=offset-1;
        if offset<2
            offset=l+1;
        end
    end
    b[n]=b[n]-b[n-1]*L[n,n-1];

    x=zeros(Float64,n);
    for i in n:-1:1
        sum=0.0;
        j=1
        while j<2*l && (i+j)<=n
            sum=sum+U[i,i+j]*x[i+j];
            j=j+1;
        end
        x[i]=(b[i]-sum)/U[i,i]
    end
    return x;
end
#=
    Function that using Gauss algoritm eleminates k rows that are below 'row'-th row.
    In every row functions changes up to M indexes.
    Input:
        A -SparseMatrixCSC
        b - right side vector
        row - A[row,row] is the cell from where functions starts.
        m - how much columns we taking care of.
        k - how much rows are we taking care of.
=#
function gaussKdownMright(A::SparseMatrixCSC{Float64, Int64},b::Vector{Float64},row::Int64,k::Int64,m::Int64)
    for i in (row+1):(row+k)
        factor=Float64(A[i,row]/A[row,row]);
        for g in row:(row+m-1)
            A[i,g]=A[i,g]-factor*A[row,g];
        end
        b[i]=b[i]-factor*b[row];
    end
end
#=
    After applying gauss elimination onto every row this function computes x_i
    (assuming that x_j for every j that j>i is already computed) from equasion Ax=b.
    Input:
        A -SparseMatrixCSC
        b- right side's vector
        row - which row we take care of
        depth -  how much cells in given row we should check for actual values.
=#
function findXi(A::SparseMatrixCSC{Float64, Int64},b::Vector{Float64},row::Int64,depth::Int64)
    for i in (row+1):(row+depth-1)
        b[row]=b[row]-A[row,i]*b[i];
        A[row,i]=0;
    end
    b[row]=b[row]/A[row,row];
    A[row,row]=1;
end
#=
    Function solves equasion Ax=b using Guass elimination algorithm and then solving
    equasion A'x=b', where A' is upper square matrix. Its uses
    'gaussKdownMright' and 'findXi'.
    Input:
        A  - SparseMatrixCSC from left side.
        b - right's side vector.
        n -size of A
        l - size of smaller blocks of A.
    Output:
        Solution of equasion
=#
function gauss(A::SparseMatrixCSC{Float64, Int64},b::Vector{Float64},n::Int64,l::Int64)
    x=zeros(Float64,n);
    temp=l;
    for i in 1:(n-l)
        if(A[i,i]<=eps(Float64))
            error("Use second method! Element close to zero.");
        end
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
#=
    Function solves equasion Ax=b using Guass elimination algorithm and then solving
    equasion A'x=b', where A' is upper square matrix. During Gauss elimination algorithm
    picks dominant element for every row that's why vectors p,s are needed. Its uses
    'gaussDomKdownMright' and 'findXiDom'.
    Input:
        A  - SparseMatrixCSC from left side.
        b - right's side vector.
        n -size of A
        l - size of smaller blocks of A.
    Variables:
        p - permutation vector
        s - vector that is used for picking the most dominant element.
    Output:
            Solution of equasion that was converted with vector p.
=#
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

#=
    Function that using Gauss algoritm eleminates k rows that are below 'row'-th row.
    In every row functions changes up to M indexes. This functions takes care of
    the most dominant matrix A's element aswell.
    Input:
        A -SparseMatrixCSC
        b - right side vector
        row - A[row,row] is the cell from where functions starts.
        m - how much columns we taking care of.
        k - how much rows are we taking care of.
=#
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

#=
    After applying gauss elimination onto every row this function computes x_i
    (assuming that x_j for every j that j>i is already computed) from equasion Ax=b.
    This function is used when picking dominant element is taken care of. Function needs
    additional parameter 'p'.
    Input:
        A -SparseMatrixCSC
        b- right side's vector
        row - which row we take care of
        depth -  how much cells in given row we should check for actual values.
        p - permutation vector.
=#
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
#=
    Functiona applies Gauss elimination algorithm on matrix U (from LU) and returns
    rows, columns and values of matrix L's cells that were created during that.
    Input:
        A - SparseMatrixCSC where we want to create U matrix
        row - A[row,row] is the cell from where functions starts.
        m - how much columns we taking care of.
        k - how much rows are we taking care of.
    Output :
        (r,c,v) - corridantes  and values of the cells of L's matrix that should be added
            to  L matrix.
=#
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
#=
    Function generates L and U matrxes from matrix A using method 'gaussKdownMrightLU'
    that returns rows,column and values of L matrix and changing matrix from input
    into U matrix at the same time.
    Input:
        A - matrix that later L*U=A
        n,l -size of A
    Output:
        (L,U) - computed matrixes.
=#
function LUmatrixes(A::SparseMatrixCSC{Float64, Int64},n::Int64,l::Int64)
    row=Int64[];
    column=Int64[];
    value=Float64[];
    temp=l;
    for i in 1:(n-l)
        if(A[i,i]<=eps(Float64))
            error("Use second method! Element close to zero.");
        end
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
#=
    Function generates L and U matrxes from matrix A using method 'gaussKdownMrightLU'
    that returns rows,column and values of L matrix and changing matrix from input
    into U matrix at the same time.
    Input:
        A - matrix that later L*U=A
        n,l -size of A
    Output:
        (L,U) - computed matrixes.
=#
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
    #println(p);
    U=convertSparseMatrixWithPermutation(A,p,n);
    p2=zeros(n);
    for i in 1:n
        p2[p[i]]=i;
    end
    sparse_row_p=Int64[];
    for i in 1:length(sparse_row)
        append!(sparse_row_p,p2[sparse_row[i]]);
    end
    L=sparse(sparse_row_p,sparse_column,sparse_value);
    properb=convertBWithPermutation(vector_b,p,n);
    return (L,U,properb);
end

#=
    Functiona applies Gauss elimination algorithm on matrix U (from LU) and returns
    rows, columns and values of matrix L's cells that were created during that.
    This function search for dominant elment before applying Gauss elimination algorithm,
     and because of that it needs p and s vector .
    Input:
        A - SparseMatrixCSC where we want to create U matrix
        row - A[row,row] is the cell from where functions starts.
        m - how much columns we taking care of.
        k - how much rows are we taking care of.
        p - permuation vector
        s - vector with values max(abs(A[i,j]))
    Output :
        (r,c,v) - corridantes  and values of the cells of L's matrix that should be added
            to  L matrix.
=#
function gaussDomKdownMrightLU(A::SparseMatrixCSC{Float64, Int64},row::Int64,k::Int64,m::Int64,p::Vector{Int64},s::Vector{Float64})
    s_r=Int64[];
    s_c=Int64[];
    s_v=Float64[];

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

    push!(s_r,p[row]);
    push!(s_c,row);
    push!(s_v,1.0);

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
        push!(s_r,p[i]);
        push!(s_c,row);
        push!(s_v,factor);
    end
    return (s_r,s_c,s_v);
end

end #module
