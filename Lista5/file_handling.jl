module FileHandling
#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#
export loadA;
export loadB;
export saveSolution;
export saveSolutionWithError;
export compute_b;
export convertSparseMatrixWithPermutation;
export convertBWithPermutation;

#=
    Function that creates SparseMatrixCSC from file named 'filename'.
    Input:
        filename- name of file from which Matrix,size n and size l will be read
    Variables:
        row,column,value - arrays that store rows,columns and values fo the non-zero
            fields in SparseMatrixCSC.
    Output:
        (a,b,c) - where
            a -SparseMatrixCSC from file,
            b - size of the returned matrix
            c - size of smaller blocks of matrix a
=#
function loadA(filename::String)
    open(filename) do file
        (n,l)=split(readline(file));
        (n,l)=(parse(Int64,n),parse(Int64,l));
        #println("n: ",n);
        #println("l: ",l);
        row=Int64[]
        column=Int64[];
        value=Float64[];
        while !eof(file)
            line=split(readline(file));
            push!(row,parse(Int64,line[1]));
            push!(column,parse(Int64,line[2]));
            push!(value,parse(Float64,line[3]));
        end
        return (sparse(row,column,value),n,l);
    end
end
#=
    Function that reads vector b from file named 'filename'.
    Input:
        filename- name of file from which Matrix,size n and size l will be read
=#
function loadB(filename::String)
    open(filename) do file
        n=parse(Int64,readline(file));
        b=Vector{Float64}();
        while !eof(file)
            push!(b,parse(Float64,readline(file)));
        end
        return (b,n);
    end
end
#=
    Function that saves vector b into file named 'filename'.
    Input:
        b - vector that procedure saves.
        filename- name of file from which Matrix,size n and size l will be read
=#
function saveSolution(b::Vector{Float64},filename::String)
    open(filename,"w") do file
        for i in 1:size(b,1)
            println(file,b[i])
        end
    end
end
#=
    Function that saves vector b and err (at the beginning) into file named 'filename'.
    Input:
        err - value that will be saved in 1 st line.
        b - vector that procedure saves.
        filename- name of file from which Matrix,size n and size l will be read
=#
function saveSolutionWithError(b::Vector{Float64},err::Float64,filename::String)
    open(filename,"w") do file
        println(file,err);
        for i in 1:size(b,1)
            println(file,b[i])
        end
    end
end
#=
    Function computes vector b, from SparseMatrixCSC A.
    Input:
        A- matrix
        n- size of matrix
    Output:
        b - vector that A*[1...,1]^{T}=b
=#
function compute_b(A::SparseMatrixCSC{Float64, Int64},n::Int64)
    b = zeros(Float64, n)
    vals = nonzeros(A)
    rows = rowvals(A)
    for i = 1:n
        for j in nzrange(A, i)
            b[rows[j]] = b[rows[j]] + vals[j]
        end
    end
    return b;
end
#=
    Function that creates new SparseMatrixCSC using old SparseMatrixCSC A and A's row's permutation vector.
    Input:
        A - matrix
        p - A's rows' permutation vector
    Output:
        Matrix with swapped rows after applying given permutation.
=#
function convertSparseMatrixWithPermutation(A::SparseMatrixCSC{Float64,Int64},p::Vector{Int64},n::Int64)
    n_row=Int64[];
    n_column=Int64[];
    n_value=Float64[];

    p2=zeros(n);
    for i in 1:n
        p2[p[i]]=i;
    end
    rows = rowvals(A)
    vals = nonzeros(A)
    for i in 1:n
        for j in nzrange(A, i)
              append!(n_row,p2[rows[j]]);
              append!(n_column,i);
              append!(n_value,vals[j]);
          end
    end
    return sparse(n_row,n_column,n_value);
end
#=
    Function permuates vector b using vector p.
=#
function convertBWithPermutation(b::Vector{Float64},p::Vector{Int64},n::Int64)
    p2=zeros(Int64,n);
    for i in 1:n
        p2[p[i]]=i;
    end
    retB=zeros(Float64,n);
    for i in 1:n
        retB[i]=b[p[i]];
    end
    return retB;
end

end #end module
