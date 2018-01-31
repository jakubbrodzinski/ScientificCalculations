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

function loadA(filename::String)
    open(filename) do file
        (n,l)=split(readline(file));
        (n,l)=(parse(Int64,n),parse(Int64,l));
        println("n: ",n);
        println("l: ",l);
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

function saveSolution(b::Vector{Float64},filename::String)
    open(filename,"w") do file
        for i in 1:size(b,1)
            println(file,b[i])
        end
    end
end

function saveSolutionWithError(b::Vector{Float64},err::Float64,filename::String)
    open(filename,"w") do file
        println(err);
        for i in 1:size(b,1)
            println(file,b[i])
        end
    end
end

function compute_b(A::SparseMatrixCSC{Float64, Int64},n::Int64,l::Int64)
    b=Vector{Float64}();
    temp=1;
    for i in 1:l
        sum=0;
        for j in 1:l
            sum=sum+A[i,j];
        end
        sum=sum+A[i,l+temp];
        temp=temp+1;
        push!(b,sum);
    end
    temp=1;
    offset=l-1;
    for i in (l+1):(n-l)
        sum=0;
        for j in offset:(offset+l+1)
            sum=sum+A[i,j];
        end
        sum=sum+A[i,offset+l+1+temp];
        temp=temp+1;
        if(temp>l)
            temp=1
        end
        push!(b,sum);

        if i%l==0
            offset=offset+l;
        end
    end
    for i in (n-l+1):n
        sum=0;
        for j in (n-l-1):n
            sum=sum+A[i,j];
        end
        push!(b,sum);
    end
    return b;
end

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
