module InterPolynomialModule
#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#
export ilorazyRoznicowe
export warNewton
export naturalna
export rysujNnfx

using PyPlot

#=
    Implementation of function that finds   divided diffrences that are needed
    to find interpoling polynomials using Newtown form.
    Input:
        x - interpolation points
        f - values of the function in interpolation points.
    Variables:
        n - length of vectors from input
        fx - vector where already calculated divided diffrences are stored,
            since we calculate them recursively we use this vector not only to store them.
    Output
        fx - vector where already calculated divided diffrences are stored
=#

function ilorazyRoznicowe(x::Vector{Float64},f::Vector{Float64})
    n=length(x);
    fx=Vector{Float64}(n);
    for i in 1:n
        fx[i]=f[i];
    end
    for j in 2:n
        for i in n:-1:j
            fx[i]=(fx[i]-fx[i-1])/(x[i]-x[i-j+1]);
        end
    end
    return fx;
end

#=
    Implementation of function that computes value in given 't' of the interpoling polynomials
    in Newton form using Horner method.
    Input:
        x - interpolation points
        fx - values of the divided diffrences of given polynomial.
            fx[i]=f[x0x1...x(i-1)]
        t - point in which we compute value of polynomial
    Variables:
        nt - acumulator, we add to it new values and at the end we return nt.
        product - product is result of multiplication following parts of polynomial
            in Newton's form
    Output
        nt -  computed value
=#

function warNewton(x::Vector{Float64},fx::Vector{Float64},t::Float64)
    nt=Float64(fx[1]);
    product=Float64(1.0);
    for i in 2:length(x)
        product*=(t-x[i-1])
        nt+=fx[i]*product;
    end
    return nt;
end

#=
    Implementation of function that computes values of coefficients of
    interpoling polynomial in Newton's form built on interpoling points and
    divided diffrences from input.
    Input:
        x - interpolation points
        fx - values of the divided diffrences of given polynomial.
            fx[i]=f[x0x1...x(i-1)]
    Variables:
        n - length of vectors from input
        b - vector where we store coefficients
        from last iteration of this algorithm. After 'z' iterations
        first z elements in vector b are coefficients that we are looking for.
    Output
        b - vector with coefficients
=#

function naturalna(x::Vector{Float64},fx::Vector{Float64})
    n=length(fx);
    b=Vector{Float64}(n);
    for i in 1:n
        b[i]=fx[i]
    end
    for i in 1:n
        for j in n-1:-1:i
            b[j]=b[j]-x[j-i+1]*b[j+1];
        end
    end
    return b;
end

#=
    Implementation of function that returns graphic representation of
    interpolating polynomial of function 'f' between' a' and 'b'
    and function itself.
    Where degree of computing polynomial equals n.
    Input:
        f - function given as anynomous function
        a,b - interval which we explore
        n - degree of the polynomial
        graph_n - ammount of points that will be used to create
            polynomial's plot.
        graph_f_x,graph_w_x - values of function and polynomial
            used in graph representation of both.
    Variables:
        x - interpolation points
        y - y[i]=f(x[i])
        fx - values of the divided diffrences of given polynomial.
            fx[i]=f[x0x1...x(i-1)]

    Output
        Graphic representation of interpolating polynomial
=#

function rysujNnfx(f,a::Float64,b::Float64,n::Int)
    x=Vector{Float64}(n+1)
    y=Vector{Float64}(n+1)
    h=Float64((b-a)/n)
    for i in 0:n
        x[i+1]=a+i*h;
        y[i+1]=f(x[i+1]);
    end
    fx=ilorazyRoznicowe(x,y)
    #graph_n=Int(floor(b-a)*20)
    graph_n=n*100;
    graph_f_x=Vector{Float64}(graph_n+1)
    graph_w_x=Vector{Float64}(graph_n+1)
    graph_h=Float64((b-a)/graph_n)
    for i in 0:graph_n
        arg=a+i*graph_h;
        graph_f_x[i+1]=f(arg);
        graph_w_x[i+1]=warNewton(x,fx,arg);
    end

    plot(linspace(a,b,graph_n+1),graph_f_x,color="black",linewidth=1.0)
    plot(linspace(a,b,graph_n+1),graph_w_x,color="red",linewidth=1.0,)
end

end #end module
