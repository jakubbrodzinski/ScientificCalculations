function f(x::Float64)
    return e^x*log(e,1+e^(-x))
end

#=
    Commands that were used to generate graph:
    n=linspace(-20,40,100000);
    e=exp(1)
    x=(e.^n).*log(1+e.^(-n))
    %x=((exp(1)).^n).*log(1+(exp(1)).^(-n));
    plot(n,x)
=#
