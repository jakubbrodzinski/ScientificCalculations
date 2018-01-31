module RootsModule

#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#
export mbisekcji;
export mstycznych;
export msiecznych;

global BIS_OK=0
global BIS_SIGN_ERROR=1
global BIS_NEGATIVE_PRECISION=2;

#=
    Implementation of bisection method that is used to find function's root.
    Input:
        f - function that's roots we are looking for
        a,b - start and the end of the interval where the root is looked for
        delta, epsilon - precision of the calculations
    Variables
        l - in every iteration it's the start of the current interval
        r - in every iteration it's the end of the current interval
        m - in every iteration it's the middle of the current interval
        l_value - in every iteration it's the value of the first element in current interval
        r_value - in every iteration it's the value of the last element in current interval
        m_value - in every iteration it's the value in the middle in current interval
        it - iterator that counts number of iterations
    Output:
        (r,v,it,err) - tuple that contains:
            r - approximation of the root of the function f in [a,b]
            v - f(r)
            it - number of interation that were made
            err -   BIS_OK - 0 - no error occured
                    BIS_SIGN_ERROR - 1 - sign(f(a))==sign(f(b))
                    BIS_NEGATIVE_PRECISION - 2 - delta <= 0 or epsilon <=0
=#
function mbisekcji(f, a::Float64, b::Float64, delta::Float64, epsilon::Float64)
    l=a;
    r=b;
    l_value=f(l);
    r_value=f(r);

    if(sign(l_value)==sign(r_value))
        return (0,0,0,BIS_SIGN_ERROR);
    elseif(delta <=0 || epsilon<=0)
        return (0,0,0,BIS_NEGATIVE_PRECISION);
    end
    dif=r-l;
    dif=dif/2.0;
    m=l+dif
    m_value=f(m);
    it=1;
    while( abs(m)>delta && abs(m_value)>epsilon)
        if(sign(l_value)!=sign(m_value))    #if this is happening then function's root is in the left part of interval
            r=m;
            r_value=m_value
        else
            l=m;
            l_value=m_value;
        end
        dif=dif/2.0;
        m=l+dif;
        m_value=f(m);

        it=it+1;
    end
    return (m,m_value,it,BIS_OK);

end

global NEW_OK=0
global NEW_ITERATION_LIMIT=1
global NEW_DERIVATIVE_ALMOST_0=2
global NEW_NEGATIVE_PRECISION=3;

#=
    Implementation of Netwon's method that is used to find function's root.
    Input:
        f, pf- function that's roots we are looking for and function's derivative
        x0 - x from which we start iterations
        delta, epsilon - precision of the calculations
        maxit - maximum number of iterations
    Variables
        x_n -  at the start of every loop's iteration it's x_(n-1), then is swaped to x_(n), current approximation.
        x_n_value - in every iteration it's f(x_n)
        x_n_dvalue - in everyiteration it's df(x_n)
        it - iterator that counts number of iterations
    Output:
        (r,v,it,err) - tuple that contains:
            r - approximation of the root of the function f in [a,b]
            v - f(r)
            it - number of interation that were made
            err -   NEW_OK - 0 - no error occured
                    NEW_ITERATION_LIMIT - 1 - during maxint iterations required precision wasnt achieved
                    NEW_DERIVATIVE_ALMOST_0 - 2 - derivative is close to 0
                    NEW_NEGATIVE_PRECISION - 3 - delta <= 0 or epsilon <=0
=#
function mstycznych(f, pf,x0::Float64, delta::Float64, epsilon::Float64,maxit::Int)
    x_n = x0;
    x_n_value=f(x_n);

    if(delta <=0 || epsilon<=0)
        return (0,0,0,NEW_NEGATIVE_PRECISION);
    end
    it=0;

    for it in 1:maxit
        x_n_dvalue=pf(x_n);
        if(abs(x_n_dvalue)<eps(Float64))    #checking if derivative is close to 0.0
            return (0,0,it,NEW_DERIVATIVE_ALMOST_0);
        end
        x_p=x_n;
        x_n=Float64(x_n-(x_n_value/x_n_dvalue));
        x_n_value=f(x_n);

        if(abs(x_n-x_p)<delta || abs(x_n_value)< epsilon)             #if we have enough big precision we can end
            return (x_n,x_n_value,it,NEW_OK);
        end
    end
    return (x_n,x_n_value,it,NEW_ITERATION_LIMIT);


end

global SEC_OK=0
global SEC_ITERATION_LIMIT=1
global SEC_NEGATIVE_PRECISION=2;

#=
    Implementation of secant method that is used to find function's root.
    Input:
        f - function that's roots we are looking for
        x0,x1 - pair of x from which we start iterations
        delta, epsilon - precision of the calculations
        maxit - maximum number of iterations
    Variables
        x_n - at the start of every loop's iteration it's x_(n), then is swaped to x_(n+1) which is the closet approximation in given iteration.
        x_p - previous the most close approximation in given iteration it's x_(n-1) then it's swapped to x_(n)
        fx_n - in every iteration it's f(x_n)
        fx_p - in everyiteration it's f(x_p)
        it - iterator that counts number of iterations
    Output:
        (r,v,it,err) - tuple that contains:
            r - approximation of the root of the function f in [a,b]
            v - f(r)
            it - number of interation that were made
            err -   SEC_OK - 0 - no error occured
                    SEC_ITERATION_LIMIT - 1 - during maxint iterations required precision wasnt achieved
                    SEC_NEGATIVE_PRECISION - 2 - delta <= 0 or epsilon <=0
=#
function msiecznych(f,x0::Float64,x1::Float64, delta::Float64, epsilon::Float64,maxit::Int)
    x_n=x1;
    x_p=x0;
    fx_n=f(x_n);
    fx_p=f(x_p);

    if(delta <=0 || epsilon<=0) #input error
        return (0,0,0,SEC_NEGATIVE_PRECISION);
    end

    for it in 1:maxit
        if(abs(fx_n)>abs(fx_p))
            x_n, x_p=x_p, x_n;
            fx_n, fx_p= fx_p, fx_n;
        end
        temp=(x_n-x_p)/(fx_n-fx_p);
        x_p=x_n;
        fx_p=fx_n;
        x_n=x_n-fx_n*temp;
        fx_n=f(x_n);
        if(abs(x_n-x_p)<delta || abs(fx_n)<epsilon) # if we got enought of precision
            return (x_n,fx_n,it,SEC_OK);
        end
    end
    return (x_n,fx_n,maxit,SEC_ITERATION_LIMIT);
end

end #end_module
