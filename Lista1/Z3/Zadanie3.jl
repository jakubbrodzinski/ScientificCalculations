#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
    Function checks if in [1,2] Float64 numbers are evenly placed and if every number from [1,2] can be represented as '1+kδ'
    where k is from [1,2^52-1]

    variables:
        a,b = start and end of interval
        x - neighbor value from the x's last value calculated adding δ
        y - neighbor value from the y's last value calculated using nextfloat method
        δ - step used to calculate variable x
        k - iterator
=#

function density()
    a=Float64(1.0)
    b=Float64(2.0)
    x = a
    y = a
    δ = Float64(2.0^-52)
    k = 0
    while x < b
        k+=1
        x = a + Float64(Float64(k)*δ)
        y = nextfloat(y)
        if x != y
            println("W Float64 krok nie jest równy ",δ)
            return
        end
    end
    println("W Float64 krok jest równy ",δ)
    if(Float64(δ(k-1))==Float64(2^52-1))
        println("Każdego Float64 z przedziału [1,2] mozna przedstawic jako x=1.0+ k*δ, gdzie δ=",δ," i k [1,",k-1,"]")
    end
end

#=
    Shows that in [1/2,1] and [2,4] the numbers aren't evenly distributed
=#
println("1/2 + δ = ",bits(Float64(1.0/2.0)+Float64(2.0^-52)))
println("nextfloat(Float64(1/2)) = ",bits(nextfloat(Float64(1.0/2.0))))

println("2.0 + δ = ",bits(Float64(2.0)+Float64(2.0^-52)))
println("nextfloat(Float64(2.0)) = ",bits(nextfloat(Float64(2.0))))

density()
