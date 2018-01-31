#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#
using Polynomials;

#=
  Polynomial that we work with is called "Wilkinson's polynomial" and after this experiment we will be able to see why it's so special.
  We have to import library to be able to call it's functions.
  Variables:
    input - array with polynomial's coefficients
    input_2 - array with slightly modified polynomial's coefficients (we repeat Wilkinson's experiment)
    x0_k - array with precise polynomial's roots
    P - polynomial created from a_n * x^n +... + a_1 * x^1 + a_0 presentation
    p - polunomial created from (x-x0_n)(x-x0_n-1)...(x_x0_1) presentation
    z_roots - array with roots calculated using library's function
=#
input=Array{Float64,1}([1, -210.0, 20615.0,-1256850.0, 53327946.0,-1672280820.0, 40171771630.0, -756111184500.0, 11310276995381.0, -135585182899530.0,
 1307535010540395.0, -10142299865511450.0, 63030812099294896.0, -311333643161390640.0, 1206647803780373360.0, -3599979517947607200.0,
  8037811822645051776.0,-12870931245150988800.0, 13803759753640704000.0, -8752948036761600000.0, 2432902008176640000.0]);

input_2=Array{Float64,1}([1, Float64(-210.0-2.0^-23), 20615.0,-1256850.0, 53327946.0,-1672280820.0, 40171771630.0, -756111184500.0, 11310276995381.0, -135585182899530.0,
   1307535010540395.0, -10142299865511450.0, 63030812099294896.0, -311333643161390640.0, 1206647803780373360.0, -3599979517947607200.0,
    8037811822645051776.0,-12870931245150988800.0, 13803759753640704000.0, -8752948036761600000.0, 2432902008176640000.0]);

x0_k=Array{Float64,1}(20);
for i in 1:20
  x0_k[i]=i;
end

P=Poly(flipdim(input,1));
p=poly(x0_k);

z_roots= roots(P);

#=
  We present results via simple for loop.
=#

for k in 1:20
  println("k = ",(20-k+1),"\tz_k= ",z_roots[k]);
  println("|P(z_k)|= ",abs(polyval(P,z_roots[k])),"\t|p(z_k)|= ",abs(polyval(p,z_roots[k])),"\t|z_k-k|= ",abs(z_roots[k]-(20-k+1)));
end

#=
  We repeat loop for little bit modified input (input_2)
=#

P_2=Poly(flipdim(input_2,1));
z_roots_2= roots(P_2);

#=
  We present results via simple for loop.
=#

for k in 1:20
  println("k = ",(20-k+1),"\tz_k= ",z_roots_2[k]);
  println("|P(z_k)|= ",abs(polyval(P_2,z_roots_2[k])),"\t|p(z_k)|= ",abs(polyval(p,z_roots_2[k])),"\t|z_k-k|= ",abs(z_roots_2[k]-(20-k+1)));
end
