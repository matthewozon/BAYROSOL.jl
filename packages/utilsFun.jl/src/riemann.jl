## horrible way to compute an integral
function riemann(f::Function, a::Real, b::Real, n::Int; method="right")
  if method == "right"
     xs = a .+ collect(0.0:n) * (b-a)/n
     # as = [meth(f, l, r) for (l,r) in zip(xs[1:end-1], xs[2:end])]
     as = [f(r)*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
  elseif method == "left"
     # meth(f,l,r) = f(l) * (r-l)
     xs = a .+ collect(0.0:n) * (b-a)/n
     as = [f(l)*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
  elseif method == "trapezoid"
     # meth(f,l,r) = (1/2) * (f(l) + f(r)) * (r-l)
     xs = a .+ collect(0.0:n) * (b-a)/n
     as = [(1.0/2.0)*(f(l) + f(r))*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
  elseif method == "simpsons"
     # meth(f,l,r) = (1.0/6.0) * (f(l) + 4.0*(f((l+r)/2.0)) + f(r)) * (r-l)
     xs = a .+ collect(0.0:n) * (b-a)/n
     as = [(1.0/6.0) * (f(l) + 4.0*(f((l+r)/2.0)) + f(r)) * (r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
 else
     throw(@sprintf "quadrature %s is not implemented" method)
  end
  sum(as)
end
