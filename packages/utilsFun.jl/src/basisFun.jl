function e_0(z::Cdouble,z_min::Cdouble,z_max::Cdouble)
   ((z_max-z)/(z_max-z_min))*(z>=z_min)*(z<=z_max)
end

function e_k(z::Cdouble,z_min::Cdouble,z_mid::Cdouble,z_max::Cdouble)
   max(0.0,min((z-z_min)/(z_mid-z_min),(z_max-z)/(z_max-z_mid)))
end

function e_M(z::Cdouble,z_min::Cdouble,z_max::Cdouble)
   ((z-z_min)/(z_max-z_min))*(z>=z_min)*(z<=z_max)
end
