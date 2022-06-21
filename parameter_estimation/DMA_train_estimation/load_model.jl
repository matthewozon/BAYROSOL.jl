# the measurement model (from density to counts, must be modified so that it goes from concentration to counts)
#WARNING:it maps from dNdlo10Dp to counts, so the modification is just a multiplcation by log10(cst_r_model)
if FLAG_1952_02
    # H_DMATRAIN     = Matrix{Cdouble}(DataFrame(CSV.File(string(input_folder,"measurement_operator_dmatrain_h2so4_v2.csv"); header=false))); # delim=",",
    df_H           = CSV.File(string(input_folder,"measurement_operator_dmatrain_1952_02.csv"); header=false) |> DataFrame
    H_DMATRAIN     = Matrix{Cdouble}(df_H);
end
if FLAG_0000_00
    df_H           = CSV.File(string(input_folder,"measurement_operator_dmatrain_0000_00.csv"); header=false) |> DataFrame
    H_DMATRAIN     = Matrix{Cdouble}(df_H);
end
if FLAG_0000_01
    df_H           = CSV.File(string(input_folder,"measurement_operator_dmatrain_0000_01.csv"); header=false) |> DataFrame
    H_DMATRAIN     = Matrix{Cdouble}(df_H);
end
if FLAG_0000_02
    df_H           = CSV.File(string(input_folder,"measurement_operator_dmatrain_0000_02.csv"); header=false) |> DataFrame
    H_DMATRAIN     = Matrix{Cdouble}(df_H);
end
if FLAG_0000_03
    df_H           = CSV.File(string(input_folder,"measurement_operator_dmatrain_0000_03.csv"); header=false) |> DataFrame
    H_DMATRAIN     = Matrix{Cdouble}(df_H);
end
if FLAG_0000_04
    df_H           = CSV.File(string(input_folder,"measurement_operator_dmatrain_0000_04.csv"); header=false) |> DataFrame
    H_DMATRAIN     = Matrix{Cdouble}(df_H);
end

if FLAG_1802_01
    df_H           = CSV.File(string(input_folder,"measurement_operator_dmatrain_1802_01.csv"); header=false) |> DataFrame
    H_DMATRAIN     = Matrix{Cdouble}(df_H);
end
if FLAG_1906_03
    df_H           = CSV.File(string(input_folder,"measurement_operator_dmatrain_1906_03.csv"); header=false) |> DataFrame
    H_DMATRAIN     = Matrix{Cdouble}(df_H);
end


if FLAG_0000_00
    df_diam_mod = CSV.File(string(input_folder,"measurement_operator_discretization_0000_00.csv"); header=false) |> DataFrame
    diameter_model = dropdims(Matrix{Cdouble}(df_diam_mod),dims=1);
end
if FLAG_0000_01
    df_diam_mod = CSV.File(string(input_folder,"measurement_operator_discretization_0000_01.csv"); header=false) |> DataFrame
    diameter_model = dropdims(Matrix{Cdouble}(df_diam_mod),dims=1);
end
if FLAG_0000_02
    df_diam_mod = CSV.File(string(input_folder,"measurement_operator_discretization_0000_02.csv"); header=false) |> DataFrame
    diameter_model = dropdims(Matrix{Cdouble}(df_diam_mod),dims=1);
end
if FLAG_0000_03
    df_diam_mod = CSV.File(string(input_folder,"measurement_operator_discretization_0000_03.csv"); header=false) |> DataFrame
    diameter_model = dropdims(Matrix{Cdouble}(df_diam_mod),dims=1);
end
if FLAG_0000_04
    df_diam_mod = CSV.File(string(input_folder,"measurement_operator_discretization_0000_04.csv"); header=false) |> DataFrame
    diameter_model = dropdims(Matrix{Cdouble}(df_diam_mod),dims=2);
end


if FLAG_1952_02
    df_diam_mod = CSV.File(string(input_folder,"measurement_operator_discretization_1952_02.csv"); header=false) |> DataFrame
    diameter_model = 1.0e-9dropdims(Matrix{Cdouble}(df_diam_mod),dims=2);
end

if FLAG_1802_01
    df_diam_mod = CSV.File(string(input_folder,"measurement_operator_discretization_1802_01.csv"); header=false) |> DataFrame
    diameter_model = 1.0e-9dropdims(Matrix{Cdouble}(df_diam_mod),dims=2);
end

if FLAG_1906_03
    df_diam_mod = CSV.File(string(input_folder,"measurement_operator_discretization_1906_03.csv"); header=false) |> DataFrame
    diameter_model = 1.0e-9dropdims(Matrix{Cdouble}(df_diam_mod),dims=2);
end
