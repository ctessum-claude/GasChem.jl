@testitem "Pollu - Case1 - 60s" begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit

    tspan = (0.0, 60.0)
    rs = mtkcompile(Pollu())

    # Initial conditions from the paper (Verwer 1994, Appendix), converted from ppm to ppb
    initial_conditions = [
        rs.NO2 => 0,
        rs.NO => 0.2 * 1e3,
        rs.O3P => 0,
        rs.O3 => 0.04 * 1e3,
        rs.HO2 => 0,
        rs.OH => 0,
        rs.CH2O => 0.1 * 1e3,
        rs.CO => 0.3 * 1e3,
        rs.ALD => 0.01 * 1e3,
        rs.MEO2 => 0,
        rs.C2O3 => 0,
        rs.CO2 => 0,
        rs.PAN => 0,
        rs.CH3O => 0,
        rs.HNO3 => 0,
        rs.O1D => 0,
        rs.SO2 => 0.007 * 1e3,
        rs.SO4 => 0,
        rs.NO3 => 0,
        rs.N2O5 => 0,
    ]

    sol = solve(
        ODEProblem(rs, initial_conditions, tspan, []),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1e-12,
        reltol = 1e-12,
    )

    # Reference: O3 at t=1 min = 0.32994065756620E-02 ppm (Verwer 1994, Table, t=1 min)
    @test sol[rs.O3][end] ≈ 0.32994065756620e-2 * 1e3 rtol = 0.001
end

@testitem "Pollu - Case2 - 3600s" begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit

    tspan = (0.0, 3600.0)
    rs = mtkcompile(Pollu())

    # Initial conditions from the paper (Verwer 1994, Appendix), converted from ppm to ppb
    initial_conditions = [
        rs.NO2 => 0,
        rs.NO => 0.2 * 1e3,
        rs.O3P => 0,
        rs.O3 => 0.04 * 1e3,
        rs.HO2 => 0,
        rs.OH => 0,
        rs.CH2O => 0.1 * 1e3,
        rs.CO => 0.3 * 1e3,
        rs.ALD => 0.01 * 1e3,
        rs.MEO2 => 0,
        rs.C2O3 => 0,
        rs.CO2 => 0,
        rs.PAN => 0,
        rs.CH3O => 0,
        rs.HNO3 => 0,
        rs.O1D => 0,
        rs.SO2 => 0.007 * 1e3,
        rs.SO4 => 0,
        rs.NO3 => 0,
        rs.N2O5 => 0,
    ]

    sol = solve(
        ODEProblem(rs, initial_conditions, tspan, []),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1e-12,
        reltol = 1e-12,
    )

    # Reference: O3 at t=60 min = 0.552314020747798E-02 ppm (Verwer 1994, Table, t=60 min)
    @test sol[rs.O3][end] ≈ 0.552314020747798e-2 * 1e3 rtol = 0.001
end

@testitem "Pollu - Structural Verification" begin
    using ModelingToolkit

    rs = mtkcompile(Pollu())
    # 20 species from the paper
    @test length(unknowns(rs)) == 20
end

@testitem "Pollu - Reference solution all species at t=60s" begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit

    tspan = (0.0, 60.0)
    rs = mtkcompile(Pollu())

    # Initial conditions from the paper (Verwer 1994, Appendix), converted from ppm to ppb
    initial_conditions = [
        rs.NO2 => 0,
        rs.NO => 0.2 * 1e3,
        rs.O3P => 0,
        rs.O3 => 0.04 * 1e3,
        rs.HO2 => 0,
        rs.OH => 0,
        rs.CH2O => 0.1 * 1e3,
        rs.CO => 0.3 * 1e3,
        rs.ALD => 0.01 * 1e3,
        rs.MEO2 => 0,
        rs.C2O3 => 0,
        rs.CO2 => 0,
        rs.PAN => 0,
        rs.CH3O => 0,
        rs.HNO3 => 0,
        rs.O1D => 0,
        rs.SO2 => 0.007 * 1e3,
        rs.SO4 => 0,
        rs.NO3 => 0,
        rs.N2O5 => 0,
    ]

    sol = solve(
        ODEProblem(rs, initial_conditions, tspan, []),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1e-12,
        reltol = 1e-12,
    )

    # Reference values from Verwer 1994, Table at t=1 min, converted from ppm to ppb
    @test sol[rs.NO2][end] ≈ 0.37326304298606E-01 * 1e3 rtol = 0.01
    @test sol[rs.NO][end] ≈ 0.16251325412704E+00 * 1e3 rtol = 0.01
    @test sol[rs.O3][end] ≈ 0.32994065756620E-02 * 1e3 rtol = 0.01
    @test sol[rs.HO2][end] ≈ 0.31151619385738E-06 * 1e3 rtol = 0.03 # Radical species: larger tolerance due to fast dynamics
    @test sol[rs.OH][end] ≈ 0.26534918180427E-06 * 1e3 rtol = 0.01
    @test sol[rs.CH2O][end] ≈ 0.99423103666384E-01 * 1e3 rtol = 0.01
    @test sol[rs.CO][end] ≈ 0.30061731277555E+00 * 1e3 rtol = 0.01
end

@testitem "Pollu - Reference solution all species at t=3600s" begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit

    tspan = (0.0, 3600.0)
    rs = mtkcompile(Pollu())

    # Initial conditions from the paper (Verwer 1994, Appendix), converted from ppm to ppb
    initial_conditions = [
        rs.NO2 => 0,
        rs.NO => 0.2 * 1e3,
        rs.O3P => 0,
        rs.O3 => 0.04 * 1e3,
        rs.HO2 => 0,
        rs.OH => 0,
        rs.CH2O => 0.1 * 1e3,
        rs.CO => 0.3 * 1e3,
        rs.ALD => 0.01 * 1e3,
        rs.MEO2 => 0,
        rs.C2O3 => 0,
        rs.CO2 => 0,
        rs.PAN => 0,
        rs.CH3O => 0,
        rs.HNO3 => 0,
        rs.O1D => 0,
        rs.SO2 => 0.007 * 1e3,
        rs.SO4 => 0,
        rs.NO3 => 0,
        rs.N2O5 => 0,
    ]

    sol = solve(
        ODEProblem(rs, initial_conditions, tspan, []),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1e-12,
        reltol = 1e-12,
    )

    # Reference values from Verwer 1994, Table at t=60 min, converted from ppm to ppb
    @test sol[rs.NO2][end] ≈ 0.56462554800124E-01 * 1e3 rtol = 0.01
    @test sol[rs.NO][end] ≈ 0.13424841304232E+00 * 1e3 rtol = 0.01
    @test sol[rs.O3][end] ≈ 0.552314020747798E-02 * 1e3 rtol = 0.01
    @test sol[rs.HO2][end] ≈ 0.20189772623092E-06 * 1e3 rtol = 0.03 # Radical species: larger tolerance due to fast dynamics
    @test sol[rs.OH][end] ≈ 0.14645418635004E-06 * 1e3 rtol = 0.01
    @test sol[rs.CH2O][end] ≈ 0.77842491190039E-01 * 1e3 rtol = 0.01
    @test sol[rs.CO][end] ≈ 0.32450753533953E+00 * 1e3 rtol = 0.01
end

@testitem "Couple Pollu and FastJX" begin
    using EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using ModelingToolkit

    sf1 = couple(Pollu(), FastJX_interpolation_troposphere(0.0))
    sf2 = couple(Pollu(), FastJX(0.0))
    sys1 = convert(System, sf1)
    sys2 = convert(System, sf2)
    tspan = (0.0, 3600 * 24)
    prob1 = ODEProblem(sys1, [], tspan, [])
    prob2 = ODEProblem(sys2, [], tspan, [])
    sol1 = solve(prob1, Rosenbrock23(), saveat = 10.0)
    sol2 = solve(prob2, Rosenbrock23(), saveat = 10.0)
    @test sol1[sys1.Pollu₊O3][4320] ≈ sol2[sys2.Pollu₊O3][4320] rtol = 1e-4
end

@testitem "Couple Pollu and NEI2016MonthlyEmis" begin
    using EarthSciMLBase, EarthSciData
    using ModelingToolkit
    using Dates

    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 4);
        lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
        latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
        levrange = 1:15,
    )

    model_3way = couple(
        FastJX(get_tref(domain)),
        Pollu(),
        NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain),
    )

    sys = convert(System, model_3way)
    @test length(unknowns(sys)) == 20

    eqs = string(equations(sys))

    wanteq = "Differential(t)(Pollu₊ALD(t)) ~ Pollu₊NEI2016MonthlyEmis_ALD2(t)"
    @test contains(string(eqs), wanteq)
end
