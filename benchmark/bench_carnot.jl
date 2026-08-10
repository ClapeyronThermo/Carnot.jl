module CarnotBench

using BenchmarkTools
using Clapeyron
using Carnot

# -------------------------
# Constants
# -------------------------
const p = 2 * 101325.0
const p_out = 5 * p
const p_out_expander = p / 5
const η_isen = 0.8

const fluid = cPR(["Propane", "butane"], idealmodel = ReidIdeal)
const z = [1.0, 1.0]

const h = collect(range(-1, 1, length = 100))
    x = copy(h)
# -------------------------
# Functions
# -------------------------
f_compression(x) =
    ThermoCycleGlides.isentropic_compressor(p, p_out, η_isen, x, z, fluid)

f_expansion(x) =
    ThermoCycleGlides.isentropic_expander(p, p_out_expander, η_isen, x, z, fluid)

# -------------------------
# Benchmark suite
# # -------------------------
suite = BenchmarkGroup()

suite["utils"] = BenchmarkGroup()

suite["utils"]["compressor"] = @benchmarkable begin
    f_compression.($x)
end

suite["utils"]["expander"] =
    @benchmarkable f_expansion.(x) setup=(x = copy(h),)


end

# return benchmark suite
CarnotBench.suite