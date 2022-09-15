using CairoMakie

f = Figure(resolution = (800, 800))

for i in 1:2, j in 1:2
    Axis(f[i, j])
end

Label(f[0, :], text = "First Supertitle", textsize = 20)
Label(f[-1, :], text = "Second Supertitle", textsize = 30)
Label(f[-2, :], text = "Third Supertitle", textsize = 40)

f