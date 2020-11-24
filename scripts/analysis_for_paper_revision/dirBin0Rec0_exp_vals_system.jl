using SymPy

x,y,a,b = symbols("x y a b")

Eq(a,(x+x^2*y)/(1 + 2*x + x^2*y) )
Eq(b,(x^2*y)/(1 + 2*x + x^2*y) )
