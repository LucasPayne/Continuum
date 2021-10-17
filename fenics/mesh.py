from fenics import *
import mshr

# Define domain
center = Point(0.2, 0.2)
radius = 0.05
L = 2.2
W = 0.41
geometry = mshr.Rectangle(Point(0.0, 0.0), Point(L, W)) \
         - mshr.Circle(center, radius, 10)

print(geometry)

f = open("thing", "w+")
f.write(geometry)

