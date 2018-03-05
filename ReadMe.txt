Assignment #3: Ray tracing

FULL NAME: Darwin Mendyke

MANDATORY FEATURES
------------------
Feature:                                 Status: finish? (yes/no)
-------------------------------------    -------------------------
1) Ray tracing triangles                  yes

2) Ray tracing sphere                     yes

3) Triangle Phong Shading                 yes

4) Sphere Phong Shading                   yes

5) Shadows rays                           yes

6) Still images                           yes

7) Extra Credit (up to 20 points)

I implemented the recursive reflection tracing.
This was done by taking the reflected ray based off the view vector and the normal for a given point
then recursively calling the tracer function with this ray. The direction of the ray is calculated using the reflection
equation R = 2 * (V * N) * N - V and the ray origin is merely the intersection point offset slightly by the normal of that point.

The tracer function will continue to recursively call itself with new reflection rays up to two times, assuming
the ray keeps reflecting with more surfaces.

To run the assignment with the recursive reflection tracing, comment out line 274, then uncomment lines 275 - 287


My screenshots correspond to the following images.
000.jpg -> test1.scene
001.jpg -> test2.scene
002.jpg -> spheres.scene
003.jpg -> table.scene
004.jpg -> siggraph.scene

005.jpg -> test2.scene but with reflection
006.jpg -> spheres.scene but with reflection
007.jpg -> custom scene with reflection
008.jpg -> custom scene with reflection
