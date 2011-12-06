This code is a simple stencil based test code for computing the hyperbolic advance component
of an S3D advance.

It contains a driver, an initialization routine and a subroutine hypterm.f

hypterm.f is the relevant piece.

At the end it computes that sum of the squares of the components.

The answer should be (to some reasonable precision) :

 component, fluxmag           1   232.449442195682     
 component, fluxmag           2   78.9403947583534     
 component, fluxmag           3   147.301037115352     
 component, fluxmag           4   85.1946982609069     
 component, fluxmag           5  2.845348677916408E+021
 component, fluxmag           6   8614670.63232723     
 component, fluxmag           7   9683368.91970733     
 component, fluxmag           8   9219846.56270665     
 component, fluxmag           9   40317366.0809986     
 component, fluxmag          10  0.000000000000000E+000
 component, fluxmag          11  0.000000000000000E+000
 component, fluxmag          12  0.000000000000000E+000
 component, fluxmag          13  0.000000000000000E+000
 component, fluxmag          14  0.000000000000000E+000

