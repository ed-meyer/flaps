# Flaps: FLutter Analysis ProgramS
Flaps is a collection of programs for analyzing aircraft
or windmills for flutter instabilities which can cause a
range of problems from ride-quality to catastrophic failure.

Flutter is an aeroelastic phenomenon caused by the interaction
of a fluid and a flexible structure moving relative to the
fluid. Analyzing aircraft for flutter is an important part
of the design and certification process.

Goals in writing Flaps include
  - give working engineers a tool for designing and certifying
    aircraft that are free from flutter instabilities,
  - is relatively easy to use and integrates well with other
    structural analysis tools,
  - do the many parameter studies necessary very efficiently.
  - treat the multiple nonlinearities present in realistic
    structures

Flutter solutions are done in the frequency domain for efficiency, using
a continuation method and automatic differentiation for flexibility in
types of solution.  A variety of methods for parameterizing matrices
are included using built-in or user-defined parameters and code.
Nonlinearities are approximated with describing functions.
