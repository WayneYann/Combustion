C     CVS $Revision: 1.1.1.1 $ repositied $Date: 2006/05/26 19:09:32 $

C///////////////////////////////////////////////////////////////////////
C
C     CHEMKIN-III file cdassl.f
C
C     a modification of DASSL
C
C
C     94
C
C        1. Original modifications due to Harry Moffat.
C
C     97/10/29  Joseph Grcar
C
C        1. Replaced SIGN by SIGN77.
C
C///////////////////////////////////////////////////////////////////////

C@(#)===================================================================
C@(#)
C@(#)
C@(#)                   FILE =  ddassl.f
C@(#)
C@(#)  ---------------  VERSION = 1.2
C@(#)  |  SCCS  FILE |
C@(#)  |   SUMMARY   |  CURRENT CHECKOUT DATE = 08/11/94
C@(#)  ---------------                           at 14:55:34
C@(#)                   DATE OF NEWEST DELTA = 08/11/94
C@(#)                                            at 14:55:33
C@(#)  SCCS file name = /users/chemkin/SCCS/s.ddassl.f
C@(#)===================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE DDASSL (RES, NEQ, T, Y, YPRIME, TOUT,
     *                   INFO, RTOL, ATOL, IDID, RWORK, LRW, IWORK,
     *                   LIW, RPAR, IPAR, JAC)
C
C-----------------------------------------------------------------------
C
C  This is the september 1982 version of dassl --
C  differential/algebraic systems solver.
C
C  This code was written by
C          linda r. petzold
C          applied mathematics division 8331
C          sandia national laboratories
C          livermore, ca       94550
C
C
C  This code solves a system of differential/
C  algebraic equations of the form
C  g(t,y,yprime) = 0.
C
C  subroutine dassl uses the backward
C  differentiation formulas of orders one
C  through five to solve a system of the above
C  form for y and yprime. values for y
C  and yprime at the initial time must
C  be given as input. these values must
C  be consistent, (that is. if t,y,yprime
C  are the given initial values, they must
C  satisfy g(t,y,yprime) = 0.)
C  the subroutine solves the system from t to tout. it is
C  easy to continue the solution to get results
C  at additional tout. this is the interval
C  mode of operation. intermediate results can
C  also be obtained easily by using the intermediate-
C  output capability.
C
C  dassl uses subroutines dastep,inider,vnorm,solve,
C  njac,intrp,d1mach,wtset
C  and the linpack routines sgefa,sgesl,
C  sgbsl,and the blas routines saxpy,sscal,
C  isamax,sdot,and routines from the error-
C  handling package xerror.  the routine xerrwv from
C  xerror is machine-dependent.
C
C  a labelled common nwm001 is used internally to
C  communicate pointers to locations in rwork and iwork.
C
C
C  ------------description of arguments to dassl------------------------
C  ------------(an overview)--------------------------------------------
C
C  the parameters are
C
C  RES -- this is a subroutine which you provide
C         to define the differential/algebraic
C         system
C
C  NEQ -- this is the number of equations
C         to be solved
C
C  T -- this is the current value of the
C       independent variable.
C
C  TOUT -- this is a point at which a solution
C      is desired.
C
C  INFO(*) -- the basic task of the code is
C             to solve the system from t to
C             tout and return an answer at tout.
C             info(*) is an integer array which is
C             used to communicate exactly how you
C             want this task to be carried out.
C
C  Y(*) -- this array contains the solution
C          components at t
C
C  YPRIME(*) -- this array contains the derivatives
C               of the solution components at t
C
C  RTOL,ATOL -- these quantities represent
C               absolute and relative error
C               tolerances which you provide to indicate
C               how accurately you wish the solution
C               to be computed. you may choose them
C               to be both scalars or else both
C               vectors.
C
C  IDID -- this scalar quantity is an indicator reporting
C          what the code did. you must monitor this
C          integer variable to decide what action to
C          take next.
C
C  RWORK(*),LRW -- rwork(*) is a real work array of
C                  length lrw which provides the code
C                  with needed storage space.
C
C  IWORK(*),LIW -- iwork(*) is an integer work array
C                  of length liw which provides the code
C                  with needed storage space.
C
C  RPAR,IPAR -- these are real and integer parameter
C               arrays which you can use for
C               communication between your calling
C               program and the res subroutine
C               (and the jac subroutine)
C
C  JAC -- this is the name of a subroutine which you
C         may choose to provide for defining
C         a matrix of partial derivatives
C         described below.
C
C  quantities which are used as input items are
C     NEQ,T,Y(*),YPRIME(*),TOUT,INFO(*),
C     RTOL,ATOL,RWORK(1),RWORK(2),RWORK(3),LRW,IWORK(1),
C     IWORK(2),IWORK(3),AND LIW.
C
C  quantities which may be altered by the code are
C     T,Y(*),YPRIME(*),INFO(1),RTOL,ATOL,
C     IDID,RWORK(*) AND IWORK(*)
C
C------------input-what to do on the first call to dassl----------------
C
C
C  the first call of the code is defined to be the start of each new
C  problem. read through the descriptions of all the following items,
C  provide sufficient storage space for designated arrays, set
C  appropriate variables for the initialization of the problem, and
C  give information about how you want the problem to be solved.
C
C
C  RES -- provide a subroutine of the form
C             subroutine res(t,y,yprime,delta,ires,rpar,ipar)
C         to define the system of differential/algebraic
C         equations which is to be solved. for the given values
C         of t,y and yprime, the subroutine should
C         return the residual of the differential/algebraic
C         system
C             delta = g(t,y,yprime)
C         (delta(*) is a vector of length neq which is
C         output for res.)
C
C         subroutine res must not alter t,y or yprime.
C         you must declare the name res in an external
C         statement in your program that calls dassl.
C         you must dimension y,yprime and delta in res.
C
C         ires is an integer flag which is always equal to
C         zero on input.  subroutine res should alter ires
C         only if it encounters an illegal value of y or
C         a stop condition.  set ires = -1 if an input value
C         is illegal, and dassl will try to solve the problem
C         without getting ires = -1.  if ires = -2, dassl
C         will return control to the calling program
C         with idid = -11.
C
C         rpar and ipar are real and integer parameter arrays which
C         you can use for communication between your calling program
C         and subroutine res. they are not altered by dassl. if you
C         do not need rpar or ipar, ignore these parameters by treat-
C         ing them as dummy arguments. if you do choose to use them,
C         dimension them in your calling program and in res as arrays
C         of appropriate length.
C
C  NEQ -- set it to the number of differential equations.
C         (neq .ge. 1)
C
C  t -- set it to the initial point of the integration.
C       t must be defined as a variable.
C
C  Y(*) -- set this vector to the initial values of the neq solution
C          components at the initial point. you must dimension y of
C          length at least neq in your calling program.
C
C  YPRIME(*) -- set this vector to the initial values of
C               the neq first derivatives of the solution
C               components at the initial point. you
C               must dimension yprime at least neq
C               in your calling program.  if you do not
C               know initial values of some of the solution
C               components, see the explanation of info(11).
C
C  TOUT - set it to the first point at which a solution
C         is desired. you can not take tout = t.
C         integration either forward in t (tout .gt. t) or
C         backward in t (tout .lt. t) is permitted.
C
C         the code advances the solution from t to tout using
C         step sizes which are automatically selected so as to
C         achieve the desired accuracy. if you wish, the code will
C         return with the solution and its derivative at
C         intermediate steps (intermediate-output mode) so that
C         you can monitor them, but you still must provide tout in
C         accord with the basic aim of the code.
C
C         the first step taken by the code is a critical one
C         because it must reflect how fast the solution changes near
C         the initial point. the code automatically selects an
C         initial step size which is practically always suitable for
C         the problem. by using the fact that the code will not step
C         past tout in the first step, you could, if necessary,
C         restrict the length of the initial step size.
C
C         for some problems it may not be permissable to integrate
C         past a point tstop because a discontinuity occurs there
C         or the solution or its derivative is not defined beyond
C         tstop. when you have declared a tstop point (see info(4)
C         and rwork(1)), you have told the code not to integrate
C         past tstop. in this case any tout beyond tstop is invalid
C         input.
C
C  info(*) - use the info array to give the code more details about
C            how you want your problem solved. this array should be
C            dimensioned of length 15, though dassl uses
C            only the first nine entries. you must respond to all of
C            the following items which are arranged as questions. the
C            simplest use of the code corresponds to answering all
C            questions as yes ,i.e. setting all entries of info to 0.
C
C       info(1) - this parameter enables the code to initialize
C              itself. you must set it to indicate the start of every
C              new problem.
C
C          **** is this the first call for this problem ...
C                yes - set info(1) = 0
C                 no - not applicable here.
C                      see below for continuation calls.  ****
C
C       info(2) - how much accuracy you want of your solution
C              is specified by the error tolerances rtol and atol.
C              the simplest use is to take them both to be scalars.
C              to obtain more flexibility, they can both be vectors.
C              the code must be told your choice.
C
C          **** are both error tolerances rtol, atol scalars ...
C                yes - set info(2) = 0
C                      and input scalars for both rtol and atol
C                 no - set info(2) = 1
C                      and input arrays for both rtol and atol ****
C
C       info(3) - the code integrates from t in the direction
C              of tout by steps. if you wish, it will return the
C              computed solution and derivative at the next
C              intermediate step (the intermediate-output mode) or
C              tout, whichever comes first. this is a good way to
C              proceed if you want to see the behavior of the solution.
C              if you must have solutions at a great many specific
C              tout points, this code will compute them efficiently.
C
C          **** do you want the solution only at
C                tout (and not at the next intermediate step) ...
C                 yes - set info(3) = 0
C                  no - set info(3) = 1 ****
C
C       info(4) - to handle solutions at a great many specific
C              values tout efficiently, this code may integrate past
C              tout and interpolate to obtain the result at tout.
C              sometimes it is not possible to integrate beyond some
C              point tstop because the equation changes there or it is
C              not defined past tstop. then you must tell the code
C              not to go past.
C
C           **** can the integration be carried out without any
C                restrictions on the independent variable t ...
C                 yes - set info(4)=0
C                  no - set info(4)=1
C                       and define the stopping point tstop by
C                       setting rwork(1)=tstop ****
C
C       info(5) - to solve differential/algebraic problems it is
C              necessary to use a matrix of partial derivatives of the
C              system of differential equations.  if you do not
C              provide a subroutine to evaluate it analytically (see
C              description of the item jac in the call list), it will
C              be approximated by numerical differencing in this code.
C              although it is less trouble for you to have the code
C              compute partial derivatives by numerical differencing,
C              the solution will be more reliable if you provide the
C              derivatives via jac. sometimes numerical differencing
C              is cheaper than evaluating derivatives in jac and
C              sometimes it is not - this depends on your problem.
C
C           **** do you want the code to evaluate the partial
C                  derivatives automatically by numerical differences
C                   yes - set info(5)=0
C                    no - set info(5)=1
C                  and provide subroutine jac for evaluating the
C                  matrix of partial derivatives ****
C
C       info(6) - dassl will perform much better if the matrix of
C              partial derivatives, dg/dy + cj*dg/dyprime,
C              (here cj is a scalar determined by dassl)
C              is banded and the code is told this. in this
C              case, the storage needed will be greatly reduced,
C              numerical differencing will be performed much cheaper,
C              and a number of important algorithms will execute much
C              faster. the differential equation is said to have
C              half-bandwidths ml (lower) and mu (upper) if equation i
C              involves only unknowns y(j) with
C                             i-ml .le. j .le. i+mu
C              for all i=1,2,...,neq. thus, ml and mu are the widths
C              of the lower and upper parts of the band, respectively,
C              with the main diagonal being excluded. if you do not
C              indicate that the equation has a banded matrix of partial
C                 derivatives
C              the code works with a full matrix of neq**2 elements
C              (stored in the conventional way). computations with
C              banded matrices cost less time and storage than with
C              full matrices if  2*ml+mu .lt. neq.  if you tell the
C              code that the matrix of partial derivatives has a banded
C              structure and you want to provide subroutine jac to
C              compute the partial derivatives, then you must be careful
C              to store the elements of the matrix in the special form
C              indicated in the description of jac.
C
C          **** do you want to solve the problem using a full
C               (dense) matrix (and not a special banded
C               structure) ...
C                yes - set info(6)=0
C                 no - set info(6)=1
C                       and provide the lower (ml) and upper (mu)
C                       bandwidths by setting
C                       iwork(1)=ml
C                       iwork(2)=mu ****
C
C
C        info(7) -- you can specify a maximum (absolute value of)
C              stepsize, so that the code
C              will avoid passing over very
C              large regions.
C
C          ****  do you want the code to decide
C                on its own maximum stepsize?
C                yes - set info(7)=0
C                 no - set info(7)=1
C                      and define hmax by setting
C                      rwork(2)=hmax ****
C
C        info(8) -- differential/algebraic problems
C              may occaisionally suffer from
C              severe scaling difficulties on the
C              first step. if you know a great deal
C              about the scaling of your problem, you can
C              help to alleviate this problem by
C              specifying an initial stepsize ho.
C
C          ****  do you want the code to define
C                its own initial stepsize?
C                yes - set info(8)=0
C                 no - set info(8)=1
C                      and define ho by setting
C                      rwork(3)=ho ****
C
C        info(9) -- if storage is a severe problem,
C              you can save some locations by
C              restricting the maximum order maxord.
C              the default value is 5. for each
C              order decrease below 5, the code
C              requires neq fewer locations, however
C              it is likely to be slower. in any
C              case, you must have 1 .le. maxord .le. 5
C          ****  do you want the maximum order to
C                default to 5?
C                yes - set info(9)=0
C                 no - set info(9)=1
C                      and define maxord by setting
C                      iwork(3)=maxord ****
C
C        info(10) --if you know that the solutions to your
C               equations will
C               always be nonnegative, it may help to set this
C               parameter.  however, it is probably best to
C               try the code without using this option first,
C               and only to use this option if that doesn't
C               work very well.
C           ****  do you want the code to solve the problem without
C                 invoking any special nonnegativity constraints?
C                  yes - set info(10)=0
C                   no - set info(10)=1
C
C        info(11) --dassl normally requires the initial t,
C               y, and yprime to be consistent.  that is,
C               you must have g(t,y,yprime) = 0 at the initial
C               time.  if you do not know the initial
C               derivative precisely, you can let dassl try
C               to compute it.
C          ****   are the initial t, y, yprime consistent?
C                 yes - set info(11) = 0
C                  no - set info(11) = 1,
C                       and set yprime to an initial approximation
C                       to yprime.  (if you have no idea what
C                       yprime should be, set it to zero. note
C                       that the initial y should be such
C                       that there must exist a yprime so that
C                       g(t,y,yprime) = 0.)
C
C   rtol, atol -- you must assign relative (rtol) and absolute (atol
C              error tolerances to tell the code how accurately
C              you want
C              the solution to be computed. they must be defined as
C              variables because the code may change them. you have two
C              choices --
C                    both rtol and atol are scalars. (info(2)=0)
C                    both rtol and atol are vectors. (info(2)=1)
C              in either case all components must be non-negative.
C
C              the tolerances are used by the code in a local error
C              test
C              at each step which requires roughly that
C                    abs(local error) .le. rtol*abs(y)+atol
C              for each vector component.
C              (more specifically, a root-mean-square norm is used to
C              measure the size of vectors, and the error test uses the
C              magnitude of the solution at the beginning of the step.)
C
C              the true (global) error is the difference between the
C              true
C              solution of the initial value problem and the computed
C              approximation. practically all present day codes.
C              including this one, control the local error at each step
C              and do not even attempt to control the global error
C              directly.
C              usually, but not always, the true accuracy of
C             the computed y is comparable to the error tolerances. this
C             code will usually, but not always, deliver a more accurate
C             solution if you reduce the tolerances and integrate again.
C              by comparing two such solutions you can get a fairly
C              reliable idea of the true error in the solution at the
C              bigger tolerances.
C
C              setting atol=0. results in a pure relative error test on
C              that component. setting rtol=0. results in a pure
C              absolute
C              error test on that component. a mixed test with non-zero
C              rtol and atol corresponds roughly to a relative error
C              test when the solution component is much bigger than atol
C              and to an absolute error test when the solution component
C              is smaller than the threshold atol.
C
C              the code will not attempt to compute a solution at an
C              accuracy unreasonable for the machine being used. it will
C              advise you if you ask for too much accuracy and inform
C              you as to the maximum accuracy it believes possible.
C
C  RWORK(*) -- dimension this real work array of length lrw in your
C               calling program.
C
C  lrw -- set it to the declared length of the rwork array.
C              you must have
C                   lrw .ge. 40+(maxord+4)*neq+neq**2
C              for the full (dense) jacobian case (when info(6)=0),  or
C                   lrw .ge. 40+(maxord+4)*neq+(2*ml+mu+1)*neq
C              for the banded user-defined jacobian case
C              (when info(5)=1 and info(6)=1), or
C                    lrw .ge. 40+(maxord+4)*neq+(2*ml+mu+1)*neq
C                          +2*(neq/(ml+mu+1)+1)
C              for the banded finite-difference-generated jacobian case
C              (when info(5)=0 and info(6)=1)
C
C  IWORK(*) -- dimension this integer work array of length liw in
C             your calling program.
C
C  liw -- set it to the declared length of the iwork array.
C               you must have liw .ge. 20+neq
C
C  RPAR, IPAR -- these are parameter arrays, of real and integer
C              type, respectively. you can use them for communication
C              between your program that calls dassl and the
C              res subroutine (and the jac subroutine). they are not
C              altered by dassl. if you do not need rpar or ipar, ignore
C              these parameters by treating them as dummy arguments. if
C              you do choose to use them, dimension them in your calling
C              program and in res (and in jac) as arrays of appropriate
C              length.
C
C  JAC -- if you have set info(5)=0, you can ignore this parameter
C              by treating it as a dummy argument. otherwise, you must
C               provide a subroutine of the form
C               jac(t,y,yprime,pd,cj,rpar,ipar)
C               to define the matrix of partial derivatives
C               pd=dg/dy+cj*dg/dyprime
C               cj is a scalar which is input to jac.
C               for the given values of t,y,yprime, the
C               subroutine must evaluate the non-zero partial
C               derivatives for each equation and each solution
C               compowent, and store these values in the
C               matrix pd. the elements of pd are set to zero
C               before each call to jac so only non-zero elements
C               need to be defined.
C
C              subroutine jac must not alter t,y,(*),yprime(*),or cj.
C               you must declare the name jac in an
C               external statement in your program that calls
C               dassl. you must dimension y, yprime and pd
C               in jac.
C
C              the way you must store the elements into the pd matrix
C              depends on the structure of the matrix which you
C              indicated by info(6).
C              *** info(6)=0 -- full (dense) matrix ***
C                  when you evaluate the (non-zero) partial derivative
C                  of equation i with respect to variable j, you must
C               store it in pd according to
C                 pd(i,j) = * dg(i)/dy(j)+cj*dg(i)/dyprime(j)*
C              *** info(6)=1 -- banded jacobian with ml lower and mu
C                  upper diagonal bands (refer to info(6) description of
C                  ml and mu) ***
C                  when you evaluate the (non-zero) partial derivative
C                  of equation i with respect to variable j, you must
C                  store it in pd according to
C                  irow = i - j + ml + mu + 1
C                  pd(irow,j) = *dg(i)/dy(j)+cj*dg(i)/dyprime(j)*
C              rpar and ipar are real and integer parameter arrays which
C              you can use for communication between your calling
C              program and your jacobian subroutine jac. they are not
C              altered by dassl. if you do not need rpar or ipar, ignore
C              these parameters by treating them as dummy arguments. if
C              you do choose to use them, dimension them in your calling
C              program and in jac as arrays of appropriate length.
C
C
C
C  optionally replaceable norm routine:
C  dassl uses a weighted norm vnorm to measure the size
C  of vectors such as the estimated error in each step.
C  a function subprogram
C    double precision function vnorm(neq,v,wt,rpar,ipar)
C    dimension v(neq),wt(neq)
C  is used to define this norm.  here, v is the vector
C  whose norm is to be computed, and wt is a vector of
C  weights.  a vnorm routine has been included with dassl
C  which computes the weighted root-mean-square norm
C  given by
C    vnorm=sqrt((1/neq)*sum(v(i)/wt(i))**2)
C  this norm is suitable for most problems.  in some
C  special cases, it may be more convenient and/or
C  efficient to define your own norm by writing a function
C  subprogram to be called instead of vnorm.  this should
C  however, be attempted only after careful thought and
C  consideration.
C
C
C------output-after any return from dassl----
C
C  the principal aim of the code is to return a computed solution at
C  tout, although it is also possible to obtain intermediate results
C  along the way. to find out whether the code achieved its goal
C  or if the integration process was interrupted before the task was
C  completed, you must check the idid parameter.
C
C
C   T -- the solution was successfully advanced to the
C               output value of t.
C
C   Y(*) -- contains the computed solution approximation at t.
C
C   yprime(*) -- contains the computed derivative
C               approximation at t
C
C   idid -- reports what the code did
C
C                     *** task completed ***
C                reported by positive values of idid
C
C           idid = 1 -- a step was successfully taken in the
C                   intermediate-output mode. the code has not
C                   yet reached tout.
C
C           idid = 2 -- the integration to tout was successfully
C                   completed (t=tout) by stepping exactly to tout.
C
C           idid = 3 -- the integration to tout was successfully
C                   completed (t=tout) by stepping past tout.
C                   y(*) is obtained by interpolation.
C                   yprime(*) is obtained by interpolation.
C
C                    *** task interrupted ***
C                reported by negative values of idid
C
C           idid = -1 -- a large amount of work has been expended.
C                   (about 500 steps)
C
C           idid = -2 -- the error tolerances are too stringent.
C
C           idid = -3 -- the local error test cannot be satisfied
C                   because you specified a zero component in atol
C                   and the corresponding computed solution
C                   component is zero. thus, a pure relative error
C                   test is impossible for this component.
C
C           idid = -6 -- dassl had repeated error test
C                   failures on the last attempted step.
C
C           idid = -7 -- the corrector could not converge.
C
C           idid = -8 -- the matrix of partial derivatives
C                   is singular.
C
C           idid = -9 -- the corrector could not converge.
C                   there were repeated error test failures
C                   in this step.
C
C           idid =-10 -- the corrector could not converge
C                   because ires was equal to minus one.
C
C           idid =-11 -- ires equal to -2 was encountered
C                   and control is being returned to the
C                   calling program.
C
C           idid =-12 -- dassl failed to compute the initial
C                   yprime.
C
C
C
C           idid = -13,..,-32 -- not applicable for this code
C
C                    *** task terminated ***
C                reported by the value of idid=-33
C
C           idid = -33 -- the code has encountered trouble from which
C                   it cannot recover. a message is printed
C                   explaining the trouble and control is returned
C                   to the calling program. for example, this occurs
C                   when invalid input is detected.
C
C   rtol, atol -- these quantities remain unchanged except when
C               idid = -2. in this case, the error tolerances have been
C               increased by the code to values which are estimated
c               to be
C               appropriate for continuing the integration. however, the
C               reported solution at t was obtained using the input
c               values of rtol and atol.
C
C   rwork, iwork -- contain information which is usually of no
C               interest to the user but necessary for subsequent calls.
C               however, you may find use for
C
C               rwork(3)--which contains the step size h to be
C                       attempted on the next step.
C
C               rwork(4)--which contains the current value of the
C                       independent variable, i.e. the farthest point
C                       integration has reached. this will be different
C                       from t only when interpolation has been
C                       performed (idid=3).
C
C               rwork(7)--which contains the stepsize used
C                       on the last successful step.
C
C               iwork(7)--which contains the order of the method to
C                       be attempted on the next step.
C
C               iwork(8)--which contains the order of the method used
C                       on the last step.
C
C               iwork(11)--which contains the number of steps taken
c                          so far.
C
C               iwork(12)--which contains the number of calls to res
C                        so far.
C
C               iwork(13)--which contains the number of evaluations of
C                        the matrix of partial derivatives needed
c                        so far.
C
C               iwork(14)--which contains the total number
C                        of error test failures so far.
C
C               iwork(15)--which contains the total number
C                        of convergence test failures so far.
C                        (includes singular iteration matrix
C                        failures.)
C
C
C
C   input -- what to do to continue the integration
C            (calls after the first)                **
C
C     this code is organized so that subsequent calls to continue the
C     integration involve little (if any) additional effort on your
C     part. you must monitor the idid parameter in order to determine
C     what to do next.
C
C     recalling that the principal task of the code is to integrate
C     from t to tout (the interval mode), usually all you will need
C     to do is specify a new tout upon reaching the current tout.
C
C     do not alter any quantity not specifically permitted below,
C     in particular do not alter neq,t,y(*),yprime(*),rwork(*),iwork(*)
C     or the differential equation in subroutine res. any such
C     alteration constitutes a new problem and must be treated as such,
C     i.e. you must start afresh.
C
C     you cannot change from vector to scalar error control or vice
C     versa (info(2)) but you can change the size of the entries of
C     rtol, atol. increasing a tolerance makes the equation easier
C     to integrate. decreasing a tolerance will make the equation
C     harder to integrate and should generally be avoided.
C
C     you can switch from the intermediate-output mode to the
C     interval mode (info(3)) or vice versa at any time.
C
C     if it has been necessary to prevent the integration from going
C     past a point tstop (info(4), rwork(1)), keep in mind that the
C     code will not integrate to any tout beyound the currently
C     specified tstop. once tstop has been reached you must change
C     the value of tstop or set info(4)=0. you may change info(4)
C     or tstop at any time but you must supply the value of tstop in
C     rwork(1) whenever you set info(4)=1.
C
C     do not change info(5), info(6), iwork(1), or iwork(2)
C     unless you are going to restart the code.
C
C                    *** following a completed task ***
C     if
C     idid = 1, call the code again to continue the integration
C                  another step in the direction of tout.
C
C     idid = 2 or 3, define a new tout and call the code again.
C                  tout must be different from t. you cannot change
C                  the direction of integration without restarting.
C
C                    *** following an interrupted task ***
C                  to show the code that you realize the task was
C                  interrupted and that you want to continue, you
C                  must take appropriate action and set info(1) = 1
C     if
C     idid = -1, the code has taken about 500 steps.
C                  if you want to continue, set info(1) = 1 and
C                  call the code again. an additional 500 steps
C                  will be allowed.
C
C
C     idid = -2, the error tolerances rtol, atol have been
C                  increased to values the code estimates appropriate
C                  for continuing. you may want to change them
C                  yourself. if you are sure you want to continue
C                  with relaxed error tolerances, set info(1)=1 and
C                  call the code again.
C
C     idid = -3, a solution component is zero and you set the
C                  corresponding component of atol to zero. if you
C                  are sure you want to continue, you must first
C                  alter the error criterion to use positive values
C                  for those components of atol corresponding to zero
C                  solution components, then set info(1)=1 and call
C                  the code again.
C
C     idid = -4,-5  --- cannot occur with this code
C
C     idid = -6, repeated error test failures occurred on the
C                  last attempted step in dassl. a singularity in the
C                  solution may be present. if you are absolutely
C                  certain you want to continue, you should restart
C                  the integration.(provide initial values of y and
C                  yprime which are consistent)
C
C     idid = -7, repeated convergence test failures occurred
C                  on the last attempted step in dassl. an inaccurate or
C                  illconditioned jacobian may be the problem. if you
C                  are absolutely certain you want to continue, you
C                  should restart the integration.
C
C     idid = -8, the matrix of partial derivatives is singular.
C                  some of your equations may be redundant.
C                  dassl cannot solve the problem as stated.
C                  it is possible that the redundant equations
C                  could be removed, and then dassl could
C                  solve the problem. it is also possible
C                  that a solution to your problem either
C                  does not exist or is not unique.
C
C     idid = -9, dassl had multiple convergence test
C                  failures, preceeded by multiple error
C                  test failures, on the last attempted step.
C                  it is possible that your problem
C                  is ill-posed, and cannot be solved
C                  using this code.  or, there may be a
C                  discontinuity or a singularity in the
C                  solution.  if you are absolutely certain
C                  you want to continue, you should restart
C                  the integration.
C
C    idid =-10, dassl had multiple convergence test failures
C                  because ires was equal to minus one.
C                  if you are absolutely certain you want
C                  to continue, you should restart the
C                  integration.
C
C    idid =-11, ires=-2 was encountered, and control is being
C                  returned to the calling program.
C
C    idid =-12, dassl failed to compute the initial yprime.
C               this could happen because the initial
C               approximation to yprime was not very good, or
C               if a yprime consistent with the initial y
C               does not exist.  the problem could also be caused
C               by an inaccurate or singular iteration matrix.
C
C
C
C     idid = -13,..,-32 --- cannot occur with this code
C
C                       *** following a terminated task ***
C     if idid= -33, you cannot continue the solution of this
C                  problem. an attempt to do so will result in your
C                  run being terminated.
C
C-----------------------------------------------------------------------
C
C    ADDITIONAL DEBUGGING INFORMATION ADDED BY THIS USER
C-----------------------------------------------------------------------
C
C (note that in all cases setting info parameters to 0 will cause
c  program to run unaltered)
C
C
C   info(12) = 0 : no additional debugging information is printed out.
C              1 : Information about convergence failures and time step
c                  error control failures are printed out.
c                  Additionally, a subroutine called daserr is called
c                  which the user can
c                  tailor to his/her needs to get information about the
c                  failure.  Its form is;
c
c        subroutine DASERR(T,NEQ,H,E,NFUNCT,Y,YPRIME,WT,DELTA,IPAR,RPAR)
C
c                  where nfunct refers to the reason for the failure:
C                      Nfunct = 1 : Convergence failure
c                               2 : Time step error control failure
c                               3 : Successful completion of step in
c                                   sdaini
c                                   (Yprime is the calculated initial
c                                    deriv.)
c                               4 : Convergence failure, due to
c                                   nonnegativity
c                                   constraint being violated
c                               5 : Other failures
c                       T = Current time
c                       H = Current step size
c                       E = Normalized Error vector
c                            (If NFunct = 1, the update vector delta is
c                               returned instead)
c                       Y = Current solution vector
c                       Yprime = current derivative of the solution
c                                vector
c                       Wt = current weighting vector
C                       DELTA = the current value of delta x.
c
c             2 : Even more additional information than info(12) is
c                 printed out at every time step, concerning
c                 convergence rates and time step errors.
C
c   info(13) = 1 : Set newton damping factor in daini to 1
c
c   info(14) = 1 : Do full gauss elimination in daini for 5 steps before
c                  checking for convergence.  This step is useful for
c                  checking the consistency of an analytical jacobian.
C
c   info(15) .ne.0 : This is a linear problem!  Change the algorithm!
c             The ODE logic should be changed to reflect that iteration
c             with a current jacobian is not needed, as long as the time
c             step stays constant.  The type of linear problem
c             is specified by its value:
c                = 1 : linear non-autonomous jacobian
c                       (Jacobian explicitly depends on time)
c                = 2 : linear autonomous jacobian
c                       (Jacobian is not an explicit function of time)
C
C=======================================================================
C

C Dummy Variables
      EXTERNAL RES, JAC
      INTEGER  NEQ, INFO(15), IDID, LRW, LIW, IWORK(LIW), IPAR(*)
      DOUBLE PRECISION T, Y(NEQ), YPRIME(NEQ), RWORK(LRW), RTOL(NEQ),
     1                 ATOL(NEQ), TOUT, RPAR(*)

C Dassl Common Block
      INTEGER        NPD, NTEMP, LML, LMU, LMXORD, LMTYPE,
     *               LNST, LNRE, LNJE, LETF, LCTF, LIPVT
      COMMON/NWM001/ NPD, NTEMP, LML, LMU, LMXORD, LMTYPE,
     *               LNST, LNRE, LNJE, LETF, LCTF, LIPVT

C Local Variables
      LOGICAL DONE
      INTEGER LTSTOP, LHMAX,  LH,    LTN,    LCJ,  LCJOLD, LHOLD, LS,
     *        LROUND, LALPHA, LBETA, LGAMMA, LPSI, LSIGMA, LDELTA,
     *        LJCALC, LPHASE, LK, LKOLD, LNS, LNSTL, LNPD,
     *        LIWM, I, MXORD, LENPD, MBAND, LENIW, LENRW, NZFLG,
     *        LE, LPHI, LWT, LWM, LPD, ITEMP, MSAVE

      DOUBLE PRECISION HMAX, RTOLI, ATOLI, TN, UROUND, TDIST, HO,
     *                 YPNORM, RH, TSTOP, HMIN, X, H, R, TNEXT

      CHARACTER version*64
      DOUBLE PRECISION ZERO

C Externals
      DOUBLE PRECISION D1MACH, VNORM, SIGN77
      EXTERNAL         D1MACH, VNORM, SIGN77
C=======================================================================

      DATA LTSTOP, LHMAX,  LH,    LTN,    LCJ,  LCJOLD, LHOLD, LS,
     *     LROUND, LALPHA, LBETA, LGAMMA, LPSI, LSIGMA, LDELTA
     *   /      1,     2,   3,      4,      5,       6,     7,  8,
     *          9,    11,  17,     23,     29,      35,    41/
      DATA version /'@(#) ddassl.f 1.2 08/11/94'/
C=======================================================================
      ZERO = 0.0
      LML=1
      LMU=2
      LMXORD=3
      LMTYPE=4
      LJCALC=5
      LPHASE=6
      LK=7
      LKOLD=8
      LNS=9
      LNSTL=10
      LNST=11
      LNRE=12
      LNJE=13
      LETF=14
      LCTF=15
      LNPD = 16
      LIPVT=21
      LIWM=1
      IF(INFO(1).NE.0)GO TO 100
C
C-----------------------------------------------------------------------
C     this block is executed for the initial call only.
C     it contains checking of inputs and initializations.
C-----------------------------------------------------------------------
C
C     first check info array to make sure all elements of info
C     are either zero or one.
      DO 10 I=2,11
         IF(INFO(I).NE.0.AND.INFO(I).NE.1)GO TO 701
10       CONTINUE
C
      IF(NEQ.LE.0)GO TO 702
C
C     set pointers into iwork
C
C     check and compute maximum order
      MXORD=5
      IF(INFO(9).EQ.0)GO TO 20
         MXORD=IWORK(LMXORD)
         IF(MXORD.LT.1.OR.MXORD.GT.5)GO TO 703
20       IWORK(LMXORD)=MXORD
C
C     compute mtype,lenpd,lenrw.check ml and mu.
      IF(INFO(6).NE.0)GO TO 40
         LENPD=NEQ**2
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
         IF(INFO(5).NE.0)GO TO 30
            IWORK(LMTYPE)=2
            GO TO 60
30          IWORK(LMTYPE)=1
            GO TO 60
40    IF(IWORK(LML).LT.0.OR.IWORK(LML).GE.NEQ)GO TO 717
      IF(IWORK(LMU).LT.0.OR.IWORK(LMU).GE.NEQ)GO TO 718
      LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NEQ
      IF(INFO(5).NE.0)GO TO 50
         IWORK(LMTYPE)=5
         MBAND=IWORK(LML)+IWORK(LMU)+1
         MSAVE=(NEQ/MBAND)+1
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD+2*MSAVE
         GO TO 60
50       IWORK(LMTYPE)=4
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
C
C     check lengths of rwork and iwork
60    LENIW=20+NEQ
      IWORK(LNPD) = LENPD
      IF(LRW.LT.LENRW)GO TO 704
      IF(LIW.LT.LENIW)GO TO 705
C
C     check to see that tout is different from t
      IF(TOUT .EQ. T)GO TO 719
C
C     check hmax
      IF(INFO(7).EQ.0)GO TO 70
         HMAX=RWORK(LHMAX)
         IF(HMAX.LE.0.0D0)GO TO 710
70    CONTINUE
C
C     initialize counters
      IWORK(LNST)=0
      IWORK(LNRE)=0
      IWORK(LNJE)=0
C
      IWORK(LNSTL)=0
      IDID=1
      GO TO 200
C
C-----------------------------------------------------------------------
C     this block is for continuation calls
C     only. here we check info(1),and if the
C     last step was interrupted we check whether
C     appropriate action was taken.
C-----------------------------------------------------------------------
C
100   CONTINUE
      IF(INFO(1).EQ.1)GO TO 110
      IF(INFO(1).NE.-1)GO TO 701
C     if we are here, the last step was interrupted
C     by an error condition from dastep,and
C     appropriate action was not taken. this
C     is a fatal error.
      CALL XERRWV('DASSL--  THE LAST STEP TERMINATED WITH A NEGATIVE',
     *             49,201,0,0,0,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  VALUE (=I1) OF IDID AND NO APPROPRIATE',
     *             47,202,0,1,IDID,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  ACTION WAS TAKEN. RUN TERMINATED',
     *             41,203,1,0,0,0,0,ZERO,ZERO)
      RETURN
110   CONTINUE
      IWORK(LNSTL)=IWORK(LNST)
C
C-----------------------------------------------------------------------
C     this block is executed on all calls.
C     the error tolerance parameters are
C     checked, and the work array pointers
C     are set.
C-----------------------------------------------------------------------
C
200   CONTINUE
C     check rtol,atol
      NZFLG=0
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 210 I=1,NEQ
         IF(INFO(2).EQ.1)RTOLI=RTOL(I)
         IF(INFO(2).EQ.1)ATOLI=ATOL(I)
         IF(RTOLI.GT.0.0D0.OR.ATOLI.GT.0.0D0)NZFLG=1
         IF(RTOLI.LT.0.0D0)GO TO 706
         IF(ATOLI.LT.0.0D0)GO TO 707
210      CONTINUE
      IF(NZFLG.EQ.0)GO TO 708
C
C     set up rwork storage.iwork storage is fixed
C     in data statement.
      LE=LDELTA+NEQ
      LWT=LE+NEQ
      LPHI=LWT+NEQ
      LPD=LPHI+(IWORK(LMXORD)+1)*NEQ
      LWM=LPD
      NPD=1
      NTEMP = NPD + IWORK(LNPD)
      IF(INFO(1).EQ.1)GO TO 400
C
C-----------------------------------------------------------------------
C     this block is executed on the initial call
C     only. set the initial step size, and
C     the error weight vector, and phi.
C     compute initial yprime, if necessary.
C-----------------------------------------------------------------------
C
300   CONTINUE
      TN=T
      IDID=1
C
C     set error weight vector wt
      CALL WTSET(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      DO 305 I = 1,NEQ
         IF(RWORK(LWT+I-1) .GT. 0.0D0) THEN
         ELSE
           PRINT *,'I = ',I,'ATOL(I) = ',ATOL(I),'RTOL(I) = ',
     1     RTOL(I),'WT(I) = ',RWORK(LWT+I-1),'Y(I) = ',Y(I)
           GO TO 713
         END IF
305      CONTINUE
C
C     compute unit roundoff and hmin
      UROUND = D1MACH(4)
      RWORK(LROUND) = UROUND
      HMIN = 4.0D0*UROUND*DMAX1(DABS(T),DABS(TOUT))
C
C     check initial interval to see that it is long enough
      TDIST = DABS(TOUT - T)
      IF(TDIST .LT. HMIN) GO TO 714
C
C     check ho, if this was input
      IF (INFO(8) .EQ. 0) GO TO 310
         HO = RWORK(LH)
         IF ((TOUT - T)*HO .LT. 0.0D0) GO TO 711
         IF (HO .EQ. 0.0D0) GO TO 712
         GO TO 320
310    CONTINUE
C
C     compute initial stepsize, to be used by either
C     dastep or inider, depending on info(11)
      HO = 0.001D0*TDIST
      YPNORM = VNORM(NEQ,YPRIME,RWORK(LWT),RPAR,IPAR)
      IF (YPNORM .GT. 0.5D0/HO) HO = 0.5D0/YPNORM
      HO = SIGN77(HO,TOUT-T)
C     adjust ho if necessary to meet hmax bound
320   IF (INFO(7) .EQ. 0) GO TO 330
         RH = DABS(HO)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) HO = HO/RH
C     compute tstop, if applicable
330   IF (INFO(4) .EQ. 0) GO TO 340
         TSTOP = RWORK(LTSTOP)
         IF ((TSTOP - T)*HO .LT. 0.0D0) GO TO 715
         IF ((T + HO - TSTOP)*HO .GT. 0.0D0) HO = TSTOP - T
         IF ((TSTOP - TOUT)*HO .LT. 0.0D0) GO TO 709
C
C     compute initial derivative, if applicable
340   IF (INFO(11) .EQ. 0) GO TO 350
      CALL INIDER(T,Y,YPRIME,NEQ,
     *  RES,JAC,HO,RWORK(LWT),IDID,RPAR,IPAR,
     *  RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
     *  RWORK(LWM),IWORK(LIWM),HMIN,RWORK(LROUND),INFO(10)
     *  ,INFO(12),INFO(15),INFO(13),INFO(14))
      IF (IDID .LT. 0) THEN
        IF(INFO(12) .NE. 0) CALL DASERR(X,NEQ,HO,RWORK(LE)
     *                   ,5,Y,YPRIME,RWORK(LWT),
     *                    RWORK(LDELTA),IPAR,RPAR)
        GOTO 390
      ELSE
        IF(INFO(12) .NE. 0) CALL DASERR(X,NEQ,HO,RWORK(LE)
     *                   ,3,Y,YPRIME,RWORK(LWT),
     *                    RWORK(LDELTA),IPAR,RPAR)
      END IF
C
C     load h with ho.  store h in rwork(lh)
350   H = HO
      RWORK(LH) = H
C
C     load y and h*yprime into phi(*,1) and phi(*,2)
360   ITEMP = LPHI + NEQ
      DO 370 I = 1,NEQ
         RWORK(LPHI + I - 1) = Y(I)
370      RWORK(ITEMP + I - 1) = H*YPRIME(I)
C
390   GO TO 500
C
C-----------------------------------------------------------------------
C     this block is for continuation calls only. its
C     purpose is to check stop conditions before
C     taking a step.
C     adjust h if necessary to meet hmax bound
C-----------------------------------------------------------------------
C
400   CONTINUE
      UROUND = RWORK(LROUND)
      DONE = .FALSE.
      TN=RWORK(LTN)
      H=RWORK(LH)
      IF(INFO(7) .EQ. 0) GO TO 410
         RH = DABS(H)/RWORK(LHMAX)
         IF(RH .GT. 1.0D0) H = H/RH
410   CONTINUE
      IF(T .EQ. TOUT) GO TO 719
      IF((T - TOUT)*H .GT. 0.0D0) GO TO 711
      IF(INFO(4) .EQ. 1) GO TO 430
      IF(INFO(3) .EQ. 1) GO TO 420
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 490
      CALL INTRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
420   IF((TN-T)*H .LE. 0.0D0) GO TO 490
      IF((TN - TOUT)*H .GT. 0.0D0) GO TO 425
      CALL INTRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
425   CONTINUE
      CALL INTRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
430   IF(INFO(3) .EQ. 1) GO TO 440
      TSTOP=RWORK(LTSTOP)
      IF((TN-TSTOP)*H.GT.0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H.LT.0.0D0)GO TO 709
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 450
      CALL INTRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *   RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
440   TSTOP = RWORK(LTSTOP)
      IF((TN-TSTOP)*H .GT. 0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H .LT. 0.0D0) GO TO 709
      IF((TN-T)*H .LE. 0.0D0) GO TO 450
      IF((TN - TOUT)*H .GT. 0.0D0) GO TO 445
      CALL INTRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
445   CONTINUE
      CALL INTRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
450   CONTINUE
C     check whether we are with in roundoff of tstop
      IF(DABS(TN-TSTOP).GT.100.0D0*UROUND*
     *   (DABS(TN)+DABS(H)))GO TO 460
      IDID=2
      T=TSTOP
      DONE = .TRUE.
      GO TO 490
460   TNEXT=TN+H*(1.0D0+4.0D0*UROUND)
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 490
      H=(TSTOP-TN)*(1.0D0-4.0D0*UROUND)
      RWORK(LH)=H
C
490   IF (DONE) GO TO 590
C
C-----------------------------------------------------------------------
C     the next block contains the call to the
C     one-step integrator dastep.
C     this is a looping point for the integration
C     steps.
C     check for too many steps.
C     update wt.
C     check for too much accuracy requested.
C     compute minimum stepsize.
C-----------------------------------------------------------------------
C
500   CONTINUE
C     check for failure to compute initial yprime
      IF (IDID .EQ. -12) GO TO 527
C
C     check for too many steps
      IF((IWORK(LNST)-IWORK(LNSTL)).LT.500)
     *   GO TO 510
           IDID=-1
           GO TO 527
C
C     update wt
510   CALL WTSET(NEQ,INFO(2),RTOL,ATOL,RWORK(LPHI),
     *  RWORK(LWT),RPAR,IPAR)
C      print *,'i = ',1,'rtol(1) = ',rtol(1),
C     1           'atol(1) = ',atol(1),'wt(1) = '
C     1        ,rwork(1+lwt-1)
      DO 520 I=1,NEQ
         IF(RWORK(I+LWT-1).GT. 0.0D0) THEN
         ELSE
           IDID=-3
           PRINT *,'I = ',I,'RTOL(I) = ',RTOL(I)
     1          ,'ATOL(I) = '
     1      ,ATOL(I),'WT(I) = ',RWORK(I+LWT-1),'Y(I) = '
     1      ,RWORK(LPHI+I-1)
           GO TO 527
         END IF
520   CONTINUE
C
C     test for too much accuracy requested.
      R=VNORM(NEQ,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)*
     *   100.0D0*UROUND
      IF(R.LE.1.0D0)GO TO 525
C     multiply rtol and atol by r and return
      IF(INFO(2).EQ.1)GO TO 523
           RTOL(1)=R*RTOL(1)
           ATOL(1)=R*ATOL(1)
           IDID=-2
           GO TO 527
523   DO 524 I=1,NEQ
           RTOL(I)=R*RTOL(I)
524        ATOL(I)=R*ATOL(I)
      IDID=-2
      GO TO 527
525   CONTINUE
C
C     compute minimum stepsize
      HMIN=4.0D0*UROUND*DMAX1(DABS(TN),DABS(TOUT))
C
      CALL DASTEP(TN,Y,YPRIME,NEQ,
     *   RES,JAC,H,RWORK(LWT),INFO(1),IDID,RPAR,IPAR,
     *   RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
     *   RWORK(LWM),IWORK(LIWM),
     *   RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),
     *   RWORK(LPSI),RWORK(LSIGMA),
     *   RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),
     *   RWORK(LS),HMIN,RWORK(LROUND),
     *   IWORK(LPHASE),IWORK(LJCALC),IWORK(LK),
     *   IWORK(LKOLD),IWORK(LNS),INFO(10),
     *   INFO(12),INFO(15))
527   IF(IDID.LT.0)GO TO 600
C
C-----------------------------------------------------------------------
C     this block handles the case of a successful
C     return from dastep (idid=1) test for
C     stop conditions.
C-----------------------------------------------------------------------
C
      IF(INFO(4).NE.0)GO TO 540
           IF(INFO(3).NE.0)GO TO 530
             IF((TN-TOUT)*H.LT.0.0D0)GO TO 500
             CALL INTRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
530          IF((TN-TOUT)*H.GE.0.0D0)GO TO 535
             T=TN
             IDID=1
             GO TO 580
535          CALL INTRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
540   IF(INFO(3).NE.0)GO TO 550
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 542
         CALL INTRP(TN,TOUT,Y,YPRIME,NEQ,
     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
         T=TOUT
         IDID=3
         GO TO 580
542   IF(DABS(TN-TSTOP).LE.100.0D0*UROUND*
     *   (DABS(TN)+DABS(H)))GO TO 545
      TNEXT=TN+H*(1.0D0+4.0D0*UROUND)
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 500
      H=(TSTOP-TN)*(1.0D0-4.0D0*UROUND)
      GO TO 500
545   IDID=2
      T=TSTOP
      GO TO 580
550   IF((TN-TOUT)*H.GE.0.0D0)GO TO 555
      IF(DABS(TN-TSTOP).LE.100.0D0*UROUND*(DABS(TN)+DABS(H)))GO TO 552
      T=TN
      IDID=1
      GO TO 580
552   IDID=2
      T=TSTOP
      GO TO 580
555   CALL INTRP(TN,TOUT,Y,YPRIME,NEQ,
     *   IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID=3
580   CONTINUE
C
C--------------------------------------------------------
C     all successful returns from dassl are made from
C     this block.
C--------------------------------------------------------
C
590   CONTINUE
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C
C-----------------------------------------------------------------------
C     this block handles all unsuccessful
C     returns other than for illegal input.
C-----------------------------------------------------------------------
C
600   CONTINUE
      ITEMP=-IDID
      GO TO (610,620,630,690,690,640,650,660,670,675,
     *  680,685), ITEMP
C
C     the maximum number of steps was taken before
C     reaching tout
610   CALL XERRWV('DASSL--  AT CURRENT T (=R1)  500 STEPS',
     *             38,610,0,0,0,0,1,TN,ZERO)
      CALL XERRWV('DASSL--  TAKEN ON THIS CALL BEFORE REACHING TOUT',
     *             48,611,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     too much accuracy for machine precision
620   CALL XERRWV('DASSL--  AT T (=R1) TOO MUCH ACCURACY REQUESTED',
     *             47,620,0,0,0,0,1,TN,ZERO)
      CALL XERRWV('DASSL--  FOR PRECISION OF MACHINE. RTOL AND ATOL',
     *             48,621,0,0,0,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  WERE INCREASED TO APPROPRIATE VALUES',
     *             45,622,0,0,0,0,0,ZERO,ZERO)
C
      GO TO 690
C     wt(i) .le. 0.0d0 for some i (not at start of problem)
630   CALL XERRWV('DASSL--  AT T (=R1) SOME ELEMENT OF WT',
     *             38,630,0,0,0,0,1,TN,ZERO)
      CALL XERRWV('DASSL--  HAS BECOME .LE. 0.0',
     *             28,631,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     error test failed repeatedly or with h=hmin
640   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *             44,640,0,0,0,0,2,TN,H)
      CALL XERRWV(
     *'DASSL--  ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN',
     * 57,641,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     corrector convergence failed repeatedly or with h=hmin
650   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *             44,650,0,0,0,0,2,TN,H)
      CALL XERRWV('DASSL--  CORRECTOR FAILED TO CONVERGE REPEATEDLY',
     *             48,651,0,0,0,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  OR WITH ABS(H)=HMIN',
     *             28,652,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     the iteration matrix is singular
660   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *             44,660,0,0,0,0,2,TN,H)
      CALL XERRWV('DASSL--  ITERATION MATRIX IS SINGULAR',
     *             37,661,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     corrector failure preceeded by error test failures.
670   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *             44,670,0,0,0,0,2,TN,H)
      CALL XERRWV('DASSL--  CORRECTOR COULD NOT CONVERGE.  ALSO, THE',
     *             49,671,0,0,0,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  ERROR TEST FAILED REPEATEDLY.',
     *             38,672,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     corrector failure because ires = -1
675   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *             44,675,0,0,0,0,2,TN,H)
      CALL XERRWV('DASSL--  CORRECTOR COULD NOT CONVERGE BECAUSE',
     *             45,676,0,0,0,0,0,ZERO,ZERO)
      CALL XERRWV('DASSL--  IRES WAS EQUAL TO MINUS ONE',
     *             36,677,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     failure because ires = -2
680   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2)',
     *             40,680,0,0,0,0,2,TN,H)
      CALL XERRWV('DASSL--  IRES WAS EQUAL TO MINUS TWO',
     *             36,681,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
C
C     failed to compute initial yprime
685   CALL XERRWV('DASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE',
     *             44,685,0,0,0,0,2,TN,HO)
      CALL XERRWV('DASSL--  INITIAL YPRIME COULD NOT BE COMPUTED',
     *             45,686,0,0,0,0,0,ZERO,ZERO)
      GO TO 690
690   CONTINUE
      INFO(1)=-1
      T=TN
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C-----------------------------------------------------------------------
C     this block handles all error returns due
C     to illegal input, as detected before calling
C     dastep. first the error message routine is
C     called. if this happens twice in
C     succession, execution is terminated
C
C-----------------------------------------------------------------------
701   CALL XERRWV(
     *'DASSL--  SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE',
     * 55,1,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
702   CALL XERRWV('DASSL--  NEQ (=I1) .LE. 0',
     *             25,2,0,1,NEQ,0,0,ZERO,ZERO)
      GO TO 750
703   CALL XERRWV('DASSL--  MAXORD (=I1) NOT IN RANGE',
     *             34,3,0,1,MXORD,0,0,ZERO,ZERO)
      GO TO 750
704   CALL XERRWV(
     *'DASSL--  RWORK LENGTH NEEDED, LENRW (=I1), EXCEEDS LRW (=I2)',
     * 60,4,0,2,LENRW,LRW,0,ZERO,ZERO)
      GO TO 750
705   CALL XERRWV(
     *'DASSL--  IWORK LENGTH NEEDED, LENIW (=I1), EXCEEDS LIW (=I2)',
     * 60,5,0,2,LENIW,LIW,0,ZERO,ZERO)
      GO TO 750
706   CALL XERRWV('DASSL--  SOME ELEMENT OF RTOL IS .LT. 0',
     *             39,6,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
707   CALL XERRWV('DASSL--  SOME ELEMENT OF ATOL IS .LT. 0',
     *             39,7,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
708   CALL XERRWV('DASSL--  ALL ELEMENTS OF RTOL AND ATOL ARE ZERO',
     *             47,8,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
709   CALL XERRWV(
     *'DASSL--  INFO(4) = 1 AND TSTOP (=R1) BEHIND TOUT (=R2)',
     * 54,9,0,0,0,0,2,TSTOP,TOUT)
      GO TO 750
710   CALL XERRWV('DASSL--  HMAX (=R1) .LT. 0.0',
     *             28,10,0,0,0,0,1,HMAX,ZERO)
      GO TO 750
711   CALL XERRWV('DASSL--  TOUT (=R1) BEHIND T (=R2)',
     *             34,11,0,0,0,0,2,TOUT,T)
      GO TO 750
712   CALL XERRWV('DASSL--  INFO(8)=1 AND H0=0.0',
     *             29,12,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
713   CALL XERRWV('DASSL--  SOME ELEMENT OF WT IS .LE. 0.0',
     *             39,13,0,0,0,0,0,ZERO,ZERO)
      GO TO 750
714   CALL XERRWV(
     * 'DASSL--  TOUT (=R1) TOO CLOSE TO T (=R2) TO START INTEGRATION',
     * 61,14,0,0,0,0,2,TOUT,T)
      GO TO 750
715   CALL XERRWV('DASSL--  INFO(4)=1 AND TSTOP (=R1) BEHIND T (=R2)',
     *             49,15,0,0,0,0,2,TSTOP,T)
      GO TO 750
717   CALL XERRWV(
     *'DASSL--  ML (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ',
     * 52,17,0,1,IWORK(LML),0,0,ZERO,ZERO)
      GO TO 750
718   CALL XERRWV(
     *'DASSL--  MU (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ',
     * 52,18,0,1,IWORK(LMU),0,0,ZERO,ZERO)
      GO TO 750
719   CALL XERRWV('DASSL--  TOUT (=R1) IS EQUAL TO T (=R2)',
     *             39,19,0,0,0,0,2,TOUT,T)
      GO TO 750
750   IF(INFO(1).EQ.-1) GO TO 760
      INFO(1)=-1
      IDID=-33
      RETURN
760   CALL XERRWV('DASSL--  REPEATED OCCURRENCES OF ILLEGAL INPUT',
     *             46,801,0,0,0,0,0,ZERO,ZERO)
770   CALL XERRWV('DASSL--  RUN TERMINATED. APPARENT INFINITE LOOP',
     *             47,802,1,0,0,0,0,ZERO,ZERO)
      RETURN
C-----------end of subroutine dassl-------------------------------------
      END
      SUBROUTINE WTSET(NEQ, IWT, RTOL, ATOL, Y, WT, RPAR, IPAR)
C-----------------------------------------------------------------------
C     this subroutine sets the error weight vector
C     wt according to wt(i)=rtol(i)*abs(y(i))+atol(i),
C     i=1,-,n.
C     rtol and atol are scalars if iwt = 0,
C     and vectors if iwt = 1.
C-----------------------------------------------------------------------
C
C Dummy Variables
      INTEGER NEQ, IWT, IPAR(*)
      DOUBLE PRECISION RTOL(NEQ), ATOL(NEQ), Y(NEQ), WT(NEQ),
     1                 RPAR(*)

C Local Variables
      INTEGER I
      DOUBLE PRECISION RTOLI, ATOLI
C-----------------------------------------------------------------------
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      IF (IWT .EQ. 0) THEN
        DO 20 I=1,NEQ
           WT(I)=RTOLI*DABS(Y(I))+ATOLI
20         CONTINUE
      ELSE
        DO 30 I=1,NEQ
           WT(I)=RTOL(I)*DABS(Y(I))+ATOL(I)
30         CONTINUE
      END IF
      RETURN
C-----------end of subroutine wtset-------------------------------------
      END
      SUBROUTINE DASTEP (X, Y, YPRIME, NEQ,  RES, JAC, H, WT, JSTART,
     *                   IDID, RPAR, IPAR, PHI, DELTA, E, WM, IWM,
     *                   ALPHA, BETA, GAMMA, PSI, SIGMA, CJ, CJOLD,
     *                   HOLD, S, HMIN, UROUND, IPHASE, JCALC, K, KOLD,
     *                   NS, NONNEG, IDEBUG, LINEAR)
C
C-----------------------------------------------------------------------
C     dastep solves a system of differential/
C     algebraic equations of the form
C     g(x,y,yprime) = 0,  for one step (normally
C     from x to x+h).
C
C     the methods used are modified divided
C     difference,fixed leading coefficient
C     forms of backward differentiation
C     formulas. the code adjusts the stepsize
C     and order to control the local error per
C     step.
C
C
C     the parameters represent
C     x  --        independent variable
C     y  --        solution vector at x
C     yprime --    derivative of solution vector
C                  after successful step
C     neq --       number of equations to be integrated
C     res --       external user-supplied subroutine
C                  to evaluate the residual.  the call is
C                  call res(x,y,yprime,delta,ires,rpar,ipar)
C                  x,y,yprime are input.  delta is output.
C                  on input, ires=0.  res should alter ires only
C                  if it encounters an illegal value of y or a
C                  stop condition.  set ires=-1 if an input value
C                  of y is illegal, and dastep will try to solve
C                  the problem without getting ires = -1.  if
C                  ires=-2, dastep returns control to the calling
C                  program with idid = -11.
C     jac --       external user-supplied routine to evaluate
C                  the iteration matrix (this is optional)
C                  the call is of the form
C                  call jac(x,y,yprime,pd,cj,rpar,ipar)
C                  pd is the matrix of partial derivatives,
C                  pd=dg/dy+cj*dg/dyprime
C     h --         appropriate step size for next step.
C                  normally determined by the code
C     wt --        vector of weights for error criterion.
C     jstart --    integer variable set 0 for
C                  first step, 1 otherwise.
C     idid --      completion code with the following meanings%
C                  idid= 1 -- the step was completed successfully
C                  idid=-6 -- the error test failed repeatedly
C                  idid=-7 -- the corrector could not converge
C                  idid=-8 -- the iteration matrix is singular
C                  idid=-9 -- the corrector could not converge.
C                             there were repeated error test
C                             failures on this step.
C                  idid=-10-- the corrector could not converge
C                             because ires was equal to minus one
C                  idid=-11-- ires equal to -2 was encountered,
C                             and control is being returned to
C                             the calling program
C     RPAR,IPAR -- real and integer parameter arrays that
C                  are used for communication between the
C                  calling program and external user routines
C                  they are not altered by dastep
C     phi --       array of divided differences used by
C                  dastep. the length is neq*(k+1),where
C                  k is the maximum order
C     delta,e --   work vectors for dastep of length neq
C     wm,iwm --    real and integer arrays storing
C                  matrix information such as the matrix
C                  of partial derivatives,permutation
C                  vector,and various other information.
C
C     the other parameters are information
C     which is needed internally by dastep to
C     continue from step to step.
C
C----------------------------------------------------------------------
C         Aditional info
C            idebug = info(12)
C            Linear = info(15)
C
C-----------------------------------------------------------------------
C
C
C
C Dummy Variables
      EXTERNAL RES, JAC
      INTEGER NEQ, JSTART, IDID, IPHASE, JCALC, K, KOLD, IWM(*),
     *        NS, NONNEG, IDEBUG, LINEAR, IPAR(*)
      DOUBLE PRECISION X, Y(NEQ), YPRIME(NEQ), H, WT(NEQ), PHI(NEQ,*),
     *                 DELTA(NEQ), E(NEQ), WM(*), PSI(*), ALPHA(*),
     *                 BETA(*), GAMMA(*), SIGMA(*), RPAR(*), CJ, CJOLD,
     *                 HOLD, S, HMIN, UROUND

C Dassl Common Block
      INTEGER        NPD, NTEMP, LML, LMU, LMXORD, LMTYPE,
     *               LNST, LNRE, LNJE, LETF, LCTF, LIPVT
      COMMON/NWM001/ NPD, NTEMP, LML, LMU, LMXORD, LMTYPE,
     *               LNST, LNRE, LNJE, LETF, LCTF, LIPVT

C Local Variables
      LOGICAL CONVGD
      INTEGER MAXIT, NCF, NSF, NEF, KP1, KP2, KM1, NSP1, I, J, M, IRES,
     *        IER, KNEW, KDIFF, J1
      DOUBLE PRECISION XRATE, XOLD, DELNRM, TEMP1, TEMP2, ALPHAS,
     *                 ALPHA0, CJLAST, CK, ERK, TERK, ENORM, ERKM2,
     *                 ERKP1, ERKM1, TERKM1, TERKP1, PNORM, OLDNRM,
     *                 RATE, EST, TERKM2, HNEW, ERR, R

C Externals
      DOUBLE PRECISION VNORM
      EXTERNAL VNORM
C-----------------------------------------------------------------------
      DATA MAXIT/20/
      DATA XRATE/0.25D0/
C
C-----------------------------------------------------------------------
C     block 1.
C     initialize. on the first call,set
C     the order to 1 and initialize
C     other variables.
C-----------------------------------------------------------------------
C
C     initializations for all calls
      IDID=1
      XOLD=X
      NCF=0
      NSF=0
      NEF=0
      IF(JSTART .NE. 0) GO TO 120
C
C     if this is the first step,perform
C     other initializations
      IWM(LETF) = 0
      IWM(LCTF) = 0
      K=1
      KOLD=0
      HOLD=0.0D0
      JSTART=1
      PSI(1)=H
      CJOLD = 1.0D0/H
      CJ = CJOLD
      S = 100.D0
      JCALC = -1
      DELNRM=1.0D0
      IPHASE = 0
      NS=0
120   CONTINUE
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 2
C     compute coefficients of formulas for
C     this step.
C-----------------------------------------------------------------------
200   IF (IDEBUG .EQ. 2) WRITE(6,1201) X,H,K
205   CONTINUE
      KP1=K+1
      KP2=K+2
      KM1=K-1
      XOLD=X
      IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
      NS=MIN0(NS+1,KOLD+2)
      NSP1=NS+1
      IF(KP1 .LT. NS)GO TO 230
C
      BETA(1)=1.0D0
      ALPHA(1)=1.0D0
      TEMP1=H
      GAMMA(1)=0.0D0
      SIGMA(1)=1.0D0
      DO 210 I=2,KP1
         TEMP2=PSI(I-1)
         PSI(I-1)=TEMP1
         BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
         TEMP1=TEMP2+H
         ALPHA(I)=H/TEMP1
         SIGMA(I)=DFLOAT(I-1)*SIGMA(I-1)*ALPHA(I)
         GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
210      CONTINUE
      PSI(KP1)=TEMP1
230   CONTINUE
C
C     compute alphas, alpha0
      ALPHAS = 0.0D0
      ALPHA0 = 0.0D0
      DO 240 I = 1,K
        ALPHAS = ALPHAS - 1.0D0/DFLOAT(I)
        ALPHA0 = ALPHA0 - ALPHA(I)
240     CONTINUE
C
C     compute leading coefficient cj
      CJLAST = CJ
      CJ = -ALPHAS/H
C
C     compute variable stepsize error coefficient ck
      CK = DABS(ALPHA(KP1) + ALPHAS - ALPHA0)
      CK = DMAX1(CK,ALPHA(KP1))
C
C     decide whether new jacobian is needed
      TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
      TEMP2 = 1.0D0/TEMP1
      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
      IF (CJ .NE. CJLAST) S = 100.D0
C
C Is the matrix current for a linear non-autonomous problem?
      IF (LINEAR .EQ. 2) THEN
        IF (CJ .EQ. CJLAST) JCALC = 0
      END IF
C
C     change phi to phi star
      IF(KP1 .LT. NSP1) GO TO 280
      DO 270 J=NSP1,KP1
         DO 260 I=1,NEQ
260         PHI(I,J)=BETA(J)*PHI(I,J)
270      CONTINUE
280   CONTINUE
C
C     update time
      X=X+H
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 3
C     predict the solution and derivative,
C     and solve the corrector equation
C-----------------------------------------------------------------------
C
C     first,predict the solution and derivative
300   CONTINUE
      DO 310 I=1,NEQ
         Y(I)=PHI(I,1)
310      YPRIME(I)=0.0D0
      DO 330 J=2,KP1
         DO 320 I=1,NEQ
            Y(I)=Y(I)+PHI(I,J)
320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
330   CONTINUE
      PNORM = VNORM (NEQ,Y,WT,RPAR,IPAR)
C
C
C
C     solve the corrector equation using a
C     modified newton scheme.
      CONVGD= .TRUE.
      M=0
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
C
C
C     if indicated,reevaluate the
C     iteration matrix pd = dg/dy + cj*dg/dyprime
C     (where g(x,y,yprime)=0). set
C     jcalc to 0 as an indicator that
C     this has been done.
      IF(JCALC .NE. -1)GO TO 340
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      IF (IDEBUG .NE. 0 .AND. LINEAR .EQ. 0) THEN
        WRITE(6,1333) X , H
      ELSE IF (IDEBUG .EQ. 2) THEN
        WRITE(6,1334) X, K, KOLD, H, HOLD
      END IF
      CALL NJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,
     * IER,WT,E,WM,IWM,RES,IRES,UROUND,JAC,RPAR,IPAR)
      CJOLD=CJ
      S = 100.D0
      IF (IRES .LT. 0) GO TO 380
      IF(IER .NE. 0)GO TO 380
      NSF=0
C
C
C     initialize the error accumulation vector e.
340   CONTINUE
      DO 345 I=1,NEQ
345      E(I)=0.0D0
C
      S = 100.E0
C
C
C     corrector loop.
350   CONTINUE
C
C     multiply residual by temp1 to accelerate convergence
      TEMP1 = 2.0D0/(1.0D0 + CJ/CJOLD)
      DO 355 I = 1,NEQ
355     DELTA(I) = DELTA(I) * TEMP1
C
C     compute a new iterate (back-substitution).
C     store the correction in delta.
      CALL SOLVE(NEQ,DELTA,WM,IWM)
C
C     update y,e,and yprime
      DO 360 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
         E(I)=E(I)-DELTA(I)
360      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C
C     test for convergence of the iteration
      DELNRM=VNORM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (IDEBUG .EQ. 2) THEN
        WRITE(6,1362) M, DELNRM
      END IF
      IF (DELNRM .LE. 100.D0*UROUND*PNORM) GO TO 375
      IF (M .GT. 0) GO TO 365
         OLDNRM = DELNRM
         GO TO 367
365   RATE = (DELNRM/OLDNRM)**(1.0D0/DFLOAT(M))
      IF (RATE.GT.0.90D0) THEN
        IF (M .GT. 3 .OR.
     1      (DABS(DELNRM/OLDNRM) .GT. 4.0D0)) THEN
          GO TO 370
        ELSE
          IF (IDEBUG .EQ. 2) WRITE(6,1367) RATE
          GO TO 368
        END IF
      END IF
      S = RATE/(1.0D0 - RATE)
367   IF(LINEAR .GT. 0 .AND. JCALC .EQ. 0) GO TO 375
      IF (S*DELNRM .LE. 0.33D0) THEN
        IF (IDEBUG .EQ. 2) WRITE(6,1368) S*DELNRM
        GO TO 375
      ELSE
        IF (IDEBUG .EQ. 2) WRITE(6,1369) S*DELNRM
      END IF
C
C     the corrector has not yet converged.
C     update m and test whether the
C     maximum number of iterations have
C     been tried.
368   M=M+1
      IF(M.GE.MAXIT)GO TO 370
C
C     evaluate the residual
C     and go back to do another iteration
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL RES(X,Y,YPRIME,DELTA,IRES,
     *  RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
      GO TO 350
C
C
C     the corrector failed to converge in maxit
C     iterations. if the iteration matrix
C     is not current,re-do the step with
C     a new iteration matrix.
370   CONTINUE
      IF(JCALC.EQ.0)GO TO 380
      JCALC=-1
      GO TO 300
C
C
C     the iteration has converged.  if nonnegativity of solution is
C     required, set the solution nonnegative, if the perturbation
C     to do it is small enough.  if the change is too large, then
C     consider the corrector iteration to have failed.
375   IF(NONNEG .EQ. 0) GO TO 390
      DO 377 I = 1,NEQ
377      DELTA(I) = DMIN1(Y(I),0.0D0)
      DELNRM = VNORM(NEQ,DELTA,WT,RPAR,IPAR)
      IF(DELNRM .GT. 0.33D0) THEN
         IF (IDEBUG .NE. 0) THEN
           WRITE(6,1377) DELNRM
           CALL DASERR(X,NEQ,H,E,4,Y,YPRIME,WT,DELTA,IPAR,RPAR)
         END IF
         GO TO 380
      ELSE
        IF (IDEBUG .EQ. 2) THEN
          WRITE(6,1378) DELNRM
        END IF
      END IF
      DO 378 I = 1,NEQ
378      E(I) = E(I) - DELTA(I)
      GO TO 390
C
C
C     exits from block 3
C     no convergence with current iteration
C     matrix,or singular iteration matrix
380   IF (IDEBUG .NE. 0)
     *  CALL DASERR(X,NEQ,H,E,1,Y,YPRIME,WT,DELTA,IPAR,RPAR)
385   CONVGD= .FALSE.
390   JCALC = 1
      IF(.NOT.CONVGD)GO TO 600
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 4
C     estimate the errors at orders k,k-1,k-2
C     as if constant stepsize was used. estimate
C     the local error at order k and test
C     whether the current step is successful.
C-----------------------------------------------------------------------
C
C     estimate errors at orders k,k-1,k-2
      ENORM = VNORM(NEQ,E,WT,RPAR,IPAR)
      ERK = SIGMA(K+1)*ENORM
      TERK = FLOAT(K+1)*ERK
      EST = ERK
      KNEW=K
      IF(K .EQ. 1)GO TO 430
      DO 405 I = 1,NEQ
405     DELTA(I) = PHI(I,KP1) + E(I)
      ERKM1=SIGMA(K)*VNORM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM1 = FLOAT(K)*ERKM1
      IF(K .GT. 2)GO TO 410
      IF(TERKM1 .LE. 0.5*TERK)GO TO 420
      GO TO 430
410   CONTINUE
      DO 415 I = 1,NEQ
415     DELTA(I) = PHI(I,K) + DELTA(I)
      ERKM2=SIGMA(K-1)*VNORM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM2 = FLOAT(K-1)*ERKM2
      IF(DMAX1(TERKM1,TERKM2).GT.TERK)GO TO 430
C     lower the order
420   CONTINUE
      KNEW=K-1
      EST = ERKM1
C
C
C     calculate the local error for the current step
C     to see if the step was successful
430   CONTINUE
      ERR = CK * ENORM
      IF (IDEBUG .EQ. 2) WRITE(6,1430) ERR
      IF(ERR .GT. 1.0D0) THEN
        IF (IDEBUG .NE. 0)
     *    CALL DASERR(X,NEQ,H,E,2,Y,YPRIME,WT,DELTA,IPAR,RPAR)
        GO TO 600
      END IF
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 5
C     the step is successful. determine
C     the best order and stepsize for
C     the next step. update the differences
C     for the next step.
C-----------------------------------------------------------------------
      IDID=1
      IWM(LNST)=IWM(LNST)+1
      KDIFF=K-KOLD
      KOLD=K
      HOLD=H
C
C
C     estimate the error at order k+1 unless%
C        already decided to lower order, or
C        already using maximum order, or
C        stepsize not constant, or
C        order raised in previous step
      IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1
      IF(IPHASE .EQ. 0)GO TO 545
      IF(KNEW.EQ.KM1)GO TO 540
      IF(K.EQ.IWM(LMXORD)) GO TO 550
      IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550
      DO 510 I=1,NEQ
510      DELTA(I)=E(I)-PHI(I,KP2)
      ERKP1 = (1.0D0/DFLOAT(K+2))*VNORM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKP1 = FLOAT(K+2)*ERKP1
      IF(K.GT.1)GO TO 520
      IF(TERKP1.GE.0.5D0*TERK)GO TO 550
      GO TO 530
520   IF(TERKM1.LE.DMIN1(TERK,TERKP1))GO TO 540
      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550
C
C     raise order
530   K=KP1
      EST = ERKP1
      GO TO 550
C
C     lower order
540   K=KM1
      EST = ERKM1
      GO TO 550
C
C     if iphase = 0, increase order by one and multiply stepsize by
C     factor two
545   K = KP1
      HNEW = H*2.0D0
      H = HNEW
      GO TO 575
C
C
C     determine the appropriate stepsize for
C     the next step.
550   HNEW=H
      TEMP2=K+1
      R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      IF(R .LT. 2.0D0) GO TO 555
      HNEW = 2.0D0*H
      GO TO 560
555   IF(R .GT. 1.0D0) GO TO 560
      R = DMAX1(0.5D0,DMIN1(0.9D0,R))
      HNEW = H*R
560   H=HNEW
C
C
C     update differences for next step
575   CONTINUE
      IF(KOLD.EQ.IWM(LMXORD))GO TO 585
      DO 580 I=1,NEQ
580      PHI(I,KP2)=E(I)
585   CONTINUE
      DO 590 I=1,NEQ
590      PHI(I,KP1)=PHI(I,KP1)+E(I)
      DO 595 J1=2,KP1
         J=KP1-J1+1
         DO 595 I=1,NEQ
595      PHI(I,J)=PHI(I,J)+PHI(I,J+1)
      RETURN
C
C
C
C
C
C-----------------------------------------------------------------------
C     block 6
C     the step is unsuccessful. restore x,psi,phi
C     determine appropriate stepsize for
C     continuing the integration, or exit with
C     an error flag if there have been many
C     failures.
C-----------------------------------------------------------------------
600   IPHASE = 1
C
C     restore x,phi,psi
      X=XOLD
      IF(KP1.LT.NSP1)GO TO 630
      DO 620 J=NSP1,KP1
         TEMP1=1.0D0/BETA(J)
         DO 610 I=1,NEQ
610         PHI(I,J)=TEMP1*PHI(I,J)
620      CONTINUE
630   CONTINUE
      DO 640 I=2,KP1
640      PSI(I-1)=PSI(I)-H
C
C
C     test whether failure is due to corrector iteration
C     or error test
      IF(CONVGD)GO TO 660
      IWM(LCTF)=IWM(LCTF)+1
C
C
C     the newton iteration failed to converge with
C     a current iteration matrix.  determine the cause
C     of the failure and take appropriate action.
      IF(IER.EQ.0)GO TO 650
C
C     the iteration matrix is singular. reduce
C     the stepsize by a factor of 4. if
C     this happens three times in a row on
C     the same step, return with an error flag
      NSF=NSF+1
      IF (IDEBUG .NE. 0) WRITE(6,1450) X , H
      R = 0.25D0
      H=H*R
      IF (NSF .LT. 3 .AND. DABS(H) .GE. HMIN) GO TO 690
      IDID=-8
      GO TO 675
C
C
C     the newton iteration failed to converge for a reason
C     other than a singular iteration matrix.  if ires = -2, then
C     return.  otherwise, reduce the stepsize and try again, unless
C     too many failures have occured.
650   CONTINUE
      IF (IRES .GT. -2) GO TO 655
      IDID = -11
      GO TO 675
655   NCF = NCF + 1
      IF (IDEBUG .NE. 0) WRITE(6,1656) X , H , H*0.25
      R = 0.25D0
      H = H*R
      IF (NCF .LT. 15 .AND. DABS(H) .GE. HMIN) GO TO 690
      IDID = -7
      IF (IRES .LT. 0) IDID = -10
      IF (NEF .GE. 3) IDID = -9
      GO TO 675
C
C
C     the newton scheme converged,and the cause
C     of the failure was the error estimate
C     exceeding the tolerance.
660   NEF=NEF+1
      IWM(LETF)=IWM(LETF)+1
      IF (NEF .GT. 1) THEN
        GO TO 665
      ELSE
        IF (IDEBUG .NE. 0) WRITE(6,1662) X , H , K , KNEW
      END IF
C
C     on first error test failure, keep current order or lower
C     order by one.  compute new stepsize based on differences
C     of the solution.
      K = KNEW
      TEMP2 = K + 1
      R = 0.90D0*(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      R = DMAX1(0.25D0,DMIN1(0.9D0,R))
      H = H*R
      IF (DABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     on second error test failure, use the current order or
C     decrease order by one.  reduce the stepsize by a factor of
C     one quarter.
665   IF (NEF .GT. 2) THEN
        IF (IDEBUG .NE. 0) WRITE(6,1665) X , H , K , 1
        GO TO 670
      ELSE
        IF (IDEBUG .NE. 0) WRITE(6,1666) X , H , K , KNEW
      END IF
      K = KNEW
      H = 0.25D0*H
      IF (DABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     on third and subsequent error test failures, set the order to
C     one and reduce the stepsize by a factor of one quarter
670   K = 1
      H = 0.25D0*H
      IF (DABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C
C
C
C     for all crashes, restore y to its last value,
C     interpolate to find yprime at last x, and return
675   CONTINUE
      CALL INTRP(X,X,Y,YPRIME,NEQ,K,PHI,PSI)
      RETURN
C
C
C     go back and try this step again
690   IF (IDEBUG .NE. 0) WRITE(6,1691) X , H , K
      GO TO 205
C--------------------------------------------
C           Added format statements
1201  FORMAT(/10X,'DASTP: Old Time =',g11.4
     * ,' Attempt step with deltaT =',G11.4,' with K =',i3)
1333  FORMAT(10X,'DASTP:     New Jacobian Requested, New Time ='
     *,G14.4
     *   ,' deltaT =',G14.4)
1334  FORMAT(10X,'DASTP:     Linear problem needed a new Jacobian:'
     *,' New Time ='
     *      ,G14.4/20x,'K = ',I3,' KOLD = ',I3,' H = ',G14.4
     *      ,' HOLD = ',G14.4)
1362  FORMAT(10X,'DASTP:',15X,'Iteration # = ',I3
     *          ,' Norm of deltaY = ',G14.4)
1367  FORMAT(10X,'DASTP:',15X
     *  ,'Convergence not achieved, Solution blowing up, rate ='
     *  ,G14.4)
1368  FORMAT(10X,'DASTP:',15X
     *   ,'Successful convergence, Convergence error =',G14.4)
1369  FORMAT(10X,'DASTP:',15X
     *  ,'Convergence not achieved, Convergence error =',G14.4)
1377  FORMAT(10X,'DASTP: Nonnegativity constraint induced '
     *   ,'convergence failure, DELNRM = ',G11.4)
1378  FORMAT(10X,'DASTP:',15X,'Norm of negative part of solution = '
     *    ,G11.4)
1430  FORMAT(10x,'DASTP:',10X,'Time truncation error = ',g11.4)
1450  FORMAT(10X,'DASTP: Matrix is singular, XOLD = ',G11.4,' H = '
     *    ,G11.3)
1656  FORMAT(10X,'DASTP: Convergence test failure at X = ',g11.4
     *    ,' H = ',g11.4,', New stepsize = ',G11.4)
1662  FORMAT(10X,'DASTP: First time step error control failure, X = '
     *    ,g11.4,' Current H = ',G11.4,' K = ',I3,' New K = ',I3)
1665  FORMAT(10X,'DASTP: Third time step error control failure, X = '
     *    ,g11.4,' Current H = ',G11.4,' K = ',I3,' New K = ',I3)
1666  FORMAT(10X,'DASTP: Second time step error control failure, X = '
     *    ,g11.4,' Current H = ',G11.4,' K = ',I3,' New K = ',I3)
1691  FORMAT(/10X,'DASTP: Retry at old Time = ',g11.4
     * ,' with new deltaT = '
     *   ,g11.4,' and k = ',i3)
C------end of subroutine dastep------
      END
      SUBROUTINE SOLVE( NEQ, DELTA, WM, IWM)
C
C-----------------------------------------------------------------------
C     this routine manages the solution of the linear
C     system arising in the newton iteration.
C     matrices and real temporary storage and
C     real information are stored in the array wm.
C     integer matrix information is stored in
C     the array iwm.
C     for a dense matrix, the linpack routine
C     dgesl is called.
C     for a banded matrix,the linpack routine
C     dgbsl is called.
C-----------------------------------------------------------------------
C
C Dummy Variables
      INTEGER NEQ, IWM(*)
      DOUBLE PRECISION DELTA(NEQ),WM(*)

C Dassl Common Block
      INTEGER        NPD, NTEMP, LML, LMU, LMXORD, LMTYPE,
     *               LNST, LNRE, LNJE, LETF, LCTF, LIPVT
      COMMON/NWM001/ NPD, NTEMP, LML, LMU, LMXORD, LMTYPE,
     *               LNST, LNRE, LNJE, LETF, LCTF, LIPVT

C Local Variables
      INTEGER MEBAND, MTYPE

C*****xlinpk double precision
C      INTEGER MAXEQ, MAXBD, INFO
C      PARAMETER (MAXEQ=400, MAXBD=140)
C      DOUBLE PRECISION RELERR
C      DOUBLE PRECISION AA, DR, DC, XTEMP, R, RCOND, ANORM
C      COMMON /DEXTRA/ AA(MAXBD,MAXEQ),DR(MAXEQ),DC(MAXEQ),
C     1                XTEMP(MAXEQ), R(MAXEQ), RCOND, ANORM
C*****END xlinpk double precision
C-----------------------------------------------------------------------
      MTYPE=IWM(LMTYPE)
      GO TO(100,100,300,400,400),MTYPE
C
C     dense matrix
100   CALL DGESL(WM(NPD),NEQ,NEQ,IWM(LIPVT),DELTA,0)
      RETURN
C
C     dummy section for mtype=3
300   CONTINUE
      RETURN
C
C     banded matrix
400   MEBAND=2*IWM(LML)+IWM(LMU)+1
C*****no xlinpk
      CALL DGBSL(WM(NPD),MEBAND,NEQ,IWM(LML),
     *           IWM(LMU),IWM(LIPVT),DELTA,0)
C*****END no xlinpk
C*****xlinpk double precision
C      CALL XDGBSL(WM(NPD), MEBAND, AA, MAXBD, NEQ, IWM(LML), IWM(LMU),
C     *            IWM(LIPVT), DR, DC, DELTA, RCOND, XTEMP, RELERR,
C     *            INFO, R, ANORM)
C      IF (INFO .NE. 0) THEN
C        PRINT *,'INFO NOT EQUAL TO ZERO: results may contain '
C     1         ,'round-off error'
C      END IF
C*****END xlinpk double precision
      RETURN
C------end of subroutine solve------
      END
      SUBROUTINE NJAC(NEQ, X, Y, YPRIME, DELTA, CJ, H,
     *                IER, WT, E, WM, IWM, RES, IRES, UROUND, JAC,
     *                RPAR, IPAR)
C
C-----------------------------------------------------------------------
C     this routine computes the iteration matrix
C     pd=dg/dy+cj*dg/dyprime (where g(x,y,yprime)=0).
C     here pd is computed by the user-supplied
C     routine jac if iwm(mtype) is 1 or 4, and
C     it is computed by numerical finite differencing
C     if iwm(mtype)is 2 or 5
C     the parameters have the following meanings.
C     y        = array containing predicted values
C     yprime   = array containing predicted derivatives
C     delta    = residual evaluated at (x,y,yprime)
C                (used only if iwm(mtype)=2 or 5)
C     cj       = scalar parameter defining iteration matrix
C     h        = current stepsize in integration
C     ier      = variable which is .ne. 0
C                if iteration matrix is singular,
C                and 0 otherwise.
C     wt       = vector of weights for computing norms
C     e        = work space (temporary) of length neq
C     wm       = real work space for matrices. on
C                output it contains the lu decomposition
C                of the iteration matrix.
C     iwm      = integer work space containing
C                matrix information
C     res      = name of the external user-supplied routine
C                to evaluate the residual function g(x,y,yprime)
C     ires     = flag which is equal to zero if no illegal values
C                in res, and less than zero otherwise.  (if ires
C                is less than zero, the matrix was not completed)
C                in this case (if ires .lt. 0), then ier = 0.
C     uround   = the unit roundoff error of the machine being used.
C     jac      = name of the external user-supplied routine
C                to evaluate the iteration matrix (this routine
C                is only used if iwm(mtype) is 1 or 4)
C-----------------------------------------------------------------------
C
C Dummy Variables
      EXTERNAL RES, JAC, SIGN77
      INTEGER NEQ, IER, IWM(*), IRES, IPAR(*)
      DOUBLE PRECISION Y(NEQ), YPRIME(NEQ), DELTA(NEQ), WT(NEQ), E(NEQ),
     1                 WM(*), RPAR(*), X, CJ, H, UROUND, SIGN77

C Dassl Common Block
      INTEGER        NPD, NTEMP, LML, LMU, LMXORD, LMTYPE,
     *               LNST, LNRE, LNJE, LETF, LCTF, LIPVT
      COMMON/NWM001/ NPD, NTEMP, LML, LMU, LMXORD, LMTYPE,
     *               LNST, LNRE, LNJE, LETF, LCTF, LIPVT


C*****xlinpk double precision
C      INTEGER MAXEQ, MAXBD
C      PARAMETER (MAXEQ=400, MAXBD=140)
C      DOUBLE PRECISION AA, DR, DC, XTEMP, R, RCOND, ANORM
C      COMMON /DEXTRA/ AA(MAXBD,MAXEQ),DR(MAXEQ), DC(MAXEQ),
C     1                XTEMP(MAXEQ), R(MAXEQ), RCOND, ANORM
C*****END xlinpk double precision

C Local Variables
      INTEGER I, MTYPE, NPDM1, LENPD, NROW, L, MEBAND, MBAND, MBA,
     *        ISAVE, J, N, IPSAVE, K, I1, I2, II, MEB1, MSAVE
      DOUBLE PRECISION DEL, YSAVE, DELINV, SQUR, YPSAVE

C-----------------------------------------------------------------------
C
      IER = 0
      NPDM1=NPD-1
      MTYPE=IWM(LMTYPE)
      GO TO (100,200,300,400,500),MTYPE
C
C
C     dense user-supplied matrix
100   LENPD=NEQ*NEQ
      DO 110 I=1,LENPD
110      WM(NPDM1+I)=0.0D0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      GO TO 230
C
C
C     dense finite-difference-generated matrix
200   IRES=0
      NROW=NPDM1
      SQUR = DSQRT(UROUND)
      DO 210 I=1,NEQ
         DEL=SQUR*DMAX1(DABS(Y(I)),DABS(H*YPRIME(I)),
     *     DABS(WT(I)))
         DEL=SIGN77(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         Y(I)=Y(I)+DEL
         YPRIME(I)=YPRIME(I)+CJ*DEL
         CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0D0/DEL
         DO 220 L=1,NEQ
220      WM(NROW+L)=(E(L)-DELTA(L))*DELINV
      NROW=NROW+NEQ
      Y(I)=YSAVE
      YPRIME(I)=YPSAVE
210   CONTINUE
C
C
C     do dense-matrix lu decomposition on pd
230      CALL DGEFA(WM(NPD),NEQ,NEQ,IWM(LIPVT),IER)
      RETURN
C
C
C     dummy section for iwm(mtype)=3
300   RETURN
C
C
C     banded user-supplied matrix
400   LENPD=(2*IWM(LML)+IWM(LMU)+1)*NEQ
      DO 410 I=1,LENPD
410      WM(NPDM1+I)=0.0D0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 550
C
C
C     banded finite-difference-generated matrix
500   MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN0(MBAND,NEQ)
      MEBAND=MBAND+IWM(LML)
      MEB1=MEBAND-1
      MSAVE=(NEQ/MBAND)+1
      ISAVE=NTEMP-1
      IPSAVE=ISAVE+MSAVE
      IRES=0
      SQUR=SQRT(UROUND)
C--------------------------------------------------------------
      IRES = 1
C--------------------------------------------------------------
      DO 540 J=1,MBA
         DO 510 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          WM(ISAVE+K)=Y(N)
          WM(IPSAVE+K)=YPRIME(N)
          DEL=SQUR*DMAX1(DABS(Y(N)),DABS(H*YPRIME(N)),
     *      DABS(WT(N)))
          DEL=SIGN77(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          Y(N)=Y(N)+DEL
510       YPRIME(N)=YPRIME(N)+CJ*DEL
      CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) RETURN
      DO 530 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*DMAX1(DABS(Y(N)),DABS(H*YPRIME(N)),
     *      DABS(WT(N)))
          DEL=SIGN77(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0D0/DEL
          I1=MAX0(1,(N-IWM(LMU)))
          I2=MIN0(NEQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)+NPDM1
          DO 520 I=I1,I2
520         WM(II+I)=(E(I)-DELTA(I))*DELINV
530      CONTINUE
540   CONTINUE
C-----------------------------------------------------------
      IRES = 0
C----------------------------------------------------------
C
C
C     do lu decomposition of banded pd
550   CONTINUE
C*****no xlinpk
      CALL DGBFA(WM(NPD),MEBAND,NEQ,
     *    IWM(LML),IWM(LMU),IWM(LIPVT),IER)
C*****END no xlinpk
C*****xlinpk double precision
C      IF (NEQ .GT. MAXEQ) THEN
C        PRINT *,'DDASSL NJAC: SIZE OF PROBLEM, ',NEQ
C     1        ,' HAD EXCEEDED ALLOCATED SPACE'
C     1        ,' IN DEXTRA, ', MAXEQ
C        STOP
C      END IF
C      IF (MEBAND .GT. MAXBD) THEN
C        PRINT *,'DDASSL NJAC: BAND SIZE OF PROBLEM, ',MEBAND
C     1        ,' HAD EXCEEDED ALLOCATED SPACE'
C     1        ,' IN DEXTRA, ', MAXBD
C        STOP
C      END IF
C      CALL XDGBCO(WM(NPD), MEBAND, AA, MAXBD, NEQ, IWM(LML),
C     1            IWM(LMU), IWM(LIPVT),
C     1            RCOND, DR, DC, XTEMP, ANORM)
C*****END xlinpk double precision
C
      RETURN
C------end of subroutine njac------
      END
      SUBROUTINE INTRP(X,XOUT,YOUT,YPOUT,NEQ,KOLD,PHI,PSI)
C
C-----------------------------------------------------------------------
C     the methods in subroutine dastep use polynomials
C     to approximate the solution. intrp approximates the
C     solution and its derivative at time xout by evaluating
C     one of these polynomials,and its derivative,there.
C     information defining this polynomial is passed from
C     dastep, so intrp cannot be used alone.
C
C     the parameters are%
C     x     the current time in the integration.
C     xout  the time at which the solution is desired
C     yout  the interpolated approximation to y at xout
C           (this is output)
C     ypout the interpolated approximation to yprime at xout
C           (this is output)
C     neq   number of equations
C     kold  order used on last successful step
C     phi   array of scaled divided differences of y
C     psi   array of past stepsize history
C-----------------------------------------------------------------------
C
C Dummy Variables
      INTEGER NEQ, KOLD
      DOUBLE PRECISION YOUT(NEQ), YPOUT(NEQ), PHI(NEQ,*), PSI(*),
     1                 X, XOUT

C Local Variables
      INTEGER KOLDP1, I, J
      DOUBLE PRECISION TEMP1, C, D, GAMMA
C-----------------------------------------------------------------------
      KOLDP1=KOLD+1
      TEMP1=XOUT-X
      DO 10 I=1,NEQ
         YOUT(I)=PHI(I,1)
10       YPOUT(I)=0.0D0
      C=1.0D0
      D=0.0D0
      GAMMA=TEMP1/PSI(1)
      DO 30 J=2,KOLDP1
         D=D*GAMMA+C/PSI(J-1)
         C=C*GAMMA
         GAMMA=(TEMP1+PSI(J-1))/PSI(J)
         DO 20 I=1,NEQ
            YOUT(I)=YOUT(I)+C*PHI(I,J)
20          YPOUT(I)=YPOUT(I)+D*PHI(I,J)
30       CONTINUE
      RETURN
C
C------end of subroutine intrp------
      END
      SUBROUTINE INIDER(X,Y,YPRIME,NEQ,
     *   RES,JAC,H,WT,IDID,RPAR,IPAR,
     *   PHI,DELTA,E,WM,IWM,
     *   HMIN,UROUND,NONNEG,
     *   IDEBUG, LINEAR , IDAMP , IFULLG)
C
C
C-------------------------------------------------------
C     inider takes one step of size h or smaller
C     with the backward euler method, to
C     find yprime at the initial time x. a modified
C     damped newton iteration is used to
C     solve the corrector iteration.
C
C     the initial guess yprime is used in the
C     prediction, and in forming the iteration
C     matrix, but is not involved in the
C     error test. this may have trouble
C     converging if the initial guess is no
C     good, or if g(xy,yprime) depends
C     nonlinearly on yprime.
C
C     the parameters represent:
C     x --         independent variable
C     y --         solution vector at x
C     yprime --    derivative of solution vector
C     neq --       number of equations
C     h --         stepsize. imder may use a stepsize
C                  smaller than h.
C     wt --        vector of weights for error
C                  criterion
C     idid --      completion code with the following meanings
C                  idid= 1 -- yprime was found successfully
C                  idid=-12 -- inider failed to find yprime
C     rpar,ipar -- real and integer parameter arrays
C                  that are not altered by inider
C     phi --       work space for inider
C     delta,e --   work space for inider
C     wm,iwm --    real and integer arrays storing
C                  matrix information
C
C-----------------------------------------------------------------
C        ADDITIONAL PARAMETERS
C
C     IDEBUG = info(12)
C     LINEAR = info(15)
C     IDAMP  = info(13)
C     IFULLG = info(14)
C------------------------------------------------------------------
C
C Dummy Variables
      EXTERNAL RES,JAC
      INTEGER NEQ, IDID, IWM(*), IPAR(*), IDEBUG, LINEAR, IDAMP,
     *        IFULLG, NONNEG
      DOUBLE PRECISION X, Y(NEQ), YPRIME(NEQ), H, WT(NEQ), PHI(NEQ,*),
     *                 DELTA(NEQ), E(NEQ), WM(*), HMIN, UROUND,
     *                 RPAR(*)

C Dassl Common Block
      INTEGER        NPD, NTEMP, LML, LMU, LMXORD, LMTYPE,
     *               LNST, LNRE, LNJE, LETF, LCTF, LIPVT
      COMMON/NWM001/ NPD, NTEMP, LML, LMU, LMXORD, LMTYPE,
     *               LNST, LNRE, LNJE, LETF, LCTF, LIPVT

C Local Variables
      LOGICAL CONVGD
      INTEGER MAXIT, MJAC, NEF, NCF, NSF, I, IRES, JCALC, M, IER
      DOUBLE PRECISION DAMP, YNORM, CJ, XNEW, S, DELNRM, OLDNRM,
     *                 RATE, ERR, R

C Externals
      DOUBLE PRECISION VNORM
      EXTERNAL VNORM
C-----------------------------------------------------------------------
      DATA MAXIT/15/,MJAC/5/
      DATA DAMP/0.75D0/
C
C---------------------------------------------------
C     block 1.
C     initializations.
C---------------------------------------------------
C
      IDID=1
      NEF=0
      NCF=0
      NSF=0
      YNORM=VNORM(NEQ,Y,WT,RPAR,IPAR)
      IF (IDAMP .EQ. 1 .OR. IFULLG .EQ. 1) DAMP = 1.0
      IF (IFULLG .EQ. 1) MJAC = 1
C
C     save y and yprime in phi
      DO 100 I=1,NEQ
         PHI(I,1)=Y(I)
100      PHI(I,2)=YPRIME(I)

C
C
C----------------------------------------------------
C     block 2.
C     do one backward euler step.
C----------------------------------------------------
C
C     set up for start of corrector iteration
200   IF(IDEBUG .EQ. 2) WRITE(6,1201)X,H,1
205   CJ=1.0D0/H
      XNEW=X+H
C
C     predict solution and derivative

      DO 250 I=1,NEQ
250     Y(I)=Y(I)+H*YPRIME(I)
C
      JCALC=-1
      M=0
      CONVGD=.TRUE.
C
C
C     corrector loop.
300   IWM(LNRE)=IWM(LNRE)+1
      IRES=0

      CALL RES(XNEW,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES.LT.0) GO TO 430
C
C
C     evaluate the iteration matrix
      IF (JCALC.NE.-1) GO TO 310
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      IF (IDEBUG .EQ. 2) WRITE(6,1301)
      CALL NJAC(NEQ,XNEW,Y,YPRIME,DELTA,CJ,H,
     *   IER,WT,E,WM,IWM,RES,IRES,
     *   UROUND,JAC,RPAR,IPAR)

      S=100000.D0
      IF (IRES.LT.0) GO TO 430
      IF (IER.NE.0) GO TO 430
      NSF=0

C
C
C
C     multiply residual by damping factor
310   CONTINUE
      DO 320 I=1,NEQ
320      DELTA(I)=DELTA(I)*DAMP

C
C     compute a new iterate (back substitution)
C     store the correction in delta

      CALL SOLVE(NEQ,DELTA,WM,IWM)

C
C     update y and yprime

      DO 330 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
330      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)

C
C     test for convergence of the iteration.

      DELNRM=VNORM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.LE.100.D0*UROUND*YNORM) THEN
        IF (IDEBUG .NE. 0) WRITE(6,1331) DELNRM
        GO TO 400
      END IF

      IF (M.GT.0) GO TO 340
         OLDNRM=DELNRM
         GO TO 350

340   RATE=(DELNRM/OLDNRM)**(1.0D0/DFLOAT(M))
      IF (RATE.GT.0.90D0) THEN
        IF (M .GT. 3) THEN
          GO TO 430
        ELSE
          GO TO 355
        END IF
      END IF
      S=RATE/(1.0D0-RATE)

350   IF (S*DELNRM .LE. 0.33D0) GO TO 400
C
C
C     the corrector has not yet converged. update
C     m and and test whether the maximum
C     number of iterations have been tried.
C     every mjac iterations, get a new
C     iteration matrix.

355   M=M+1
      IF (M.GE.MAXIT) GO TO 430

      IF ((M/MJAC)*MJAC.EQ.M) JCALC=-1

      GO TO 300

C
C
C     the iteration has converged.
C     check nonnegativity constraints
400   IF (NONNEG.EQ.0) GO TO 450
      DO 410 I=1,NEQ
410      DELTA(I)=DMIN1(Y(I),0.0D0)

      DELNRM=VNORM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM .GT. 0.33D0) THEN
        IF (IDEBUG .NE. 0) THEN
          WRITE(6,1411) DELNRM
          CALL DASERR(X,NEQ,H,E,4,Y,YPRIME,WT,DELTA,IPAR,RPAR)
        END IF
        GO TO 435
      END IF

      DO 420 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
420      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
      GO TO 450
C
C
C     exits from corrector loop.
430   IF (IDEBUG .NE. 0)
     *    CALL DASERR(X,NEQ,H,DELTA,1,Y,YPRIME,WT,DELTA,IPAR,RPAR)
435   CONVGD=.FALSE.
450   IF (.NOT.CONVGD) GO TO 600
C
C
C
C-----------------------------------------------------
C     block 3.
C     the corrector iteration converged.
C     do error test.
C-----------------------------------------------------
C
      DO 510 I=1,NEQ
510      E(I)=Y(I)-PHI(I,1)
      ERR=VNORM(NEQ,E,WT,RPAR,IPAR)
      IF (ERR.LE.1.0D0) THEN
        IF (IDEBUG .EQ. 2) WRITE(6,1511) ERR
        RETURN
      ELSE
        IF (IDEBUG .NE. 0) THEN
          WRITE(6,1512) ERR
          CALL DASERR(X,NEQ,H,E,2,Y,YPRIME,WT,DELTA,IPAR,RPAR)
        END IF
      END IF
C
C
C
C--------------------------------------------------------
C     block 4.
C     the backward euler step failed. restore y
C     and yprime to their original values.
C     reduce stepsize and try again, if
C     possible.
C---------------------------------------------------------
C

600   CONTINUE
      DO 610 I=1,NEQ
         Y(I)=PHI(I,1)
610      YPRIME(I)=PHI(I,2)

      IF (CONVGD) GO TO 640
      IF (IER.EQ.0) GO TO 620
         NSF=NSF+1
         H=H*0.25D0
         IF (NSF.LT.3.AND.DABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
620   IF (IRES.GT.-2) GO TO 630
         IDID=-12
         RETURN
630   NCF=NCF+1
      H=H*0.25D0
      IF (NCF.LT.10.AND.DABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN

640   NEF=NEF+1
      R=0.90D0/(2.0D0*ERR+0.0001D0)
      R=DMAX1(0.1D0,DMIN1(0.5D0,R))
      H=H*R
      IF (DABS(H).GE.HMIN.AND.NEF.LT.10) GO TO 690
         IDID=-12
         RETURN
690   IF (IDEBUG .NE. 0) WRITE(6,1691)X, H
      GO TO 205
C-----------------------------------------------
C        ADDITIONAL FORMAT STATEMENTS
1201  FORMAT(/10X,'INIDER: Old Time =',G11.4
     * ,' Attempt step with deltaT =',G11.4,' with K =',I3)
1301  FORMAT(10X,'INIDER: Get a new Jacobian')
1331  FORMAT(10X,'INIDER: Round-off condition is triggered'
     *    /20X,'Consider step to be converged, DELNRM = ',G11.4)
1411  FORMAT(10X,'INIDER: Nonnegativity constraint induced '
     *    ,'an error, Norm of negative part of solution = '
     *    ,G11.4)
1511  FORMAT(10X,'INIDER:',10X,'Time truncation error =',G11.4)
1512  FORMAT(10X,'INIDER: Time Step truncation error '
     *  ,'failure, Error ='
     1          ,G11.4)
1691  FORMAT(/10X,'INIDER: Retry at Old Time =',g11.4
     *,'Attempt to solve problem with a'
     *       ,' new step size, H = ',g11.4)
C-------------end of subroutine inider----------------------
      END
      DOUBLE PRECISION FUNCTION VNORM(NEQ, V, WT, RPAR, IPAR)
C
C-----------------------------------------------------------------------
C     this function routine computes the weighted
C     root-mean-square norm of the vector of length
C     neq contained in the array v,with weights
C     contained in the array wt of length neq.
C        vnorm=sqrt((1/neq)*sum(v(i)/wt(i))**2)
C-----------------------------------------------------------------------
C
C Dummy Variables
      INTEGER NEQ, IPAR(*)
      DOUBLE PRECISION V(NEQ), WT(NEQ), RPAR(*)

C Local Variables
      INTEGER I
      DOUBLE PRECISION VMAX, SUM
C-----------------------------------------------------------------------
      VNORM = 0.0D0
      VMAX = 0.0D0
      DO 10 I = 1,NEQ
      if (wt(i) .le. 0.0d0 ) then
        print *,'VNORM: i = ',i,' wt = ',wt(i)
        stop
      end if
10      IF(DABS(V(I)/WT(I)) .GT. VMAX) VMAX = DABS(V(I)/WT(I))
      IF(VMAX .LE. 0.0D0) GO TO 30
      SUM = 0.0D0
      DO 20 I = 1,NEQ
20      SUM = SUM + ((V(I)/WT(I))/VMAX)**2
      VNORM = VMAX*DSQRT(SUM/DFLOAT(NEQ))
30    CONTINUE
      RETURN
C------end of function vnorm------
      END

