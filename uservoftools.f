c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             FUNC2D1                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c X,Y      = coordinates of the point where VALUE is computed         c
c On return:                                                          c
c===========                                                          c
c FUNC2D1  = VALUE of the implicit interface shape function:          c
c            > 0 (inside the interface), < 0 (outside the             c
c            interface); = 0 (on the interface)                       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      FUNCTION FUNC2D1(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION FUNC2D1,X,Y
c* Circle with radius 1.0 centered at (2.0,2.0):      
c* Exact area encolsed by the interface=PI      
      FUNC2D1=((X-2.0)**2+(Y-2.0)**2-1.0**2)
      RETURN
      END 
c------------------------- END OF FUNC2D1 ----------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             FUNC2D2                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c X,Y      = coordinates of the point where VALUE is computed         c
c On return:                                                          c
c===========                                                          c
c FUNC2D2  = VALUE of the implicit interface shape function:          c
c            > 0 (inside the interface), < 0 (outside the             c
c            interface); = 0 (on the interface)                       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      FUNCTION FUNC2D2(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION FUNC2D2,X,Y
c* Egg shape: Ellipse with semi-major axis 1.25, semi-minor axis 0.75 and
c* centered at (0.5,0.5) but with parameter on major axis
    
      FUNC2D2=((X-0.5)**2+(Y-0.75)**2-0.15**2)  
C	FUNC2D2=(-1.5*X - Y + 5.1)

      RETURN
      END 
c------------------------- END OF FUNC2D2 ----------------------------c
c---------------------------------------------------------------------c
      
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             FUNC2D3                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c X,Y    = coordinates of the point where VALUE is computed         c
c On return:                                                          c
c===========                                                          c
c FUNC2D3  = VALUE of the implicit interface shape function:          c
c            > 0 (inside the interface), < 0 (outside the             c
c            interface); = 0 (on the interface)                       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      FUNCTION FUNC2D3(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION FUNC2D3,X,Y
c* straight line   
      FUNC2D3=(X - Y + 0.45)
      RETURN
      END 
c------------------------- END OF FUNC2D3 ----------------------------c
c---------------------------------------------------------------------c

c                             FUNC3D1                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c X,Y,Z    = coordinates of the point where VALUE is computed         c
c On return:                                                          c
c===========                                                          c
c FUNC3D1  = VALUE of the implicit interface shape function:          c
c            > 0 (inside the interface), < 0 (outside the             c
c            interface); = 0 (on the interface)                       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      FUNCTION FUNC3D1(X,Y,Z)
      IMPLICIT NONE
      DOUBLE PRECISION FUNC3D1,X,Y,Z
c* Sphere with radious 0.6 centered at (0.5,0.5,0.5):      
c* Exact volume encolsed by the interface=(4/3)*PI*0.6**3      
      FUNC3D1=-1.0*((X-0.5)**2+(Y-0.5)**2+(Z-0.5)**2-0.6**2)
      RETURN
      END 
c------------------------- END OF FUNC3D1 ----------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             FUNC3D2                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c X,Y,Z    = coordinates of the point where VALUE is computed         c
c On return:                                                          c
c===========                                                          c
c FUNC3D2  = VALUE of the implicit interface shape function:          c
c            > 0 (inside the interface), < 0 (outside the             c
c            interface); = 0 (on the interface)                       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      FUNCTION FUNC3D2(X,Y,Z)
      IMPLICIT NONE
      DOUBLE PRECISION FUNC3D2,X,Y,Z
c* Torus with major radius 2/3, minor radius 1/3 and centered at (0.5,0.5,0.5):
c* Exact volume encolsed by the interface=2*PI**2*(2/3)*(1/3)**2      
      FUNC3D2=(1.0/3.0)**2-((2.0/3.0)-((X-0.5)**2+(Z-0.5)**2)**0.5)**2-
     -     (Y-0.5)**2
      RETURN
      END 
c------------------------- END OF FUNC3D2 ----------------------------c
c---------------------------------------------------------------------c
