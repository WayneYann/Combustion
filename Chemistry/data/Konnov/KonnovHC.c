



#if 0




\\
\\
\\  This is the mechanism file
\\
\\
ELEM H C O N AR
END
SPECIES
H         H2        O         O2        OH        HO2       H2O
H2O2      CO        CO2       HCO       CH3       CH4       C2H6
CH2O      C2H5      CH2       CH3O      CH2OH     CH        C2H2
C2H4      C2H3      CH3OH     CH3HCO    C2H       CH2CO     HCCO
C2H4O     ! ethylene oxide
SCH2      ! CH2(1) - singlet
C2        C2O       CH3CO     C  
CH3CO3    ! CH3C(O)OO    
CH3CO3H   ! CH3C(O)OOH 
CH3O2     CH3O2H    C2H5O2H   C2H5O2    CH3CO2    
CH3CO2H   ! CH3C(O)OH acetic acid
C2H5OH
C2H5O     ! CH3-CH2-O(.)
SC2H5O    ! CH3-CH(.)-OH
PC2H5O    ! CH2(.)-CH2-OH
CH2HCO    CN        H2CN      N         NH        HCN       NO 
HCNO      HOCN      HNCO      NCO       N2O       NH2       N2O3
HNO       NO2       C2N2      NNH       NH3       N2H2      HONO 
NO3       HNO3      N2H3      N2H4      CNN       HCNN      N2O4 
NH2OH     HNOH      H2NO      HNNO      HCNH
H2CNO     ! H2C*N=O 
CH3NO     ! Nitrosomethyl
CH2CHOW   ! CH2CHO exited
C2H3O     ! CH2-CH(.)
          ! |__O__|
CH3HCOW   ! CH3HCO exited

C3H6OH    O2C3H6OH
C3H5O2    ! CH2=CHCH2OO
C3H5O2H   ! CH2=CHCH2OOH
C3H5O     ! CH2=CHCH2O(.)
NC3H7O2   ! CH3CH2CH2OO(.)
NC3H7O2H  ! CH3CH2CH2OOH
IC3H7O2   ! CH3CH(OO.)CH3
IC3H7O2H  ! CH3CH(OOH)CH3
IC3H7O    ! CH3CH(O)CH3
NC3H7O    ! CH3CH2CH2O(.)
C3H6      C3H8      
IC3H7     ! CH3C(.)HCH3
NC3H7     ! CH3CH2C(.)H2
C3H2    
C3H3      ! propargyl CH2=C=CH
SC3H5     ! 1-propenyl or 2-methylvinyl or H3C-CH=(.)CH
PC3H4     ! propyne or methylacetylene CH3CCH
TC3H5     ! 2-propenyl (CH2=C)(.)-CH3 or 1-methylvinyl radical
C3H6O     ! propene oxide
C2H5CHO   ! propionaldehyde
C2H5CO    ! propionaldehyde radical
C3H5      ! allyl radical CH2CH=CH2
C3H4      ! allene or propadiene CH2=C=CH2

IC4H7     ! CH2=CHCHCH3
C4H2      C4H       
C4H6      ! 1,3-butadiene CH2=CHCH=CH2
H2C4O 
C4H4      ! 1-buten-3yne CH#CCH=CH2
IC4H5     ! 1,3-butadien-2-yl (two mesomers) CH2=CHC=CH2
NC4H5     ! 1,3-butadienyl CH2=CHCH=CH
C4H8      ! 1-butene CH3CH2CH=CH2
T2C4H8    ! trans-2-butene CH3CH=CHCH3 (E)-2-butene
C2C4H8    ! cis-2-butene CH3CH=CHCH3 (Z)-2-butene
IC4H3     ! 1-buten-3-yn-2yl (two mesomers) CH2=CC#CH
NC4H3     ! 1-buten-3-ynyl CH=CHC#CH

C6H6      ! benzene        
C6H5O     ! phenoxy radical
C6H5      ! phenyl        
AR        N2
END
REACTIONS

!*********************************************************************
!
! A.KONNOV's detailed reaction mechanism   VER 0.5 CHECKED 
! C6H10, C6H9 deleted 3/5/99, C3H4CY deleted 4/5/99
! NH2OH added 11/11/99 HNOH, H2NO added 12/11/99
! HNNO added 13/11/99 C2H5CHO added 16/03/00 C2H5CO added 17/03/00
! HCNH, CH3NO, H2CNO added 20/06/00
!*********************************************************************
H+H+M=H2+M                     7.000E+17     -1.0         0.0
 H2/0.0/ N2/0.0/ H/0.0/ H2O/14.3/ CO/3.0/ CO2/3.0/ 
H+H+H2=H2+H2                   1.000E+17     -0.6         0.0 
H+H+N2=H2+N2                   5.400E+18     -1.3         0.0 
H+H+H=H2+H                     3.200E+15      0.0         0.0 
O+O+M=O2+M                     1.000E+17     -1.0         0.0 
 O/71.0/ O2/20.0/ NO/5.0/ N2/5.0/ N/5.0/ H2O/5.0/
O+H+M=OH+M                     6.200E+16     -0.6         0.0 
 H2O/5.0/ 
H2+O2=OH+OH                    2.500E+12      0.0     39000.0 
O+H2=OH+H                      5.060E+04      2.67     6290.0 
H+O2=OH+O                      9.750E+13      0.0     14850.0 
H+O2(+M)=HO2(+M)               1.480E+12      0.6         0.0 
    LOW /3.50E+16 -0.41 -1116.0/ 
    TROE /0.5 100000 10/
 AR/0.0/ H2O/10.6/ H2/1.5/ CO2/2.4/ 
H+O2(+AR)=HO2(+AR)             1.480E+12      0.6         0.0
    LOW /7.00E+17 -0.8  0.0/
    TROE /0.45 10 100000/ 
H+OH+M=H2O+M                   2.200E+22     -2.0         0.0 
 H2O/6.4/ AR/0.38/ CO2/1.9/ 
H2+OH=H2O+H                    1.000E+08      1.6      3300.0 
OH+OH=H2O+O                    1.500E+09      1.14      100.0 
HO2+OH=H2O+O2                  2.890E+13      0.0      -500.0 
HO2+O=OH+O2                    1.630E+13      0.0      -445.0 
H+HO2=H2+O2                    4.280E+13      0.0      1411.0 
H+HO2=OH+OH                    1.700E+14      0.0       875.0 
H+HO2=H2O+O                    3.000E+13      0.0      1720.0 
HO2+HO2=H2O2+O2                4.200E+14      0.0     12000.0 
    DUPLICATE
HO2+HO2=H2O2+O2                1.300E+11      0.0     -1640.0 
    DUPLICATE
OH+OH(+M)=H2O2(+M)             7.200E+13     -0.37        0.0 
    LOW /2.2E+19 -0.76 0.0/ 
    TROE /0.5 100000 10/
    H2O/0.0/
OH+OH(+H2O)=H2O2(+H2O)         7.200E+13     -0.37        0.0 
    LOW /1.45E+18 0.0 0.0/ 
H2O2+OH=HO2+H2O                1.000E+12      0.0         0.0 
    DUPLICATE
H2O2+OH=HO2+H2O                5.800E+14      0.0      9560.0 
    DUPLICATE 
H2O2+H=HO2+H2                  1.700E+12      0.0      3755.0 
H2O2+H=H2O+OH                  1.000E+13      0.0      3575.0 
H2O2+O=HO2+OH                  2.800E+13      0.0      6400.0 

N2+O=NO+N                      1.800E+14      0.0     76100.0
N+O2=NO+O                      9.000E+09      1.0      6500.0
NO+M=N+O+M                     9.640E+14      0.0    148300.0
    N2 /1.5/ NO /3.0/ CO2/2.5/
NO+NO=N2+O2                    3.000E+11      0.0     65000.0
N2O(+M)=N2+O(+M)               1.260E+12      0.0     62620.0
    LOW / 4.000E+14 0.0 56640.0/ 
   O2/1.4/ N2/1.7/ H2O/12.0/ NO/3.0/ N2O/3.5/
N2O+O=N2+O2                    1.000E+14      0.0     28200.0
N2O+O=NO+NO                    6.920E+13      0.0     26630.0
N2O+N=N2+NO                    1.000E+13      0.0     20000.0
N2O+NO=N2+NO2                  2.750E+14      0.0     50000.0
NO+O(+M)=NO2(+M)               1.300E+15     -0.75        0.0 
    LOW /4.72E+24 -2.87 1551.0/
    TROE /0.962 10.0 7962.0 /
    AR /0.6/ NO2 /6.2/ NO /1.8/ O2 /0.8/ N2O /4.4/ CO2/0/
    H2O /10.0/
NO+O(+CO2)=NO2(+CO2)           1.300E+15     -0.75        0.0
    LOW /4.0E+22 -2.16 1051.0/
    TROE /0.962 10.0 7962.0 /          
NO2+O=NO+O2                    3.910E+12      0.0      -238.0 
NO2+N=N2O+O                    8.400E+11      0.0         0.0
NO2+N=NO+NO                    1.000E+12      0.0         0.0
NO2+NO=N2O+O2                  1.000E+12      0.0     60000.0
NO2+NO2=NO+NO+O2               3.950E+12      0.0     27590.0 
NO2+NO2=NO3+NO                 1.130E+04      2.58    22720.0 
NO2+O(+M)=NO3(+M)              1.330E+13      0.0         0.0
    LOW / 1.49E+28 -4.08 2467.0 /
    TROE /0.86 10.0 2800.0 /
    H2O/10.0/ O2/0.8/ H2/2.0/ CO2 /0/ 
NO2+O(+CO2)=NO3(+CO2)          1.330E+13      0.0         0.0
    LOW / 1.34E+28 -3.94 2277.0 /
    TROE /0.86 10.0 2800.0 /    
NO3=NO+O2                      2.500E+06      0.0     12120.0
NO3+NO2=NO+NO2+O2              1.200E+11      0.0      3200.0 
NO3+O=NO2+O2                   1.020E+13      0.0         0.0
NO3+NO3=NO2+NO2+O2             5.120E+11      0.0      4870.0
N2O4(+M)=NO2+NO2(+M)           4.050E+18     -1.1     12840.0
     LOW /1.96E+28  -3.8  12840./
     AR/0.8/ N2O4/2.0/ NO2/2.0/
N2O4+O=N2O3+O2                 1.210E+12      0.0         0.0
NO2+NO(+M)=N2O3(+M)            1.600E+09      1.4         0.0
     LOW /1.0E+33  -7.7  0.0/
     N2/1.36/
N2O3+O=NO2+NO2                 2.710E+11      0.0         0.0
N2+M=N+N+M                     1.000E+28     -3.33   225000.0
    N/5/ O/2.2/
NH+M=N+H+M                     2.650E+14      0.0     75500.0
NH+H=N+H2                      3.200E+13      0.0       325.0
NH+N=N2+H                      9.000E+11      0.5         0.0
NH+NH=NNH+H                    5.100E+13      0.0         0.0
NH+NH=NH2+N                    5.950E+02      2.89    -2030.0
NH+NH=N2+H2                    1.000E+08      1.0         0.0
NH2+M=NH+H+M                   3.160E+23     -2.0     91400.0
NH+H2=NH2+H                    1.000E+14      0.0     20070.0
NH2+N=N2+H+H                   6.900E+13      0.0         0.0
NH2+NH=N2H2+H                  1.500E+15     -0.5         0.0
NH2+NH=NH3+N                   1.000E+13      0.0      2000.0
NH3+NH=NH2+NH2                 3.160E+14      0.0     26770.0
NH2+NH2=N2H2+H2                1.000E+13      0.0      1500.0
N2H3+H=NH2+NH2                 5.000E+13      0.0      2000.0
NH3+M=NH2+H+M                  2.200E+16      0.0     93470.0
NH3+M=NH+H2+M                  6.300E+14      0.0     93390.0
NH3+H=NH2+H2                   5.420E+05      2.4      9920.0
NH3+NH2=N2H3+H2                1.000E+11      0.5     21600.0
NNH=N2+H                       3.000E+08      0.0         0.0
     DUPLICATE
NNH+M=N2+H+M                   1.000E+13      0.5      3060.0
     DUPLICATE
NNH+H=N2+H2                    1.000E+14      0.0         0.0
NNH+N=NH+N2                    3.000E+13      0.0      2000.0
NNH+NH=N2+NH2                  2.000E+11      0.5      2000.0
NNH+NH2=N2+NH3                 1.000E+13      0.0         0.0
NNH+NNH=N2H2+N2                1.000E+13      0.0      4000.0
N2H2+M=NNH+H+M                 5.000E+16      0.0     50000.0
     H2O/15.0/ O2/2.0/ N2/2.0/ H2/2.0/
N2H2+M=NH+NH+M                 3.160E+16      0.0     99400.0
     H2O/15.0/ O2/2.0/ N2/2.0/ H2/2.0/
N2H2+H=NNH+H2                  8.500E+04      2.63     -230.0
N2H2+N=NNH+NH                  1.000E+06      2.0         0.0 
N2H2+NH=NNH+NH2                1.000E+13      0.0      6000.0
N2H2+NH2=NH3+NNH               8.800E-02      4.05    -1610.0
N2H3+NH=N2H2+NH2               2.000E+13      0.0         0.0
N2H3+NNH=N2H2+N2H2             1.000E+13      0.0      4000.0
N2H3+M=NH2+NH+M                5.000E+16      0.0     60000.0
N2H3+M=N2H2+H+M                1.000E+16      0.0     37000.0
N2H3+H=N2H2+H2                 1.000E+13      0.0         0.0
N2H3+H=NH+NH3                  1.000E+11      0.0         0.0
N2H3+N=N2H2+NH                 1.000E+06      2.0         0.0 
N2H3+NH2=NH3+N2H2              1.000E+11      0.5         0.0
N2H3+N2H2=N2H4+NNH             1.000E+13      0.0      6000.0
N2H3+N2H3=NH3+NH3+N2           3.000E+12      0.0         0.0
N2H3+N2H3=N2H4+N2H2            1.200E+13      0.0         0.0 
N2H4(+M)=NH2+NH2(+M)           5.000E+14      0.0     60000.0
    LOW/1.50E+15 0.0 39000.0 /
  N2/2.4/ NH3/3.0/ N2H4/4.0/
N2H4+M=N2H3+H+M                1.000E+15      0.0     63600.0
  N2/2.4/ NH3/3.0/ N2H4/4.0/ 
N2H4+H=N2H3+H2                 7.000E+12      0.0      2500.0
N2H4+H=NH2+NH3                 2.400E+09      0.0      3100.0
N2H4+N=N2H3+NH                 1.000E+10      1.0      2000.0
N2H4+NH=NH2+N2H3               1.000E+09      1.5      2000.0
N2H4+NH2=N2H3+NH3              1.800E+06      1.71    -1380.0

N+OH=NO+H                      2.800E+13      0.0         0.0 
N2O+H=N2+OH                    2.200E+14      0.0     16750.0 
N2O+H=NH+NO                    6.700E+22     -2.16    37155.0
N2O+H=NNH+O                    5.500E+18     -1.06    47290.0 
N2O+H=HNNO                     8.000E+24     -4.39    10530.0 
N2O+OH=N2+HO2                  1.000E+14      0.0     30000.0 
HNO+NO=N2O+OH                  8.500E+12      0.0     29580.0 
HNO+NO+NO=HNNO+NO2             1.600E+11      0.0      2090.0 
NH+NO+M=HNNO+M                 1.630E+23     -2.6      1820.0 
HNNO+H=N2O+H2                  2.000E+13      0.0         0.0 
HNNO+H=NH2+NO                  1.000E+12      0.0         0.0 
HNNO+O=N2O+OH                  2.000E+13      0.0         0.0 
HNNO+OH=H2O+N2O                2.000E+13      0.0         0.0 
HNNO+OH=HNOH+NO                1.000E+12      0.0         0.0 
HNNO+NO=N2+HONO                2.600E+11      0.0      1610.0 
HNNO+NO=NNH+NO2                3.200E+12      0.0       540.0 
HNNO+NO=N2O+HNO                1.000E+12      0.0         0.0 
HNNO+NO2=N2O+HONO              1.000E+12      0.0         0.0 
HNNO+NO2=NNH+NO3               1.000E+13      0.0     17000.0 
NO2+H=NO+OH                    1.320E+14      0.0       362.0 
NO2+OH=HO2+NO                  1.810E+13      0.0      6676.0 
NO2+HO2=HONO+O2                4.640E+11      0.0      -479.0
NO2+H2=HONO+H                  7.330E+11      0.0     28800.0  
NO2+NH=N2O+OH                  8.650E+10      0.0     -2270.0
NO2+NH=NO+HNO                  1.245E+11      0.0     -2270.0 
NO3+H=NO2+OH                   6.620E+13      0.0         0.0
NO3+OH=NO2+HO2                 1.210E+13      0.0         0.0
NO3+HO2=HNO3+O2                5.550E+11      0.0         0.0
NO3+HO2=NO2+OH+O2              1.510E+12      0.0         0.0
N2O4+H2O=HONO+HNO3             2.520E+14      0.0     11590.0
N2O3+H2O=HONO+HONO             3.790E+13      0.0      8880.0
H+NO(+M)=HNO(+M)               1.520E+15     -0.41        0.0
    LOW /4.00E+20 -1.75 0.0 /
    H2O/10.0/ O2/1.5/ AR/0.75/ H2/2.0/ CO2/3.0/
HNO+H=NO+H2                    4.460E+11      0.72      655.0 
HNO+OH=NO+H2O                  1.300E+07      1.88     -956.0 
HNO+O=OH+NO                    5.000E+11      0.5      2000.0
HNO+O=NO2+H                    5.000E+10      0.0      2000.0
HNO+O2=NO+HO2                  2.200E+10      0.0      9140.0
HNO+N=NO+NH                    1.000E+11      0.5      2000.0 
HNO+N=H+N2O                    5.000E+10      0.5      3000.0
HNO+NH=NH2+NO                  5.000E+11      0.5         0.0
HNO+NH2=NH3+NO                 2.000E+13      0.0      1000.0 
HNO+HNO=N2O+H2O                3.630E-03      3.98     1190.0 
HNO+HNO=HNOH+NO                2.000E+08      0.0      4170.0 
HNO+NO2=HONO+NO                6.020E+11      0.0      2000.0
NO+OH(+M)=HONO(+M)             2.000E+12     -0.05     -721.0 
    LOW / 5.08E+23 -2.51 -67.6 /
    TROE /0.62 10.0 100000.0 /
    H2O/10.0/ O2/2.0/ AR/0.75/ H2/2.0/ CO2/0.0/  
NO+OH(+CO2)=HONO(+CO2)         2.000E+12     -0.05     -721.0
    LOW / 1.70E+23 -2.3 -246.0 /
    TROE /0.62 10.0 100000.0 /    
NO2+H+M=HONO+M                 1.400E+18     -1.5       900.0
HONO+H=HNO+OH                  5.640E+10      0.86     4970.0 
HONO+H=NO+H2O                  8.120E+06      1.89     3840.0 
HONO+O=OH+NO2                  1.200E+13      0.0      5960.0
HONO+OH=H2O+NO2                1.690E+12      0.0      -517.0 
HONO+NH=NH2+NO2                1.000E+13      0.0         0.0 
HONO+HONO=H2O+NO2+NO           1.000E+13      0.0      8540.0
HONO+NH2=NO2+NH3               5.000E+12      0.0         0.0
NO2+OH(+M)=HNO3(+M)            2.410E+13      0.0         0.0
    LOW / 6.42E+32 -5.49 2350.0 /
    TROE /1.0 10.0 1168.0 /
    H2O/10.0/ O2/2.0/ AR/0.75/ H2/2.0/ CO2/0.0/ 
NO2+OH(+CO2)=HNO3(+CO2)        2.410E+13      0.0         0.0
    LOW / 5.80E+32 -5.4 2186.0 /
    TROE /1.0 10.0 1168.0 /       
NO+HO2+M=HNO3+M                1.500E+24     -3.5      2200.0
HNO3+H=H2+NO3                  5.560E+08      1.53    16400.0 
HNO3+H=H2O+NO2                 6.080E+01      3.29     6290.0 
HNO3+H=OH+HONO                 3.820E+05      2.3      6980.0 
HNO3+OH=NO3+H2O                1.030E+10      0.0     -1240.0
NH3+O=NH2+OH                   1.100E+06      2.1      5210.0 
NH3+OH=NH2+H2O                 5.000E+07      1.6       950.0 
NH3+HO2=NH2+H2O2               3.000E+11      0.0     22000.0
NH2+HO2=NH3+O2                 1.650E+04      1.55     2027.0 
NH2+O=H2+NO                    5.000E+12      0.0         0.0 
NH2+O=HNO+H                    4.500E+13      0.0         0.0 
NH2+O=NH+OH                    7.000E+12      0.0         0.0 
NH2+OH=NH+H2O                  9.000E+07      1.5      -460.0 
NH2+OH=NH2OH                   1.790E+13      0.2         0.0 
NH2+HO2=HNO+H2O                5.680E+15     -1.12      707.0 
NH2+HO2=H2NO+OH                2.910E+17     -1.32     1248.0 
NH2+O2=HNO+OH                  1.000E+13      0.0     26290.0 
NH2+O2=H2NO+O                  6.000E+13      0.0     29880.0 
NH2+NO=NNH+OH                  2.290E+10      0.425    -814.0
NH2+NO=N2+H2O                  2.770E+20     -2.65     1258.0 
NH2+NO=H2+N2O                  1.000E+13      0.0     33700.0 
NH2+NO2=N2O+H2O                1.620E+16     -1.44      270.0
NH2+NO2=H2NO+NO                6.480E+16     -1.44      270.0 
NH+O=NO+H                      7.000E+13      0.0         0.0 
NH+O=N+OH                      7.000E+12      0.0         0.0 
NH+OH=HNO+H                    2.000E+13      0.0         0.0 
NH+OH=N+H2O                    2.000E+09      1.2         0.0 
NH+OH=NO+H2                    2.000E+13      0.0         0.0 
NH+HO2=HNO+OH                  1.000E+13      0.0      2000.0
NH+O2=HNO+O                    4.000E+13      0.0     17880.0 
NH+O2=NO+OH                    4.500E+08      0.79     1190.0 
NH+H2O=HNO+H2                  2.000E+13      0.0     13850.0 
NH+N2O=N2+HNO                  2.000E+12      0.0      6000.0
NNH+O=NH+NO                    2.000E+14      0.0      4000.0 
NH+NO=N2+OH                    6.100E+13     -0.50      120.0 
N2H4+O=N2H2+H2O                8.500E+13      0.0      1200.0
N2H4+O=N2H3+OH                 2.500E+12      0.0      1200.0 
N2H4+OH=N2H3+H2O               3.000E+10      0.68     1290.0
N2H4+OH=NH3+H2NO               3.670E+13      0.0         0.0 
N2H4+HO2=N2H3+H2O2             4.000E+13      0.0      2000.0 
N2H3+O=N2H2+OH                 2.000E+13      0.0      1000.0
N2H3+O=NNH+H2O                 3.160E+11      0.5         0.0 
N2H3+O=NH2+HNO                 1.000E+13      0.0         0.0
N2H3+OH=N2H2+H2O               3.000E+10      0.68     1290.0
N2H3+OH=NH3+HNO                1.000E+12      0.0     15000.0 
N2H3+O2=N2H2+HO2               3.000E+12      0.0         0.0
N2H3+HO2=N2H2+H2O2             1.000E+13      0.0      2000.0 
N2H3+HO2=N2H4+O2               8.000E+12      0.0         0.0
N2H3+NO=HNO+N2H2               1.000E+12      0.0         0.0 
N2H2+O=NH2+NO                  1.000E+13      0.0         0.0 
N2H2+O=NNH+OH                  2.000E+13      0.0      1000.0 
N2H2+OH=NNH+H2O                5.920E+01      3.4     -1360.0 
N2H2+HO2=NNH+H2O2              1.000E+13      0.0      2000.0 
N2H2+NO=N2O+NH2                3.000E+10      0.0         0.0 
NNH+O=N2+OH                    1.700E+16     -1.23      500.0
NNH+OH=N2+H2O                  2.400E+22     -2.88     2444.0 
NNH+O2=N2+HO2                  1.200E+12     -0.34      150.0
NNH+O2=N2O+OH                  2.900E+11     -0.34      150.0
NNH+HO2=N2+H2O2                1.000E+13      0.0      2000.0 
NNH+NO=N2+HNO                  5.000E+13      0.0         0.0 
NH2OH+OH=HNOH+H2O              2.500E+13      0.0      4250.0
H2NO+M=H2+NO+M                 7.830E+27     -4.29    60300.0 
  H2O/10.0/
H2NO+M=HNO+H+M                 2.800E+24     -2.83    64915.0 
  H2O/10.0/
H2NO+M=HNOH+M                  1.100E+29     -3.99    43980.0
  H2O/10.0/
H2NO+H=HNO+H2                  3.000E+07      2.0      2000.0 
H2NO+H=NH2+OH                  5.000E+13      0.0         0.0 
H2NO+O=HNO+OH                  3.000E+07      2.0      2000.0 
H2NO+OH=HNO+H2O                2.000E+07      2.0      1000.0 
H2NO+HO2=HNO+H2O2              2.900E+04      2.69    -1600.0 
H2NO+NH2=HNO+NH3               3.000E+12      0.0      1000.0 
H2NO+O2=HNO+HO2                3.000E+12      0.0     25000.0
H2NO+NO=HNO+HNO                2.000E+07      2.0     13000.0 
H2NO+NO2=HONO+HNO              6.000E+11      0.0      2000.0 
HNOH+M=HNO+H+M                 2.000E+24     -2.84    58935.0 
   H2O/10.0/
HNOH+H=HNO+H2                  4.800E+08      1.5       380.0 
HNOH+H=NH2+OH                  4.000E+13      0.0         0.0 
HNOH+O=HNO+OH                  7.000E+13      0.0         0.0 
    DUPLICATE
HNOH+O=HNO+OH                  3.300E+08      1.5      -360.0 
    DUPLICATE
HNOH+OH=HNO+H2O                2.400E+06      2.0     -1190.0
HNOH+HO2=HNO+H2O2              2.900E+04      2.69    -1600.0
HNOH+NH2=HNO+NH3               1.800E+06      1.94    -1150.0
HNOH+NO2=HONO+HNO              6.000E+11      0.0      2000.0
HNOH+O2=HNO+HO2                3.000E+12      0.0     25000.0
HNOH+HNO=NH2OH+NO              1.000E+12      0.0      3000.0 

!END
CO+HO2=CO2+OH                  1.500E+14      0.0     23650.0 
CO+OH=CO2+H                    1.170E+07      1.354    -725.0 
CO+O+M=CO2+M                   6.160E+14      0.0      3000.0 
 H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/0.88/
 CH4/3.2/ CH3OH/7.5/ 
CO+O2=CO2+O                    2.500E+12      0.0     47800.0 
HCO+M=H+CO+M                   1.560E+14      0.0     15760.0 
 H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/1.0/
 CH4/3.2/ CH3OH/7.5/ 
HCO+OH=CO+H2O                  1.000E+14      0.0         0.0 
HCO+O=CO+OH                    3.000E+13      0.0         0.0 
HCO+O=CO2+H                    3.000E+13      0.0         0.0 
HCO+H=CO+H2                    9.000E+13      0.0         0.0 
HCO+O2=CO+HO2                  2.700E+13      0.0      1190.0 
HCO+CH3=CO+CH4                 1.200E+14      0.0         0.0 
HCO+HO2=CO2+OH+H               3.000E+13      0.0         0.0 
HCO+HCO=CH2O+CO                3.000E+13      0.0         0.0 
HCO+HCO=H2+CO+CO               2.200E+13      0.0         0.0 
CH4(+M)=CH3+H(+M)              2.400E+16      0.0    104913.0 
    LOW /4.5E+17 0.0 90800/ 
    TROE /1.0 10.0 1350.0 7830.0/ 
    CH4/0.0/ H2/2.0/ CO/2.0/ CO2/3.0/ H2O/5.0/ 
CH4(+CH4)=CH3+H(+CH4)          2.400E+16      0.0    104913.0 
    LOW /8.4E+18 0.0 90800/ 
    TROE /0.31 2210.0 90/ 
CH4+HO2=CH3+H2O2               9.000E+12      0.0     24641.0 
CH4+OH=CH3+H2O                 1.548E+07      1.83     2774.0 
CH4+O=CH3+OH                   7.200E+08      1.56     8485.0 
CH4+H=CH3+H2                   1.300E+04      3.0      8050.0 
CH4+CH2=CH3+CH3                4.300E+12      0.0     10038.0 
CH4+O2=CH3+HO2                 4.000E+13      0.0     56900.0 
CH3+M=CH2+H+M                  2.720E+36     -5.31   117100.0 
 H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/1.0/
 CH4/3.2/ CH3OH/7.5/ 
CH3+M=CH+H2+M                  1.000E+16      0.0     85240.0 
 H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/0.88/
 CH4/3.2/ CH3OH/7.5/ 
CH3+HO2=CH3O+OH                1.800E+13      0.0         0.0 
CH3+OH=CH2OH+H                 2.640E+19     -1.8      8068.0 
CH3+OH=CH3O+H                  5.740E+12     -0.23    13931.0 
CH3+OH=CH2+H2O                 8.900E+18     -1.8      8067.0 
CH3+OH=CH2O+H2                 3.190E+12     -0.53    10810.0 
CH3+O=H+CH2O                   8.430E+13      0.0         0.0 
CH3+O2=CH2O+OH                 3.400E+11      0.0      8940.0 
CH3+O2=CH3O+O                  1.320E+14      0.0     31400.0 
CH3+CH3=C2H5+H                 5.000E+12      0.099   10600.0
CH3+CH3(+M)=C2H6(+M)           9.210E+16     -1.174     636.0 
    LOW /1.13E+36 -5.246 1705/ 
    TROE /0.405 1120.0 69.6/
    H2/2.0/ CO/2.0/ CO2/3.0/ H2O/5.0/ 
CH3+CH3O=CH4+CH2O              2.409E+13      0.0         0.0 
CH3+CH2OH=CH4+CH2O             8.500E+13      0.0         0.0 
CH3+H=SCH2+H2                  6.000E+13      0.0     15100.0 
CH3+O2(+M)=CH3O2(+M)           7.800E+08      1.2         0.0 
    LOW/5.8E+25 -3.30 0.0/
    TROE /0.495 2325.5 10/
CH3+CH3=C2H4+H2                1.000E+14      0.0     32000.0 
CH3+OH=SCH2+H2O                7.200E+13      0.0      2780.0 
CH2+OH=CH2O+H                  2.500E+13      0.0         0.0 
CH2+O=CO+H2                    4.800E+13      0.0         0.0 
CH2+O=CO+H+H                   7.200E+13      0.0         0.0 
CH2+O=CH+OH                    3.000E+14      0.0     11920.0 
CH2+O=HCO+H                    3.000E+13      0.0         0.0 
CH2+H=CH+H2                    3.120E+13      0.0     -1340.0 
CH2+O2=HCO+OH                  4.300E+10      0.0      -500.0 
CH2+O2=CO2+H2                  6.900E+11      0.0       500.0 
CH2+O2=CO2+H+H                 1.600E+12      0.0      1000.0 
CH2+O2=CO+H2O                  1.900E+10      0.0     -1000.0 
CH2+O2=CO+OH+H                 8.600E+10      0.0      -500.0 
CH2+O2=CH2O+O                  5.000E+13      0.0      9000.0 
CH2+CO2=CH2O+CO                1.100E+11      0.0      1000.0 
CH2+CH2=C2H2+H2                1.580E+15      0.0     11950.0 
CH2+CH2=C2H2+H+H               2.000E+14      0.0     11000.0 
CH2+CH2=CH3+CH                 2.400E+14      0.0      9940.0
CH2+CH2=C2H3+H                 2.000E+13      0.0         0.0
CH2+CH3=C2H4+H                 4.200E+13      0.0         0.0 
CH2+CH=C2H2+H                  4.000E+13      0.0         0.0 
CH2+C=CH+CH                    1.620E+12      0.67    46800.0 
CH2+M=C+H2+M                   1.600E+14      0.0     64000.0 
CH2+M=CH+H+M                   5.600E+15      0.0     89600.0 
SCH2+M=CH2+M                   6.000E+12      0.0         0.0 
 H2/2.5/ H2O/5.0/ CO/1.875/ CO2/3.75/ AR/0.6/ CH4/1.2/
 C2H2/8.0/ C2H4/4.0/ C2H6/3.6/ H/33.3/ 
SCH2+O2=CO+OH+H                3.000E+13      0.0         0.0 
SCH2+H=CH+H2                   3.000E+13      0.0         0.0 
SCH2+O=CO+H+H                  1.500E+13      0.0         0.0 
SCH2+O=CO+H2                   1.500E+13      0.0         0.0 
SCH2+OH=CH2O+H                 3.000E+13      0.0         0.0 
SCH2+HO2=CH2O+OH               3.000E+13      0.0         0.0 
SCH2+H2O2=CH3O+OH              3.000E+13      0.0         0.0 
SCH2+H2O=>CH3OH                1.800E+13      0.0         0.0 
SCH2+CH2O=CH3+HCO              1.200E+12      0.0         0.0 
SCH2+HCO=CH3+CO                1.800E+13      0.0         0.0 
SCH2+CH3=C2H4+H                1.800E+13      0.0         0.0 
SCH2+CH4=CH3+CH3               4.000E+13      0.0         0.0 
SCH2+C2H6=CH3+C2H5             1.200E+14      0.0         0.0 
SCH2+CO2=CH2O+CO               3.000E+12      0.0         0.0 
SCH2+CH2CO=C2H4+CO             1.600E+14      0.0         0.0 
CH+OH=HCO+H                    3.000E+13      0.0         0.0 
CH+O=CO+H                      4.000E+13      0.0         0.0 
CH+O=C+OH                      1.520E+13      0.0      4730.0 
H2O+C=CH+OH                    7.800E+11      0.67    39300.0 
CH+O2=HCO+O                    4.900E+13      0.0         0.0 
CH+O2=CO+OH                    4.900E+13      0.0         0.0 
CH+CO2=HCO+CO                  3.220E-02      4.44    -3530.0 
CH+CH4=C2H4+H                  3.900E+14     -0.4         0.0 
CH+CH3=C2H3+H                  3.000E+13      0.0         0.0 
CH2+OH=CH+H2O                  1.130E+07      2.0      3000.0 
CH+H=C+H2                      7.900E+13      0.0       160.0 
CH+H2O=CH2O+H                  1.170E+15     -0.75        0.0 
CH+H2O=CH2OH                   5.700E+12      0.0      -760.0 
CH+CH2O=CH2CO+H                1.000E+14      0.0      -515.0 
CH3O+M=CH2O+H+M                5.400E+13      0.0     13500.0 
 H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/0.88/
 CH4/3.2/ CH3OH/7.5/ 
CH3O+HO2=CH2O+H2O2             3.000E+11      0.0         0.0 
CH3O+OH=CH2O+H2O               1.800E+13      0.0         0.0 
CH3O+O=CH2O+OH                 1.800E+12      0.0         0.0 
CH3O+H=CH2O+H2                 1.800E+13      0.0         0.0 
CH3O+O2=CH2O+HO2               2.200E+10      0.0      1750.0 
CH3O+CH2O=CH3OH+HCO            1.000E+11      0.0      2980.0 
CH3O+CO=CH3+CO2                6.810E-18      9.2     -2850.0 
CH3O+HCO=CH3OH+CO              9.000E+13      0.0         0.0 
CH3O+C2H5=CH2O+C2H6            2.410E+13      0.0         0.0 
CH3O+C2H3=CH2O+C2H4            2.410E+13      0.0         0.0 
CH3O+C2H4=CH2O+C2H5            1.200E+11      0.0      6750.0 
CH3O+H=CH2OH+H                 3.400E+06      1.6         0.0 
CH3O+H=SCH2+H2O                1.000E+12      0.0         0.0 
CH2O+M=HCO+H+M                 5.000E+35     -5.54    96680.0 
 H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/1.0/
 CH4/3.2/ CH3OH/7.5/ 
CH2O+M=CO+H2+M                 1.100E+36     -5.54    96680.0 
 H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/1.0/
 CH4/3.2/ CH3OH/7.5/ 
CH2O+HO2=HCO+H2O2              4.110E+04      2.5     10210.0 
CH2O+OH=HCO+H2O                3.433E+09      1.18     -447.0 
CH2O+O=HCO+OH                  4.100E+11      0.57     2760.0 
CH2O+H=HCO+H2                  1.260E+08      1.62     2166.0 
CH2O+O2=HCO+HO2                6.000E+13      0.0     40650.0 
CH2O+CH3=HCO+CH4               7.800E-08      6.1      1970.0 
C2H6(+M)=C2H5+H(+M)            8.850E+20     -1.228  102210.0 
    LOW /6.90E+42 -6.431 107175.0/
    SRI /47.61 16182.0 3371.0/ 
 H2/2.0/ H2O/6.0/ CH4/2.0/ CO/1.5/ CO2/2.0/ C2H6/3.0/ AR/0.7/
C2H6+HO2=C2H5+H2O2             1.330E+13      0.0     20535.0 
C2H6+OH=C2H5+H2O               7.200E+06      2.0       870.0 
C2H6+O=C2H5+OH                 1.000E+09      1.5      5800.0 
C2H6+H=C2H5+H2                 1.400E+09      1.5      7400.0 
C2H6+H=CH3+CH4                 5.400E+04      0.0     11630.0 
C2H6+O2=C2H5+HO2               6.000E+13      0.0     52000.0 
C2H6+CH3=C2H5+CH4              1.470E-07      6.0      6060.0 
C2H6+CH2=CH3+C2H5              6.500E+12      0.0      7911.0 
C2H6+C2H3=C2H4+C2H5            8.566E-02      4.14     2543.0 
C2H6+HCO=CH2O+C2H5             4.700E+04      2.72    18235.0 
C2H5(+M)=C2H4+H(+M)            1.110E+10      1.037   36767.0 
    LOW /4.0E+33 -4.99 40000.0/
    TROE /0.832 10 1203.0/
 H2/2.0/ CO/2.0/ CO2/3.0/ H2O/5.0/ CH4/2.0/ C2H6/0.0/ AR/0.7/
C2H5(+C2H6)=C2H4+H(+C2H6)      8.200E+13      0.0     39880.0 
    LOW /1.0E+18 0.0 33380.0/
    TROE /0.75 97.0 1379.0/
C2H5+HO2=C2H4+H2O2             1.800E+12      0.0         0.0 
C2H5+OH=C2H4+H2O               2.409E+13      0.0         0.0 
C2H5+OH=>CH3+CH2O+H            2.409E+13      0.0         0.0 
C2H5+O=CH2O+CH3                4.240E+13      0.0         0.0 
C2H5+O=CH3HCO+H                5.300E+13      0.0         0.0 
C2H5+O=C2H4+OH                 3.460E+13      0.0         0.0 
C2H5+H=C2H4+H2                 1.700E+12      0.0         0.0 
C2H5+O2=C2H4+HO2               2.560E+19     -2.77     1980.0 
C2H5+CH3=C2H4+CH4              1.100E+12      0.0         0.0 
C2H5+C2H5=C2H4+C2H6            1.400E+12      0.0         0.0 
C2H5+HO2=C2H5O+OH              3.000E+13      0.0         0.0 
C2H4+M=C2H2+H2+M               3.500E+16      0.0     71530.0 
    H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/1.0/
    CH4/3.2/ CH3OH/7.5/ 
C2H4+M=C2H3+H+M                2.600E+17      0.0     96570.0 
    H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/1.0/
    CH4/3.2/ CH3OH/7.5/ 
C2H4+OH=C2H3+H2O               5.530E+05      2.31     2900.0 
C2H4+O=CH3+HCO                 8.100E+06      1.88      180.0 
C2H4+H=C2H3+H2                 4.490E+07      2.12    13366.0 
C2H4+O2=C2H3+HO2               4.000E+13      0.0     61500.0 
C2H4+C2H4=C2H5+C2H3            1.860E+14      0.0     64200.0 
C2H4+CH3=C2H3+CH4              4.200E+12      0.0     11100.0 
C2H4+O=CH2HCO+H                4.700E+06      1.88      180.0 
C2H4+O=CH2O+CH2                3.000E+04      1.88      180.0 
C2H4+O=CH2CO+H2                6.700E+05      1.88      180.0 
C2H4+O=C2H3+OH                 1.510E+07      1.91     3790.0 
C2H4+OH=CH2O+CH3               2.000E+12      0.0       960.0 
C2H4+OH(+M)=PC2H5O(+M)         5.420E+12      0.0         0.0
      LOW /1.19E+27 -3.1 0.0/ 
C2H4+HO2=C2H3+H2O2             1.120E+13      0.0     30400.0 
C2H4+CH3O=C2H3+CH3OH           1.000E+11      0.0     10000.0 
C2H3(+M)=C2H2+H(+M)            2.100E+14      0.0     39740.0 
    LOW /4.15E+41 -7.5 45500.0/
    TROE /0.65 100000 10/
 H2/2.0/ CO/2.0/ CO2/3.0/ H2O/5.0/ CH4/2.0/ C2H6/3.0/ AR/0.7/
C2H3+HO2=>CH3+CO+OH            3.000E+13      0.0         0.0 
C2H3+OH=C2H2+H2O               3.000E+13      0.0         0.0 
C2H3+H=C2H2+H2                 1.200E+13      0.0         0.0 
C2H3+O=CH3+CO                  1.000E+13      0.0         0.0 
C2H3+O2=CH2O+HCO               1.700E+29     -5.312    6500.0 
C2H3+CH=CH2+C2H2               5.000E+13      0.0         0.0 
C2H3+CH3=C2H2+CH4              2.050E+13      0.0         0.0 
C2H3+C2H=C2H2+C2H2             3.000E+13      0.0         0.0 
C2H3+HCO=C2H4+CO               9.034E+13      0.0         0.0 
C2H3+CH2O=C2H4+HCO             5.420E+03      2.81     5862.0 
C2H3+C2H3=C2H2+C2H4            1.450E+13      0.0         0.0 
C2H3+O=C2H2+OH                 1.000E+13      0.0         0.0 
C2H3+O=CH2+HCO                 1.000E+13      0.0         0.0 
C2H3+O=CH2CO+H                 1.000E+13      0.0         0.0 
C2H3+OH=CH3HCO                 3.000E+13      0.0         0.0 
C2H3+O2=C2H2+HO2               5.190E+15     -1.26     3310.0 
      DUPLICATE
C2H3+O2=C2H2+HO2               2.120E-06      6.0      9484.0
      DUPLICATE
C2H3+O2=CH2HCO+O               3.500E+14     -0.61     5260.0
C2H3+CH2=C2H2+CH3              3.000E+13      0.0         0.0 
C2H2=C2H+H                     2.373E+32     -5.28   130688.0 
C2H2+O2=HCCO+OH                2.000E+08      1.5     30100.0 
C2H2+O2=C2H+HO2                1.200E+13      0.0     74520.0 
C2H2+OH=C2H+H2O                3.385E+07      2.0     14000.0 
C2H2+OH=CH2CO+H                1.100E+13      0.0      7170.0 
C2H2+O=CH2+CO                  1.200E+06      2.1      1570.0 
C2H2+O=HCCO+H                  5.000E+06      2.1      1570.0 
C2H2+CH3=C2H+CH4               1.800E+11      0.0     17290.0 
C2H2+O=C2H+OH                  3.000E+14      0.0     25000.0 
C2H2+OH=CH3+CO                 4.830E-04      4.0     -2000.0 
C2H2+HO2=CH2CO+OH              6.100E+09      0.0      7950.0 
C2H2+O2=HCO+HCO                4.000E+12      0.0     28000.0 
C2H+OH=HCCO+H                  2.000E+13      0.0         0.0 
C2H+OH=C2+H2O                  4.000E+07      2.0      8000.0 
C2H+O=CO+CH                    1.450E+13      0.0       460.0 
C2H+O2=HCO+CO                  9.000E+12      0.0         0.0 
C2H+H2=C2H2+H                  7.880E+05      2.39      346.0 
C2H+O2=CO+CO+H                 9.000E+12      0.0         0.0 
C2H+O2=HCCO+O                  6.000E+11      0.0         0.0 
CH2CO(+M)=CH2+CO(+M)           3.000E+14      0.0     71000.0 
    LOW /2.300E+15   0.0  57600.0/ 
    H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/1.0/
    CH4/3.2/ CH3OH/7.5/ 
CH2CO+O2=CH2O+CO2              2.000E+13      0.0     61500.0 
CH2CO+HO2=>CH2O+CO+OH          6.000E+11      0.0     12738.0 
CH2CO+O=HCCO+OH                1.000E+13      0.0      8000.0 
CH2CO+OH=CH2OH+CO              1.000E+13      0.0         0.0 
CH2CO+H=CH3+CO                 3.280E+10      0.851    2840.0 
CH2CO+CH3=C2H5+CO              2.400E+12      0.0      8000.0 
CH2CO+CH2=C2H4+CO              2.900E+12      0.0      3800.0 
CH2CO+CH2=HCCO+CH3             3.600E+13      0.0     11000.0
CH2CO+CH3=HCCO+CH4             7.500E+12      0.0     13000.0
CH2CO+OH=CH2O+HCO              2.800E+13      0.0         0.0 
CH2CO+H=HCCO+H2                1.800E+14      0.0      8600.0 
CH2CO+O=HCO+HCO                7.500E+11      0.0      1350.0 
CH2CO+O=HCO+CO+H               7.500E+11      0.0      1350.0 
CH2CO+O=CH2O+CO                7.500E+11      0.0      1350.0 
CH2CO+OH=HCCO+H2O              7.500E+12      0.0      2000.0 
HCCO+M=CH+CO+M                 6.000E+15      0.0     58821.0 
    H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/0.88/
    CH4/3.2/ CH3OH/7.5/ 
HCCO+OH=HCO+CO+H               1.000E+13      0.0         0.0 
HCCO+OH=C2O+H2O                3.000E+13      0.0         0.0 
HCCO+O=CO+CO+H                 1.000E+14      0.0         0.0 
HCCO+O=CH+CO2                  2.950E+13      0.0      1110.0 
HCCO+H=CH2+CO                  1.500E+14      0.0         0.0 
HCCO+O2=CO2+CO+H               5.400E+11      0.0       850.0 
HCCO+CH2=C2H+CH2O              1.000E+13      0.0      2000.0 
HCCO+CH2=C2H3+CO               3.000E+13      0.0         0.0 
HCCO+CH3=C2H4+CO               2.000E+12      0.0         0.0
HCCO+CH=CO+C2H2                5.000E+13      0.0         0.0 
HCCO+HCCO=CO+C2H2+CO           1.000E+13      0.0         0.0 
HCCO+OH=HCO+HCO                1.000E+13      0.0         0.0 
HCCO+O2=CO+CO+OH               5.400E+11      0.0       850.0 
HCCO+O2=CO2+HCO                5.400E+11      0.0       850.0 
CH3OH(+M)=CH3+OH(+M)           1.700E+16      0.0     90885.0 
    LOW /6.60E+16 0.0 65730.0/ 
    TROE /0.82 200.0 1438.0/
CH3OH+HO2=CH2OH+H2O2           9.640E+10      0.0     12580.0 
CH3OH+OH=CH2OH+H2O             1.440E+06      2.0      -840.0 
CH3OH+OH=CH3O+H2O              1.000E+13      0.0      1700.0 
CH3OH+O=CH2OH+OH               1.630E+13      0.0      5030.0 
CH3OH+H=CH2OH+H2               1.640E+07      2.0      4520.0 
CH3OH+CH3=CH2OH+CH4            3.190E+01      3.17     7172.0 
CH3OH+CH3=CH3O+CH4             1.450E+01      3.1      6935.0 
CH3OH+C2H5=C2H6+CH3O           1.440E+01      3.1      8942.0 
CH3OH+H=CH3+H2O                2.000E+14      0.0      5300.0 
CH3OH+O=CH3O+OH                1.000E+13      0.0      4680.0 
CH3OH+CH3=C2H6+OH              2.000E+12      0.0     15000.0 
CH3OH+CH3O=CH2OH+CH3OH         3.000E+11      0.0      4070.0 
CH3OH(+M)=CH2OH+H(+M)          1.380E+16      0.0     95950.0 
    LOW /5.35E+16 0.0 70800.0/ 
    TROE /0.82 200.0 1438.0/
CH3OH+H=H2+CH3O                4.000E+13      0.0      6095.0 
CH3OH+O2=CH2OH+HO2             2.050E+13      0.0     44900.0 
CH3OH+C2H5=C2H6+CH2OH          3.190E+01      3.2      9161.0 
CH2OH+M=CH2O+H+M               1.140E+43     -8.0     43000.0 
    H2O/16.0/ CH4/3.0/ CO2/3.75/  CO/1.875/  H2/2.5/   CH3OH/6.0/
CH2OH+H=CH2O+H2                1.000E+13      0.0         0.0 
CH2OH+O2=CH2O+HO2              1.500E+15     -1.0         0.0 
    duplicate
CH2OH+O2=CH2O+HO2              7.200E+13      0.0      3570.0 
    duplicate
H+CH2OH=SCH2+H2O               1.000E+12      0.0         0.0 
CH2OH+O=CH2O+OH                9.000E+13      0.0         0.0 
CH2OH+OH=CH2O+H2O              1.000E+13      0.0         0.0 
CH2OH+HO2=CH2O+H2O2            1.210E+13      0.0         0.0 
CH2OH+CH2OH=CH3OH+CH2O         4.820E+12      0.0         0.0 
CH2OH+CH2OH=CH2O+CH2O+H2       1.000E+15     -0.7         0.0 
CH2OH+HCO=CH3OH+CO             1.210E+14      0.0         0.0 
CH2OH+CH2O=CH3OH+HCO           5.490E+03      2.8      5900.0 
CH2OH+CH3O=CH3OH+CH2O          2.400E+13      0.0         0.0 
CH3O+CH3O=CH3OH+CH2O           2.320E+13      0.0         0.0
CH3HCO=CH3+HCO                 7.100E+15      0.0     81790.0 
CH3HCO+HO2=CH3CO+H2O2          3.000E+12      0.0     12000.0 
CH3HCO+OH=CH3CO+H2O            2.300E+10      0.73    -1100.0 
CH3HCO+O=CH3CO+OH              5.800E+12      0.0      1800.0 
CH3HCO+H=CH3CO+H2              4.100E+09      1.16     2400.0 
CH3HCO+O2=CH3CO+HO2            3.000E+13      0.0     39200.0 
CH3HCO+CH3=CH3CO+CH4           7.600E+00      3.4      3740.0 
CH3HCO+H=CH2HCO+H2             7.000E+08      1.5      7400.0 
CH3HCO+O=CH2HCO+OH             5.000E+08      1.5      5800.0 
CH3HCO+OH=CH2HCO+H2O           2.000E+14      0.0      6000.0 
CH3HCO+HO2=CH2HCO+H2O2         3.000E+13      0.0     15000.0 
CH3HCO+CH2=CH3CO+CH3           1.660E+12      0.0      3510.0 
CH3HCO+CH3=CH2HCO+CH4          1.580E+00      4.0      7720.0 
CH3HCO+CH3O=CH3CO+CH3OH        5.000E+12      0.0         0.0 
CH3HCO+C2H5=CH3CO+C2H6         1.260E+12      0.0      8500.0 
CH3HCO+C2H3=CH3CO+C2H4         8.130E+10      0.0      3680.0 
CH2HCO=CH3CO                   1.600E+11      0.0     21600.0 
CH3HCO+CH2HCO=CH3CO+CH3HCO     3.000E+12      0.0     11200.0 
CH3CO(+M)=CH3+CO(+M)           2.800E+13      0.0     17150.0 
    LOW /6.0E+15 0.0 14070.0/ 
    TROE /0.5 100000 10/
    H2/2.5/ H2O/6.2/ CO/1.875/ CO2/3.75/ AR/0.88/
    CH4/3.2/ CH3OH/7.5/ 
CH3CO+H=CH2CO+H2               1.150E+13      0.0         0.0 
CH3CO+H=CH3+HCO                2.150E+13      0.0         0.0 
CH3CO+O=CH2CO+OH               4.000E+13      0.0         0.0 
CH3CO+O=CH3+CO2                1.500E+14      0.0         0.0 
CH3CO+CH3=C2H6+CO              3.300E+13      0.0         0.0 
CH3CO+CH3=CH4+CH2CO            6.100E+12      0.0         0.0
CH2HCO+H=CH2CO+H2              2.000E+13      0.0         0.0 
CH2HCO+O2=CH2O+OH+CO           1.800E+10      0.0         0.0 
CH2HCO+O2=CH2CO+HO2            1.500E+11      0.0         0.0 
CH2HCO=CH2CO+H                 1.580E+13      0.0     35200.0 
C2H5O=CH3+CH2O                 1.000E+15      0.0     21600.0 
C2H5O+O2=CH3HCO+HO2            3.600E+10      0.0      1090.0 
C2H5O=CH3HCO+H                 2.000E+14      0.0     23300.0 
C2H5O+OH=CH3HCO+H2O            1.000E+14      0.0         0.0 
C2H5O+H=CH3HCO+H2              1.000E+14      0.0         0.0 
C2H5O+O=CH3HCO+OH              1.210E+14      0.0         0.0 
C2H5O+HO2=CH3HCO+H2O2          1.000E+14      0.0         0.0 
C2H5O+C2H5O=C2H5OH+CH3HCO      5.000E+13      0.0         0.0
C2H5O+PC2H5O=C2H5OH+CH3HCO     5.000E+13      0.0         0.0
C2H5O+SC2H5O=C2H5OH+CH3HCO     5.000E+13      0.0         0.0
SC2H5O+M=CH3HCO+H+M            5.000E+13      0.0     21860.0 
SC2H5O+H=CH3HCO+H2             2.000E+13      0.0         0.0 
SC2H5O+OH=CH3HCO+H2O           1.500E+13      0.0         0.0 
SC2H5O+O=CH3HCO+OH             9.040E+13      0.0         0.0 
SC2H5O+O2=CH3HCO+HO2           8.400E+15     -1.20        0.0 
    duplicate
SC2H5O+O2=CH3HCO+HO2           4.800E+14      0.0      5000.0 
    duplicate
SC2H5O+HO2=CH3HCO+H2O2         1.000E+13      0.0         0.0 
SC2H5O+SC2H5O=C2H5OH+CH3HCO    3.500E+13      0.0         0.0
SC2H5O+PC2H5O=C2H5OH+CH3HCO    5.000E+13      0.0         0.0
PC2H5O=SC2H5O                  1.000E+11      0.0     27000.0 
PC2H5O+PC2H5O=C2H5OH+CH3HCO    3.400E+13      0.0         0.0
C2H5OH=CH2OH+CH3               3.100E+15      0.0     80600.0 
C2H5OH+OH=SC2H5O+H2O           3.000E+13      0.0      5960.0 
C2H5OH+OH=C2H5O+H2O            1.138E+06      2.0       914.0 
C2H5OH+OH=PC2H5O+H2O           2.563E+06      2.06      860.0 
C2H5OH+O=SC2H5O+OH             6.000E+05      2.46     1850.0 
C2H5OH+O=C2H5O+OH              4.820E+13      0.0      6856.0 
C2H5OH+O=PC2H5O+OH             5.000E+12      0.0      4411.0 
C2H5OH+H=C2H5+H2O              5.900E+11      0.0      3450.0 
C2H5OH+H=SC2H5O+H2             4.400E+12      0.0      4570.0 
C2H5OH+HO2=SC2H5O+H2O2         2.000E+13      0.0     17000.0 
C2H5OH+CH3=SC2H5O+CH4          4.000E+11      0.0      9700.0 
C2H5OH+CH3=PC2H5O+CH4          3.000E+00      4.0     10480.0 
C2H5OH+CH3=C2H5O+CH4           8.000E+10      0.0      9400.0 
C2H5OH+CH3O=SC2H5O+CH3OH       2.000E+11      0.0      7000.0 
C2H5OH+CH2O=C2H5O+CH3O         1.500E+12      0.0     79500.0 
C2H5OH+C2H5O=C2H5OH+SC2H5O     2.000E+11      0.0      7000.0 
C2H5OH=C2H5+OH                 5.000E+16      0.0     91212.0 
C2H5OH=C2H4+H2O                1.000E+14      0.0     76706.0 
C2H5OH+O2=PC2H5O+HO2           4.000E+13      0.0     50900.0 
C2H5OH+O2=SC2H5O+HO2           4.000E+13      0.0     51200.0 
C2H5OH+O2=C2H5O+HO2            2.000E+13      0.0     56000.0 
C2H5OH+H=PC2H5O+H2             2.000E+12      0.0      9500.0 
C2H5OH+H=C2H5O+H2              1.760E+12      0.0      4570.0 
C2H5OH+HO2=H2O2+C2H5O          1.000E+11      0.0     15500.0 
C2H5OH+HO2=H2O2+PC2H5O         1.000E+11      0.0     12500.0 
C2H5OH+C2H5=PC2H5O+C2H6        1.500E+12      0.0     11700.0 
C2H5OH+C2H5=SC2H5O+C2H6        4.000E+13      0.0     10000.0 
C2H5OH+CH2OH=SC2H5O+CH3OH      4.000E+11      0.0      9700.0 
C+OH=CO+H                      5.000E+13      0.0         0.0 
C+O2=CO+O                      1.200E+14      0.0      4000.0 
C+CH3=C2H2+H                   5.000E+13      0.0         0.0 
C+CH2=C2H+H                    5.000E+13      0.0         0.0 
CH2O+CH3O2=HCO+CH3O2H          2.000E+12      0.0     11660.0 
CH3O2+H=CH3O+OH                9.600E+13      0.0         0.0
CH3O2+OH=CH3OH+O2              6.000E+13      0.0         0.0
CH3O2+CH3=CH3O+CH3O            2.400E+13      0.0         0.0 
CH3O2+CH3O2=>CH2O+CH3OH+O2     2.700E+10      0.0      -780.0 
CH3O2+CH3O2=>CH3O+CH3O+O2      2.800E+10      0.0      -780.0 
CH3O2+H2O2=CH3O2H+HO2          2.400E+12      0.0     10000.0 
CH3O2H=CH3O+OH                 6.000E+14      0.0     42300.0 
CH3O2+HO2=CH3O2H+O2            2.290E+11      0.0     -1550.0 
CH3O2H+OH=CH3O2+H2O            1.150E+12      0.0      -380.0 
CH4+CH3O2=CH3+CH3O2H           1.810E+11      0.0     18600.0 
C2H6+CH3O2=C2H5+CH3O2H         2.950E+11      0.0     14940.0
CH3OH+CH3O2=CH2OH+CH3O2H       1.810E+12      0.0     13800.0 
CH3O2H+O=OH+CH3O2              2.000E+13      0.0      4750.0 
CH3CO+O2=CH3CO3                1.000E+10      0.0     -2700.0 
CH3HCO+CH3CO3=CH3CO+CH3CO3H    1.200E+11      0.0      4900.0 
CH3HCO+C2H5O2=CH3CO+C2H5O2H    1.150E+11      0.0     10000.0 
C2H5+O2(+M)=C2H5O2(+M)         2.200E+10      0.772    -570.0 
     LOW /7.10E+42 -8.24 4270.0 /
C2H5O2=C2H4+HO2                5.620E+11      0.0     28900.0
C2H5O2+HO2=C2H5O2H+O2          3.400E+11      0.0     -1300.0 
C2H5O2H=C2H5O+OH               4.000E+15      0.0     43000.0 
C2H5O2H+O=OH+C2H5O2            2.000E+13      0.0      4750.0 
C2H5O2H+OH=C2H5O2+H2O          2.000E+12      0.0      -370.0 
CH4+C2H5O2=CH3+C2H5O2H         1.140E+13      0.0     20460.0 
CH4+CH3CO3=CH3+CH3CO3H         1.140E+13      0.0     20460.0 
C2H4+C2H5O2=C2H3+C2H5O2H       1.000E+12      0.0     25000.0 
C2H4+CH3CO3=C2H3+CH3CO3H       3.000E+12      0.0     29000.0 
CH3CO3+HO2=CH3CO3H+O2          1.000E+12      0.0         0.0 
CH3CO3H=>CH3CO2+OH             1.150E+13      0.0     32550.0 
CH3CO3H=>CH3+CO2+OH            2.000E+14      0.0     40150.0 
CH3CO3+CH3O2=>CH3CO2+CH3O+O2   1.080E+15      0.0      3600.0 
CH3CO3+CH3O2=>CH3CO2H+CH2O+O2  2.470E+09      0.0     -4200.0 
CH3CO3+HO2=>CH3CO2+OH+O2       2.590E+11      0.0     -2080.0 
CH3CO3+CH3CO3=>CH3CO2+CH3CO2+O2  1.690E+12    0.0     -1060.0 
CH3CO2+M=>CH3+CO2+M            8.700E+15      0.0     14400.0 
CH3CO2H=CH4+CO2                7.080E+13      0.0     74600.0
CH3CO2H=CH2CO+H2O              4.470E+14      0.0     79800.0
CH3CO2H+OH=CH3CO2+H2O          2.400E+11      0.0      -400.0
CH3OH+C2H5O2=CH2OH+C2H5O2H     6.300E+12      0.0     19360.0 
CH3OH+CH3CO3=CH2OH+CH3CO3H     6.300E+12      0.0     19360.0 
CH2O+C2H5O2=HCO+C2H5O2H        1.300E+11      0.0      9000.0 
CH2O+CH3CO3=HCO+CH3CO3H        1.000E+12      0.0     10560.0 
C2H4+CH3O2=C2H3+CH3O2H         1.000E+13      0.0     25000.0 
CH3HCO+CH3O2=CH3CO+CH3O2H      1.150E+11      0.0     10000.0 
C2H5OH+CH3O2=SC2H5O+CH3O2H     1.000E+13      0.0     10000.0 
C2H5+CH3O2=C2H5O+CH3O          2.410E+13      0.0         0.0 
C2H4+HO2=C2H4O+OH              2.200E+12      0.0     17200.0 
C2H4+CH3O=C2H4O+CH3            1.000E+11      0.0     14500.0 
C2H4+CH3O2=C2H4O+CH3O          7.000E+11      0.0     14500.0 
C2H4O=>CH3HCOW                 1.600E+13      0.0     54300.0 
CH3HCOW+M=>CH3HCO+M            1.000E+14      0.0         0.0 
CH3HCOW=>CH3+HCO               5.000E+08      0.0         0.0 
C2H4O+H=H2+C2H3O               8.000E+13      0.0      9740.0 
C2H4O+H=H2O+C2H3               5.000E+09      0.0      5030.0 
C2H4O+H=C2H4+OH                9.510E+10      0.0      5030.0 
C2H4O+CH2HCO=CH3HCO+C2H3O      1.000E+11      0.0     14000.0 
C2H4O+CH3=CH4+C2H3O            1.070E+12      0.0     11900.0 
C2H4O+O=OH+C2H3O               1.910E+12      0.0      5300.0 
C2H4O+OH=H2O+C2H3O             1.780E+13      0.0      3600.0 
C2H3O=>CH2CHOW                 1.000E+11      0.0     10000.0 
C2H3O=>CH3+CO                  8.000E+11      0.0     10000.0 
C2H3O+H+M=>C2H4O+M             4.000E+15      0.0         0.0 
CH2CHOW+M=>CH2HCO+M            1.000E+14      0.0         0.0 
CH2CHOW=>CH3+CO                1.000E+08      0.0         0.0 
CH2CHOW=>OH+C2H2               1.000E+11      0.0     17000.0 
CH2CHOW=>CH2CO+H               1.000E+08      0.0         0.0 
C2H4O+O2=HO2+C2H3O             1.000E+14      0.0     52000.0 
C2H4O+HO2=H2O2+C2H3O           5.000E+13      0.0     18000.0 
CH3HCOW+O2=>HO2+CH3CO          1.000E+14      0.0         0.0 
CH2CHOW+O2=>HO2+CH2CO          1.000E+14      0.0         0.0 
CH2+C2H2=H+C3H3                1.200E+13      0.0      6620.0 
CH2+C2H4=C3H6                  3.160E+12      0.0      5280.0 
SCH2+C2H4=>C3H6                1.000E+14      0.0         0.0
CH2+C3H8=CH3+IC3H7             1.500E+00      3.46     7470.0 
CH2+C3H8=CH3+NC3H7             9.000E-01      3.65     7150.0 
SCH2+C2H2=C3H3+H               1.800E+14      0.0         0.0 
C2H3+CH2=C3H4+H                3.000E+13      0.0         0.0 
C2H3+C2H2=C4H4+H               1.930E+12      0.0      6000.0 
C2H3+C2H3=C4H6                 7.230E+13      0.0         0.0 
C2H2+CH3=SC3H5                 1.610E+40     -8.58    20331.0 
C2H2+CH3=C3H5                  2.610E+46     -9.82    36951.0 
C2H2+CH3=C3H4+H                6.740E+19     -2.08    31591.0 
CH2CO+C2H3=C3H5+CO             1.000E+12      0.0      3000.0 
HCCO+C2H2=C3H3+CO              1.000E+11      0.0      3000.0 
C3H8(+M)=C2H5+CH3(+M)          1.100E+17      0.0     84400.0 
    LOW /7.83E+18 0.0 65000.0/
C3H8+O2=NC3H7+HO2              4.000E+13      0.0     50870.0 
C3H8+O2=IC3H7+HO2              4.000E+13      0.0     47690.0 
C3H8+HO2=NC3H7+H2O2            4.760E+04      2.55    16490.0 
C3H8+HO2=IC3H7+H2O2            9.640E+03      2.6     13910.0 
C3H8+OH=NC3H7+H2O              3.160E+07      1.80      934.0 
C3H8+OH=IC3H7+H2O              7.060E+06      1.90     -159.0 
C3H8+O=NC3H7+OH                3.715E+06      2.4      5505.0 
C3H8+O=IC3H7+OH                5.495E+05      2.5      3140.0 
C3H8+H=NC3H7+H2                1.336E+06      2.54     6756.0 
C3H8+H=IC3H7+H2                1.300E+06      2.4      4470.0 
C3H8+CH3=NC3H7+CH4             9.000E-01      3.65     7150.0 
C3H8+CH3=IC3H7+CH4             1.500E+00      3.46     5480.0 
C3H8+C2H5=NC3H7+C2H6           9.000E-01      3.65     9140.0 
C3H8+C2H5=IC3H7+C2H6           1.200E+00      3.46     7470.0 
C3H8+C2H3=NC3H7+C2H4           6.000E+02      3.3     10502.0 
C3H8+C2H3=IC3H7+C2H4           1.000E+03      3.1      8829.0 
C3H8+IC3H7=NC3H7+C3H8          8.440E-03      4.2      8720.0 
C3H8+C3H5=NC3H7+C3H6           2.350E+02      3.3     19800.0 
C3H8+C3H5=IC3H7+C3H6           7.840E+01      3.3     18200.0 
C3H8+CH3O=NC3H7+CH3OH          4.340E+11      0.0      6460.0 
C3H8+CH3O=IC3H7+CH3OH          1.450E+11      0.0      4570.0 
NC3H7=C2H4+CH3                 1.260E+13      0.0     30404.0 
NC3H7+O2=C3H6+HO2              1.000E+12      0.0      5000.0 
IC3H7=C2H4+CH3                 1.000E+12      0.0     34500.0 
IC3H7+O2=C3H6+HO2              2.754E+10      0.0     -2151.0 
C3H6=C3H5+H                    4.570E+14      0.0     88900.0 
C3H6=SC3H5+H                   7.590E+14      0.0    101300.0 
C3H6=TC3H5+H                   1.450E+15      0.0     98060.0 
C3H6=C2H3+CH3                  1.100E+21     -1.2     97720.0 
C3H6+HO2=C3H6O+OH              1.050E+12      0.0     14210.0 
C3H6+HO2=C3H5+H2O2             9.640E+03      2.6     13910.0 
C3H6+HO2=SC3H5+H2O2            7.500E+09      0.0     12570.0 
C3H6+HO2=TC3H5+H2O2            3.000E+09      0.0      9930.0 
C3H6+OH=C3H5+H2O               3.120E+06      2.0      -300.0 
C3H6+OH=SC3H5+H2O              2.140E+06      2.0      2780.0 
C3H6+OH=TC3H5+H2O              1.110E+06      2.0      1450.0 
C3H6+O=C2H5+HCO                6.833E+06      1.57     -628.0 
C3H6+O=CH3+CH3CO               9.111E+06      1.57     -628.0 
C3H6+O=C2H4+CH2O               4.555E+06      1.57     -628.0 
NC3H7=C3H6+H                   1.000E+14      0.0     37286.0 
C3H6+H=IC3H7                   5.704E+09      1.16      874.0 
C3H6+H=C3H5+H2                 6.457E+12      0.0      4445.0 
C3H6+H=SC3H5+H2                7.810E+05      2.5     12280.0 
C3H6+O2=SC3H5+HO2              1.950E+12      0.0     39000.0 
C3H6+O2=TC3H5+HO2              1.950E+12      0.0     39000.0 
C3H6+O2=C3H5+HO2               1.950E+12      0.0     39000.0 
C3H6+CH3=C3H5+CH4              2.210E+00      3.5      5680.0 
C3H6+CH3=SC3H5+CH4             1.350E+00      3.5     12850.0 
C3H6+CH3=TC3H5+CH4             8.400E-01      3.5     11660.0 
C3H6+C2H5=C3H5+C2H6            2.230E+00      3.5      6640.0 
C3H6O=C2H5+HCO                 2.450E+13      0.0     58500.0 
C3H6O=C2H5CHO                  1.820E+14      0.0     58500.0 
C3H6O=CH3+CH3CO                4.540E+13      0.0     59900.0 
C3H6O=CH3+CH2HCO               2.450E+13      0.0     58820.0 
C3H6O=CH3+C2H3O                8.000E+15      0.0     92010.0 
C2H5CHO=C2H5+HCO               2.450E+16      0.0     73000.0 
C2H5CHO+O=C2H5CO+OH            5.680E+12      0.0      1540.0 
C2H5CHO+OH=C2H5CO+H2O          1.210E+13      0.0         0.0 
C2H5CHO+HO2=C2H5CO+H2O2        1.520E+09      0.0         0.0 
C2H5CHO+C2H5=C2H5CO+C2H6       5.000E+10      0.0      6290.0 
C2H5CO=C2H5+CO                 5.890E+12      0.0     14400.0 
C3H5+O2=>CH2O+CH2HCO           5.000E+12      0.0     19190.0 
C3H5+H=C3H4+H2                 1.800E+13      0.0         0.0 
C3H5+O=>C2H4+CO+H              1.807E+14      0.0         0.0 
C3H5+CH3=C3H4+CH4              3.000E+12     -0.32     -130.0 
C3H5+C2H5=C3H4+C2H6            9.640E+11      0.0      -130.0 
C3H5+C2H3=C3H4+C2H4            2.400E+12      0.0         0.0 
C3H5+C2H3=C3H6+C2H2            4.800E+12      0.0         0.0
SC3H5+O2=CH3HCO+HCO            4.340E+12      0.0         0.0 
SC3H5+HO2=>CH2CO+CH3+OH        4.500E+12      0.0         0.0 
SC3H5+H=C3H4+H2                3.333E+12      0.0         0.0 
SC3H5+O=>CH2CO+CH3             1.807E+14      0.0         0.0 
SC3H5+CH3=C3H4+CH4             1.000E+11      0.0         0.0 
SC3H5+C2H5=C3H4+C2H6           1.000E+11      0.0         0.0 
SC3H5+C2H3=C3H4+C2H4           1.000E+11      0.0         0.0 
TC3H5+O2=CH3CO+CH2O            4.335E+11      0.0         0.0 
TC3H5+HO2=>CH2CO+CH3+OH        4.500E+12      0.0         0.0 
TC3H5+H=C3H4+H2                3.333E+12      0.0         0.0 
TC3H5+O=>HCCO+CH3+H            1.807E+14      0.0         0.0 
TC3H5+CH3=C3H4+CH4             1.000E+11      0.0         0.0 
TC3H5+C2H5=C3H4+C2H6           1.000E+11      0.0         0.0 
TC3H5+C2H3=C3H4+C2H4           1.000E+11      0.0         0.0 
C3H4+M=C3H3+H+M                2.000E+18      0.0     80000.0 
H2O/16.0/ CO2/3.75/ CO/1.875/ H2/2.5/ CH4/3.0/ C3H6/16.0/ C2H4/16.0/ C3H8/16.0/
C3H4(+M)=PC3H4(+M)             1.070E+14      0.0     64300.0 
    LOW /3.48E+17 0.0 48390.0/
C3H4+O2=C3H3+HO2               4.000E+13      0.0     61500.0 
C3H4+HO2=>CH2CO+CH2+OH         8.000E+12      0.0     19000.0 
C3H4+OH=CH2CO+CH3              3.120E+12      0.0      -397.0 
C3H4+OH=C3H3+H2O               2.000E+07      2.0      1000.0 
C3H4+O=C2H3+HCO                1.100E-02      4.613   -4243.0 
C3H4+H=C3H5                    1.200E+11      0.69     3000.0 
C3H4+H=TC3H5                   8.500E+12      0.0      2000.0 
C3H4+H=C3H3+H2                 2.000E+07      2.0      5000.0 
C3H4+CH3=C3H3+CH4              2.000E+11      0.0      7700.0 
PC3H4+M=C3H3+H+M               4.700E+18      0.0     80000.0 
H2O/16.0/ CO2/3.75/ CO/1.875/ H2/2.5/ CH4/3.0/ C3H6/16.0/ C2H4/16.0/ C3H8/16.0/
PC3H4+O2=>HCCO+OH+CH2          2.000E+08      1.5     30100.0 
PC3H4+O2=C3H3+HO2              5.000E+12      0.0     51000.0 
PC3H4+HO2=>C2H4+CO+OH          3.000E+12      0.0     19000.0 
PC3H4+OH=C3H3+H2O              2.000E+07      2.0      1000.0 
PC3H4+OH=CH2CO+CH3             5.000E-04      4.5     -1000.0 
PC3H4+O=CH2CO+CH2              6.400E+12      0.0      2010.0 
PC3H4+O=C2H3+HCO               3.200E+12      0.0      2010.0 
PC3H4+O=HCCO+CH3               6.300E+12      0.0      2010.0 
PC3H4+O=>HCCO+CH2+H            3.200E+11      0.0      2010.0 
PC3H4+H=TC3H5                  6.500E+12      0.0      2000.0 
PC3H4+H=C3H3+H2                2.000E+07      2.0      5000.0 
PC3H4+H=C2H2+CH3               1.300E+05      2.5      1000.0 
PC3H4+CH3=C3H3+CH4             1.500E+00      3.5      5600.0 
PC3H4+C2H3=C3H3+C2H4           1.000E+12      0.0      7700.0 
PC3H4+C3H5=C3H3+C3H6           1.000E+12      0.0      7700.0 
C3H3+H=C3H2+H2                 5.000E+13      0.0      3000.0 
C3H3+O=>C2H+HCO+H              7.000E+13      0.0         0.0 
C3H3+O=>C2H2+CO+H              7.000E+13      0.0         0.0 
C3H3+OH=C3H2+H2O               1.000E+13      0.0         0.0 
C3H3+O2=CH2CO+HCO              3.010E+10      0.0      2870.0 
C3H3+CH=IC4H3+H                7.000E+13      0.0         0.0 
C3H3+CH=NC4H3+H                7.000E+13      0.0         0.0 
C3H3+CH2=C4H4+H                4.000E+13      0.0         0.0 
C3H3+C3H3=C6H5+H               2.000E+12      0.0         0.0 
CH+C2H2=C3H2+H                 1.000E+14      0.0         0.0 
C3H2+O2=HCCO+CO+H              1.000E+14      0.0      3000.0 
C3H2+OH=C2H2+HCO               5.000E+13      0.0         0.0 
C3H2+CH2=IC4H3+H               3.000E+13      0.0         0.0 
C4H8=IC4H7+H                   4.078E+18     -1.0     97350.0 
C4H8=C2C4H8                    4.000E+11      0.0     60000.0 
C4H8=T2C4H8                    4.000E+11      0.0     60000.0 
C4H8=C3H5+CH3                  1.000E+16      0.0     73000.0 
C4H8=C2H3+C2H5                 1.000E+19     -1.0     96770.0 
C4H8+O2=IC4H7+HO2              4.000E+12      0.0     33200.0 
C4H8+HO2=IC4H7+H2O2            1.000E+11      0.0     17060.0 
C4H8+OH=NC3H7+CH2O             6.500E+12      0.0         0.0 
C4H8+OH=CH3HCO+C2H5            1.000E+11      0.0         0.0 
C4H8+OH=C2H6+CH3CO             1.000E+10      0.0         0.0 
C4H8+OH=IC4H7+H2O              2.250E+13      0.0      2217.0 
C4H8+O=C3H6+CH2O               2.505E+12      0.0         0.0 
C4H8+O=CH3HCO+C2H4             1.250E+12      0.0       850.0 
C4H8+O=C2H5+CH3CO              1.625E+13      0.0       850.0 
C4H8+O=IC4H7+OH                9.600E+12      0.0      1970.0 
C4H8+O=NC3H7+HCO               1.800E+05      2.5     -1029.0
C4H8+H=IC4H7+H2                5.000E+13      0.0      3900.0 
C4H8+CH3=IC4H7+CH4             1.000E+11      0.0      7300.0 
C4H8+C2H5=IC4H7+C2H6           1.000E+11      0.0      8000.0 
C4H8+C3H5=IC4H7+C3H6           7.900E+10      0.0     12400.0 
C4H8+SC3H5=IC4H7+C3H6          8.000E+10      0.0     12400.0 
C4H8+TC3H5=IC4H7+C3H6          8.000E+10      0.0     12400.0 
C2C4H8=T2C4H8                  4.000E+13      0.0     62000.0 
C2C4H8=C4H6+H2                 1.000E+13      0.0     65500.0 
C2C4H8=IC4H7+H                 4.074E+18     -1.0     97350.0 
C2C4H8=SC3H5+CH3               2.000E+16      0.0     95000.0 
C2C4H8+OH=IC4H7+H2O            1.250E+14      0.0      3060.0 
C2C4H8+OH=CH3HCO+C2H5          1.400E+13      0.0         0.0 
C2C4H8+O=IC3H7+HCO             6.030E+12      0.0         0.0 
C2C4H8+O=CH3HCO+C2H4           1.000E+12      0.0         0.0 
C2C4H8+H=IC4H7+H2              1.000E+13      0.0      3500.0 
C2C4H8+CH3=IC4H7+CH4           1.000E+11      0.0      8200.0 
T2C4H8=IC4H7+H                 4.074E+18     -1.0     97350.0 
T2C4H8=SC3H5+CH3               2.000E+16      0.0     96000.0 
T2C4H8+OH=IC4H7+H2O            1.000E+14      0.0      3060.0 
T2C4H8+OH=CH3HCO+C2H5          1.500E+13      0.0         0.0 
T2C4H8+O=IC3H7+HCO             6.030E+12      0.0         0.0 
T2C4H8+O=CH3HCO+C2H4           1.000E+12      0.0         0.0 
T2C4H8+H=IC4H7+H2              5.000E+12      0.0      3500.0 
T2C4H8+CH3=IC4H7+CH4           1.000E+11      0.0      8200.0 
IC4H7=C4H6+H                   1.200E+14      0.0     49300.0 
IC4H7=C2H4+C2H3                1.000E+14      0.0     49000.0 
IC4H7+H=C4H6+H2                3.160E+12      0.0         0.0 
IC4H7+O2=C4H6+HO2              1.000E+11      0.0         0.0 
IC4H7+CH3=C4H6+CH4             1.000E+13      0.0         0.0 
IC4H7+C2H3=C4H6+C2H4           4.000E+12      0.0         0.0 
IC4H7+C2H5=C4H6+C2H6           4.000E+12      0.0         0.0 
IC4H7+C2H5=C4H8+C2H4           5.000E+11      0.0         0.0 
IC4H7+C2H5=T2C4H8+C2H4         5.000E+11      0.0         0.0 
IC4H7+C2H5=C2C4H8+C2H4         5.000E+11      0.0         0.0 
IC4H7+C3H5=C4H6+C3H6           4.000E+13      0.0         0.0 
IC4H7+IC4H7=C4H6+C4H8          3.160E+12      0.0         0.0 
C2H3+C2H4=C4H6+H               3.000E+12      0.0      1000.0 
C4H6+H=NC4H5+H2                3.000E+07      2.0     13000.0 
C4H6+H=IC4H5+H2                3.000E+07      2.0      6000.0 
C4H6+OH=NC4H5+H2O              2.000E+07      2.0      5000.0 
C4H6+OH=IC4H5+H2O              2.000E+07      2.0      2000.0 
C4H6+O=C2H4+CH2CO              1.000E+12      0.0         0.0 
C4H6+O=PC3H4+CH2O              1.000E+12      0.0         0.0 
C2H2+NC4H5=C6H6+H              2.800E+03      2.9      1400.0 
NC4H5+OH=C4H4+H2O              2.000E+07      2.0      1000.0 
NC4H5+H=C4H4+H2                3.000E+07      2.0      1000.0 
NC4H5+H=IC4H5+H                1.000E+14      0.0         0.0 
IC4H5=C4H4+H                   2.000E+15      0.0     45000.0 
NC4H5=C4H4+H                   1.600E+14      0.0     41400.0 
C4H4+OH=IC4H3+H2O              1.000E+07      2.0      2000.0 
C4H4+OH=NC4H3+H2O              7.500E+06      2.0      5000.0 
C4H4+H=NC4H3+H2                2.000E+07      2.0     15000.0 
NC4H3+H=IC4H3+H                1.000E+14      0.0         0.0 
IC4H3+CH2=C3H4+C2H             2.000E+13      0.0         0.0 
IC4H3+O2=CH2CO+HCCO            1.000E+12      0.0         0.0 
IC4H3+OH=C4H2+H2O              3.000E+13      0.0         0.0 
IC4H3+O=CH2CO+C2H              2.000E+13      0.0         0.0 
IC4H3+H=C4H2+H2                5.000E+13      0.0         0.0 
NC4H3+C2H2=C6H5                2.800E+03      2.9      1400.0 
NC4H3+M=C4H2+H+M               1.000E+16      0.0     59700.0 
H2O/16.0/ CO2/3.75/ CO/1.875/ H2/2.5/ CH4/3.0/ C3H6/16.0/ C2H4/16.0/ C3H8/16.0/
IC4H3+M=C4H2+H+M               4.460E+15      0.0     46516.0 
H2O/16.0/ CO2/3.75/ CO/1.875/ H2/2.5/ CH4/3.0/ C3H6/16.0/ C2H4/16.0/ C3H8/16.0/
IC4H3+O=H2C4O+H                2.000E+13      0.0         0.0 
H2C4O+H=C2H2+HCCO              5.000E+13      0.0      3000.0 
H2C4O+OH=CH2CO+HCCO            1.000E+07      2.0      2000.0 
C4H2+OH=H2C4O+H                6.660E+12      0.0      -410.0 
C2H2+C2H2=IC4H3+H              2.200E+12      0.0     64060.0 
C2H2+C2H2=NC4H3+H              1.000E+12      0.0     66000.0 
C2H2+C2H2=C4H4                 5.500E+12      0.0     37000.0
C4H2(+M)=C4H+H(+M)             2.200E+14      0.0    116740.0 
       LOW /3.50E+17  0.0  80065.0/
      H2O /16.0/ H2 /2.5/ CO /1.875/ CO2 /3.75/ CH4 /3.0/
      C3H6 /16.0/ C2H4 /16.0/ C3H8 /16.0/
C4H2+O=C3H2+CO                 2.700E+13      0.0      1720.0 
C2H2+C2H=C4H2+H                1.820E+14      0.0       467.0 
C2H2+C2H=NC4H3                 1.000E+13      0.0         0.0
C4H+O2=C2H+CO+CO               1.000E+14      0.0         0.0 
C2O+H=CH+CO                    1.320E+13      0.0         0.0 
C2O+O=CO+CO                    5.200E+13      0.0         0.0 
C2O+OH=CO+CO+H                 2.000E+13      0.0         0.0 
C2O+O2=CO+CO+O                 2.000E+13      0.0         0.0 
C2O+O2=CO+CO2                  2.000E+13      0.0         0.0
C2+H2=C2H+H                    6.600E+13      0.0      7950.0 
C2+O=C+CO                      3.600E+14      0.0         0.0
C2+O2=CO+CO                    9.000E+12      0.0       980.0 
C2+OH=C2O+H                    5.000E+13      0.0         0.0 
C6H5+OH=C6H5O+H                5.000E+13      0.0         0.0 
C6H5+O2=C6H5O+O                2.600E+13      0.0      6120.0 
C6H5+HO2=C6H5O+OH              5.000E+13      0.0      1000.0 
C6H6+H=C6H5+H2                 3.000E+12      0.0      8100.0 
C6H6+OH=C6H5+H2O               1.680E+08      1.42     1450.0 
C6H6+O=C6H5O+H                 2.780E+13      0.0      4910.0 
C6H6+O2=C6H5O+OH               4.000E+13      0.0     34000.0 
H+C6H5=C6H6                    7.800E+13      0.0         0.0 
C3H3+O=>C2H3+CO                3.800E+13      0.0         0.0 
C3H3+O=CH2O+C2H                2.000E+13      0.0         0.0 
C3H3+O2=>HCCO+CH2O             6.000E+12      0.0         0.0 
C3H3+CH3=C2H5+C2H              1.000E+13      0.0     37500.0 
C3H3+CH3=C4H6                  5.000E+12      0.0         0.0 
C3H6+C2H3=C3H5+C2H4            2.210E+00      3.5      4680.0 
C3H6+C2H3=SC3H5+C2H4           1.350E+00      3.5     10860.0 
C3H6+C2H3=TC3H5+C2H4           8.400E-01      3.5      9670.0 
C3H6+CH3O=C3H5+CH3OH           9.000E+01      2.95    12000.0 
CH2+C2H2=C3H4                  1.200E+13      0.0      6620.0 
C3H4+C3H4=C3H5+C3H3            5.000E+14      0.0     64700.0 
C3H4+OH=CH2O+C2H3              1.700E+12      0.0      -300.0 
C3H4+OH=HCO+C2H4               1.700E+12      0.0      -300.0 
C3H4+O=CH2O+C2H2               1.000E+12      0.0         0.0 
C3H4+O=>CO+C2H4                7.800E+12      0.0      1600.0 
C3H4+C3H5=C3H3+C3H6            2.000E+12      0.0      7700.0 
C3H4+C2H=C3H3+C2H2             1.000E+13      0.0         0.0 
PC3H4=C2H+CH3                  4.200E+16      0.0    100000.0 
PC3H4+C2H=C3H3+C2H2            1.000E+13      0.0         0.0 
C3H2+O2=HCO+HCCO               1.000E+13      0.0         0.0 
C2H2+C2H3=NC4H5                2.510E+05      1.9      2100.0 
C2H3+C2H3=IC4H5+H              4.000E+13      0.0         0.0 
IC4H5+H=C4H4+H2                3.000E+07      2.0      1000.0 
C4H2+H=C4H+H2                  1.000E+14      0.0     35000.0 
C4H6+OH=C3H5+CH2O              7.230E+12      0.0      -994.0 
C4H8+IC4H7=IC4H7+C2C4H8        3.980E+10      0.0     12400.0 
C4H8+IC4H7=IC4H7+T2C4H8        3.980E+10      0.0     12400.0 
C3H3+C3H3=C6H6                 3.000E+11      0.0         0.0 
C3H3+C3H4=C6H6+H               1.400E+12      0.0     10000.0 
C3H5+C2H5=C3H6+C2H4            2.600E+12      0.0      -130.0 
C3H6+OH=C2H5+CH2O              8.000E+12      0.0         0.0 
C3H6+OH=CH3+CH3HCO             3.400E+11      0.0         0.0 
C3H5+O2=C3H4+HO2               1.200E+12      0.0     13550.0 
CH2O+C3H5=HCO+C3H6             8.000E+10      0.0     12400.0 
CH3HCO+C3H5=CH3CO+C3H6         3.800E+11      0.0      7200.0 
C3H8+CH3O2=NC3H7+CH3O2H        6.030E+12      0.0     19380.0 
C3H8+CH3O2=IC3H7+CH3O2H        1.990E+12      0.0     17050.0 
C3H8+C2H5O2=NC3H7+C2H5O2H      6.030E+12      0.0     19380.0 
C3H8+C2H5O2=IC3H7+C2H5O2H      1.990E+12      0.0     17050.0 
C3H8+IC3H7O2=NC3H7+IC3H7O2H    6.030E+12      0.0     19380.0 
C3H8+IC3H7O2=IC3H7+IC3H7O2H    1.990E+12      0.0     17050.0 
C3H8+NC3H7O2=NC3H7+NC3H7O2H    6.030E+12      0.0     19380.0 
C3H8+NC3H7O2=IC3H7+NC3H7O2H    1.990E+12      0.0     17050.0 
NC3H7+O2=NC3H7O2               4.820E+12      0.0         0.0 
IC3H7+O2=IC3H7O2               6.620E+12      0.0         0.0 
NC3H7+HO2=NC3H7O+OH            3.200E+13      0.0         0.0 
IC3H7+HO2=IC3H7O+OH            3.200E+13      0.0         0.0 
NC3H7+CH3O2=NC3H7O+CH3O        3.800E+12      0.0     -1200.0 
IC3H7+CH3O2=IC3H7O+CH3O        3.800E+12      0.0     -1200.0 
NC3H7+NC3H7O2=NC3H7O+NC3H7O    3.800E+12      0.0     -1200.0 
IC3H7+NC3H7O2=IC3H7O+NC3H7O    3.800E+12      0.0     -1200.0 
NC3H7+IC3H7O2=NC3H7O+IC3H7O    3.800E+12      0.0     -1200.0 
IC3H7+IC3H7O2=IC3H7O+IC3H7O    3.800E+12      0.0     -1200.0 
NC3H7O2+HO2=NC3H7O2H+O2        4.600E+10      0.0     -2600.0 
IC3H7O2+HO2=IC3H7O2H+O2        4.600E+10      0.0     -2600.0 
CH3+NC3H7O2=CH3O+NC3H7O        3.800E+12      0.0     -1200.0 
CH3+IC3H7O2=CH3O+IC3H7O        3.800E+12      0.0     -1200.0 
NC3H7O2H=NC3H7O+OH             4.000E+15      0.0     43000.0 
IC3H7O2H=IC3H7O+OH             4.000E+15      0.0     43000.0 
NC3H7O=C2H5+CH2O               5.000E+13      0.0     15700.0 
IC3H7O=CH3+CH3HCO              4.000E+14      0.0     17200.0 
C3H6+OH(+M)=C3H6OH(+M)         1.810E+13      0.0         0.0 
    LOW /1.33E+30  -3.5  0.0/
C3H6OH=>C2H5+CH2O              1.400E+09      0.0     17200.0 
C3H6OH=>CH3+CH3HCO             1.000E+09      0.0     17200.0 
C3H6OH+O2=O2C3H6OH             1.000E+12      0.0     -1100.0 
O2C3H6OH=>CH3HCO+CH2O+OH       1.000E+16      0.0     25000.0 
C3H6+CH3O2=C3H5+CH3O2H         2.000E+12      0.0     17000.0 
C3H6+CH3O2=C3H6O+CH3O          4.000E+11      0.0     11720.0 
C3H6+C2H5O2=C3H5+C2H5O2H       3.200E+11      0.0     14900.0 
C3H6+C3H5O2=C3H5+C3H5O2H       3.200E+11      0.0     14900.0 
C3H6+C3H5O2=C3H6O+C3H5O        1.050E+11      0.0     14200.0 
C3H6+CH3CO3=C3H5+CH3CO3H       3.200E+11      0.0     14900.0 
C3H6+NC3H7O2=C3H5+NC3H7O2H     3.200E+11      0.0     14900.0 
C3H6+IC3H7O2=C3H5+IC3H7O2H     3.200E+11      0.0     14900.0 
C3H6+NC3H7O2=C3H6O+NC3H7O      1.700E+07      0.0         0.0 
C3H5+O2=C3H5O2                 1.200E+10      0.0     -2300.0 
C3H5+HO2=C3H5O+OH              9.000E+12      0.0         0.0 
C3H5+CH3O2=C3H5O+CH3O          3.800E+11      0.0     -1200.0 
C3H5O2+CH3=C3H5O+CH3O          3.800E+11      0.0     -1200.0 
C3H5O2+C3H5=C3H5O+C3H5O        3.800E+11      0.0     -1200.0 
C3H5O2+HO2=C3H5O2H+O2          4.600E+10      0.0     -2600.0 
C3H5O2+HO2=>C3H5O+OH+O2        1.000E+12      0.0         0.0 
C3H5O2+CH3O2=>C3H5O+CH3O+O2    1.700E+11      0.0     -1000.0 
C3H5O2+C3H5O2=>C3H5O+C3H5O+O2  3.700E+12      0.0      2200.0 
C3H5O=CH2O+C2H3                1.000E+14      0.0     21600.0 
C3H5O2H=C3H5O+OH               4.000E+15      0.0     43000.0 
CH2O+C3H5O2=HCO+C3H5O2H        1.300E+11      0.0     10500.0 
CH2O+NC3H7O2=HCO+NC3H7O2H      1.300E+11      0.0      9000.0 
CH2O+IC3H7O2=HCO+IC3H7O2H      1.300E+11      0.0      9000.0 
C2H4+NC3H7O2=C2H3+NC3H7O2H     7.100E+11      0.0     25000.0 
C2H4+IC3H7O2=C2H3+IC3H7O2H     7.100E+11      0.0     25000.0 
CH4+C3H5O2=CH3+C3H5O2H         1.140E+13      0.0     20460.0 
CH4+NC3H7O2=CH3+NC3H7O2H       1.140E+13      0.0     20460.0 
CH4+IC3H7O2=CH3+IC3H7O2H       1.140E+13      0.0     20460.0 
CH3OH+NC3H7O2=CH2OH+NC3H7O2H   6.300E+12      0.0     19360.0 
CH3OH+IC3H7O2=CH2OH+IC3H7O2H   6.300E+12      0.0     19360.0 
CH3HCO+C3H5O2=CH3CO+C3H5O2H    1.150E+11      0.0     10000.0 
CH3HCO+NC3H7O2=CH3CO+NC3H7O2H  1.150E+11      0.0     10000.0 
CH3HCO+IC3H7O2=CH3CO+IC3H7O2H  1.150E+11      0.0     10000.0 

C+N2+M=CNN+M                   1.120E+15      0.0         0.0
C2H+NO=HCN+CO                  6.000E+13      0.0       570.0
C2H+HCN=CN+C2H2                3.200E+12      0.0      1530.0 
CH2+NO=HCN+OH                  5.000E+11      0.0      2870.0
HCN+M=H+CN+M                   3.570E+26     -2.6    124900.0 
C2N2+M=CN+CN+M                 3.200E+16      0.0     94400.0
CN+N2O=CNN+NO                  6.000E+13      0.0     15360.0
     DUPLICATE
CN+N2O=CNN+NO                  1.800E+10      0.0      1450.0
     DUPLICATE
CH+N2(+M)=HCNN(+M)             3.100E+12      0.15        0.0
     LOW / 1.30E+25 -3.16 740.0 /
     TROE /0.667 235.0 2117.0 4536.0 /
     H2O/10.0/ O2/2.0/ AR/0.75/ H2/2.0/  
HCNN+H=H2+CNN                  5.000E+13      0.0         0.0
HCNN+H=>CH2+N2                 2.000E+13      0.0      3000.0
HCNN+O=OH+CNN                  2.000E+13      0.0     20000.0
HCNN+O=CO+H+N2                 5.000E+13      0.0     15000.0
HCNN+O=HCN+NO                  5.000E+13      0.0     15000.0
HCNN+OH=H2O+CNN                1.000E+13      0.0      8000.0
HCNN+OH=H+HCO+N2               1.000E+13      0.0     16000.0
HCNN+O2=HO2+CNN                1.000E+12      0.0      4000.0
HCNN+O2=>H+CO2+N2              4.000E+12      0.0         0.0
HCNN+O2=HCO+N2O                4.000E+12      0.0         0.0  
CNN+O=CO+N2                    1.000E+13      0.0         0.0
CNN+O=CN+NO                    1.000E+14      0.0     20000.0
CNN+OH=H+CO+N2                 1.000E+13      0.0      1000.0
CNN+H=NH+CN                    5.000E+14      0.0     40000.0
CNN+OH=HCN+NO                  1.000E+12      0.0      1000.0
CNN+H=HCN+N                    5.000E+13      0.0     25000.0
CNN+O2=NO+NCO                  1.000E+13      0.0      5000.0
HNO+CH3=NO+CH4                 8.200E+05      1.87      954.0 
HONO+CH3=NO2+CH4               8.100E+05      1.87     5504.0
H2NO+CH3=CH3O+NH2              2.000E+13      0.0         0.0 
H2NO+CH3=HNO+CH4               1.600E+06      1.87     2960.0
HNOH+CH3=HNO+CH4               1.600E+06      1.87     2096.0
NH2OH+CH3=HNOH+CH4             1.600E+06      1.87     6350.0
NH2OH+CH3=H2NO+CH4             8.200E+05      1.87     5500.0
N2H2+CH3=NNH+CH4               1.600E+06      1.87     2970.0
N2H3+CH3=N2H2+CH4              8.200E+05      1.87     1818.0
N2H4+CH3=N2H3+CH4              3.300E+06      1.87     5325.0
CH4+NH=CH3+NH2                 9.000E+13      0.0     20080.0 
CH4+NH2=CH3+NH3                1.200E+13      0.0     15150.0 
CH3+NH2=CH2+NH3                1.600E+06      1.87     7570.0
C2H6+NH=C2H5+NH2               7.000E+13      0.0     16700.0 
C2H6+NH2=C2H5+NH3              9.700E+12      0.0     11470.0 
C3H8+NH2=NC3H7+NH3             1.700E+13      0.0     10660.0 
C3H8+NH2=IC3H7+NH3             4.500E+11      0.0      6150.0 
CH3+NO(+M)=CH3NO(+M)           1.000E+13      0.0         0.0
    LOW /1.90E+18 0.0 0.0/
    SRI /0.03 -790.0 1.0/
CH3NO+H=H2CNO+H2               4.400E+08      1.5       377.0 
CH3NO+H=CH3+HNO                1.800E+13      0.0      2800.0
CH3NO+O=H2CNO+OH               3.300E+08      1.5      3615.0
CH3NO+O=CH3+NO2                1.700E+06      2.08        0.0
CH3NO+OH=H2CNO+H2O             3.600E+06      2.0     -1192.0
CH3NO+OH=CH3+HONO              2.500E+12      0.0      1000.0
CH3NO+CH3=H2CNO+CH4            7.900E+05      1.87     5415.0
CH3NO+NH2=H2CNO+NH3            2.800E+06      1.94     1073.0
H2CNO=HNCO+H                   2.300E+42     -9.11    53840.0
H2CNO+O2=CH2O+NO2              2.900E+12     -0.31    17700.0
H2CNO+H=CH3+NO                 4.000E+13      0.0         0.0
H2CNO+H=HCNO+H2                4.800E+08      1.5      -894.0
H2CNO+O=HCNO+OH                3.300E+08      1.5      -894.0
H2CNO+O=CH2O+NO                7.000E+13      0.0         0.0
H2CNO+OH=CH2OH+NO              4.000E+13      0.0         0.0
H2CNO+OH=HCNO+H2O              2.400E+06      2.0     -1192.0
H2CNO+CH3=C2H5+NO              3.000E+13      0.0         0.0
H2CNO+CH3=HCNO+CH4             1.600E+06      1.87    -1113.0
H2CNO+NH2=HCNO+NH3             1.800E+06      1.94    -1152.0
CH3+NO2=CH3O+NO                1.400E+13      0.0         0.0
CH+NO2=HCO+NO                  1.200E+14      0.0         0.0
CH2+NO2=CH2O+NO                4.200E+13      0.0         0.0
CN+NO=N2+CO                    1.000E+11      0.0         0.0
HNCO+M=H+NCO+M                 5.000E+15      0.0    120000.0
HNCO+N=NH+NCO                  4.000E+13      0.0     36000.0
CH3O+HNO=CH3OH+NO              3.160E+13      0.0         0.0
NCO+HO2=HNCO+O2                2.000E+13      0.0         0.0
N2O+CO=CO2+N2                  2.510E+14      0.0     46000.0
N2O+CH2=CH2O+N2                1.000E+12      0.0         0.0
N2O+CH3=CH3O+N2                9.000E+09      0.0         0.0
N2O+HCO=CO2+H+N2               1.700E+14      0.0     20000.0
N2O+HCCO=CO+HCO+N2             1.700E+14      0.0     25500.0
N2O+C2H2=HCCO+H+N2             6.590E+16      0.0     61200.0
N2O+C2H3=CH2HCO+N2             1.000E+11      0.0         0.0
HOCN+O=NCO+OH                  1.500E+04      2.64     4000.0
HOCN+H=NCO+H2                  2.000E+07      2.0      2000.0
HOCN+H=NH2+CO                  1.200E+08      0.61     2080.0 
HOCN+OH=NCO+H2O                6.380E+05      2.0      2560.0
HOCN+CH3=NCO+CH4               8.200E+05      1.87     6620.0
HOCN+NH2=NCO+NH3               9.200E+05      1.94     3645.0
CN+NO2=CO+N2O                  4.930E+14     -0.752     344.0
CN+NO2=CO2+N2                  3.700E+14     -0.752     344.0
CN+CO2=NCO+CO                  3.670E+06      2.16    26900.0
CN+NH3=HCN+NH2                 9.200E+12      0.0      -357.0
HNCO+CN=HCN+NCO                1.500E+13      0.0         0.0
NCO+CN=CNN+CO                  1.800E+13      0.0         0.0
HONO+NCO=HNCO+NO2              3.600E+12      0.0         0.0
NCO+CH2O=HNCO+HCO              6.000E+12      0.0         0.0
CH+N2=HCN+N                    3.680E+07      1.42    20723.0 
NH2+C=CH+NH                    5.800E+11      0.67    20900.0 
C+N2=CN+N                      5.200E+13      0.0     44700.0 
CH2+N2=HCN+NH                  4.800E+12      0.0     35850.0 
C2+N2=CN+CN                    1.500E+13      0.0     41700.0
H2CN+N=N2+CH2                  6.000E+13      0.0       400.0 
H2CN+H=HCN+H2                  2.400E+08      1.5      -894.0
H2CN+O=HCN+OH                  1.700E+08      1.5      -894.0
H2CN+O=HNCO+H                  6.000E+13      0.0         0.0 
H2CN+O=HCNO+H                  2.000E+13      0.0         0.0
H2CN+M=HCN+H+M                 3.000E+14      0.0     22000.0 
H2CN+HO2=HCN+H2O2              1.400E+04      2.69    -1610.0
H2CN+O2=CH2O+NO                3.000E+12      0.0      6000.0
H2CN+CH3=HCN+CH4               8.100E+05      1.87    -1113.0
H2CN+OH=HCN+H2O                1.200E+06      2.0     -1192.0
H2CN+NH2=HCN+NH3               9.200E+05      1.94    -1152.0
C+NO=CN+O                      2.000E+13      0.0         0.0 
CH+NO=HCN+O                    8.690E+13      0.0         0.0 
CH+NO=CN+OH                    1.680E+12      0.0         0.0 
CH+NO=CO+NH                    9.840E+12      0.0         0.0 
CH+NO=NCO+H                    1.670E+13      0.0         0.0 
CH2+NO=HNCO+H                  2.500E+12      0.0      5970.0 
CH2+NO=HCNO+H                  3.800E+13     -0.36      576.0
CH2+NO=NH2+CO                  2.300E+16     -1.43     1331.0
CH2+NO=H2CN+O                  8.100E+07      1.42     4110.0 
CH3+NO=HCN+H2O                 2.400E+12      0.0     15700.0 
CH3+NO=H2CN+OH                 5.200E+12      0.0     24240.0 
HCCO+NO=HCNO+CO                4.640E+13      0.0       700.0 
HCCO+NO=HCN+CO2                1.390E+13      0.0       700.0
SCH2+NO=HCN+OH                 1.000E+14      0.0         0.0 
HCNO=HCN+O                     4.200E+31     -6.12    61210.0
HCNO+H=HCN+OH                  1.000E+14      0.0     12000.0 
HCNO+H=HNCO+H                  2.100E+15     -0.69     2850.0
HCNO+H=HOCN+H                  1.400E+11     -0.19     2484.0
HCNO+H=NH2+CO                  1.700E+14     -0.75     2890.0
HCNO+O=HCO+NO                  7.000E+13      0.0         0.0
CH2+N=HCN+H                    5.000E+13      0.0         0.0 
CH2+N=NH+CH                    6.000E+11      0.67    40500.0 
CH+N=CN+H                      1.670E+14     -0.09        0.0 
CH+N=C+NH                      4.500E+11      0.65     2400.0 
N+CO2=NO+CO                    1.900E+11      0.0      3400.0 
N+HCCO=HCN+CO                  5.000E+13      0.0         0.0 
CH3+N=H2CN+H                   7.100E+13      0.0         0.0 
CH3+N=HCNH+H                   1.200E+11      0.52      367.6
HCNH=HCN+H                     6.100E+28     -5.69    24270.0
HCNH+H=H2CN+H                  2.000E+13      0.0         0.0
HCNH+H=HCN+H2                  2.400E+08      1.5      -894.0
HCNH+O=HNCO+H                  7.000E+13      0.0         0.0
HCNH+O=HCN+OH                  1.700E+08      1.5      -894.0 
HCNH+OH=HCN+H2O                1.200E+06      2.0     -1192.0
HCNH+CH3=HCN+CH4               8.200E+05      1.87    -1113.0 
C2H3+N=HCN+CH2                 2.000E+13      0.0         0.0 
CN+H2O=HCN+OH                  4.000E+12      0.0      7400.0 
CN+H2O=HOCN+H                  4.000E+12      0.0      7400.0 
OH+HCN=HOCN+H                  3.200E+04      2.45    12120.0 
OH+HCN=HNCO+H                  5.600E-06      4.71     -490.0 
OH+HCN=NH2+CO                  6.440E+10      0.0     11700.0 
HOCN+H=HNCO+H                  1.000E+13      0.0         0.0 
HCN+O=NCO+H                    1.380E+04      2.64     4980.0 
HCN+O=NH+CO                    3.450E+03      2.64     4980.0 
HCN+O=CN+OH                    2.700E+09      1.58    26600.0 
CN+H2=HCN+H                    2.000E+04      2.87     1600.0 
CN+O=CO+N                      1.900E+12      0.46      720.0 
CN+O2=NCO+O                    7.200E+12      0.0      -400.0 
CN+OH=NCO+H                    4.000E+13      0.0         0.0 
CN+HCN=C2N2+H                  1.510E+07      1.71     1530.0 
CN+NO2=NCO+NO                  5.320E+15     -0.752     344.0 
CN+N2O=NCO+N2                  6.000E+12      0.0     15360.0 
C2N2+O=NCO+CN                  4.570E+12      0.0      8880.0 
C2N2+OH=HNCO+CN                1.860E+11      0.0      2900.0 
C2N2+OH=HOCN+CN                2.000E+12      0.0     19000.0
HNCO+H=H2+NCO                  1.760E+05      2.41    12300.0 
HNCO+H=NH2+CO                  3.600E+04      2.49     2340.0 
HNCO+M=NH+CO+M                 1.100E+16      0.0     86000.0 
  N2/1.5/ O2/1.5/ H2O/18.6/
HNCO+O=NCO+OH                  2.200E+06      2.11    11430.0 
HNCO+O=NH+CO2                  9.800E+07      1.41     8530.0 
HNCO+O=HNO+CO                  1.500E+08      1.57    44012.0 
HNCO+OH=NCO+H2O                3.450E+07      1.5      3600.0 
HNCO+OH=NH2+CO2                6.300E+10     -0.06    11645.0
HNCO+HO2=NCO+H2O2              3.000E+11      0.0     29000.0 
HNCO+O2=HNO+CO2                1.000E+12      0.0     35000.0 
HNCO+NH2=NCO+NH3               5.000E+12      0.0      6200.0 
HNCO+NH=NCO+NH2                1.040E+15      0.0     39390.0 
NCO+H=NH+CO                    5.360E+13      0.0         0.0 
NCO+O=NO+CO                    4.200E+13      0.0         0.0 
NCO+O=N+CO2                    8.000E+12      0.0      2500.0 
NCO+N=N2+CO                    2.000E+13      0.0         0.0 
NCO+OH=NO+HCO                  5.000E+12      0.0     15000.0 
NCO+M=N+CO+M                   2.200E+14      0.0     54050.0 
NCO+NO=N2O+CO                  4.600E+18     -2.01      934.0 
NCO+NO=N2+CO2                  5.800E+18     -2.01      934.0 
NCO+O2=NO+CO2                  2.000E+12      0.0     20000.0 
NCO+HCO=HNCO+CO                3.600E+13      0.0         0.0 
NCO+NO2=CO+NO+NO               2.830E+13     -0.646    -326.0 
NCO+NO2=CO2+N2O                3.570E+14     -0.646    -326.0 
NCO+HNO=HNCO+NO                1.800E+13      0.0         0.0 
NCO+NCO=CO+CO+N2               3.000E+12      0.0         0.0 
NO+HCO=CO+HNO                  7.240E+13     -0.4         0.0 
NO2+CO=CO2+NO                  9.000E+13      0.0     33800.0 
NO2+HCO=H+CO2+NO               8.400E+15     -0.75     1930.0 
CH3O+NO2=HONO+CH2O             3.000E+12      0.0         0.0 
CH3O+NO=CH2O+HNO               1.300E+14     -0.7         0.0 
NO2+CH2O=HONO+HCO              1.000E+10      0.0     15100.0 
NO+CH2O=HNO+HCO                1.000E+13      0.0     40820.0 
NO2+HCO=HONO+CO                1.000E+13      0.0         0.0 
NO2+HCO=OH+NO+CO               1.000E+14      0.0         0.0 
NCO+N=NO+CN                    2.700E+18     -0.995   17200.0 
CN+CH4=HCN+CH3                 9.000E+04      2.64     -300.0 
C+NO=CO+N                      2.800E+13      0.0         0.0 
NH+CO2=HNO+CO                  1.000E+13      0.0     14350.0 
NCO+CH4=HNCO+CH3               1.000E+13      0.0      8130.0 
C+N2O=CN+NO                    4.800E+12      0.0         0.0 
CH+NH2=HCN+H+H                 3.000E+13      0.0         0.0 
CH+NH=HCN+H                    5.000E+13      0.0         0.0 
CH2+NH=HCN+H+H                 3.000E+13      0.0         0.0 
CH3+N=HCN+H+H                  2.000E+11      0.0         0.0 
CH3+N=HCN+H2                   7.100E+12      0.0         0.0 
CH4+N=NH+CH3                   1.000E+13      0.0     24000.0 
C3H3+N=HCN+C2H2                1.000E+13      0.0         0.0 
CH+N2O=HCN+NO                  1.340E+13      0.0      -510.0
CH+N2O=CO+H+N2                 5.200E+12      0.0      -510.0
!
END



\\
\\
\\  This is the therm file
\\
\\
THERMO ALL
 300.000  1000.000  5000.000
H                 L 6/94H   10   00   00   0G   200.000  6000.00  1000.0       1
 0.25000000E+01 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 0.25473660E+05-0.44668285E+00 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.25473660E+05-0.44668285E+00 0.26219035E+05    4
H2  REF ELEMENT   RUS 78H   20   00   00   0G   200.000  6000.00  1000.0       1
 0.29328305E+01 0.82659802E-03-0.14640057E-06 0.15409851E-10-0.68879615E-15    2
-0.81305582E+03-0.10243164E+01 0.23443029E+01 0.79804248E-02-0.19477917E-04    3
 0.20156967E-07-0.73760289E-11-0.91792413E+03 0.68300218E+00 0.00000000E+00    4
O                 L 1/90O   10   00   00   0G   200.000  6000.00  1000.0       1
 2.54363697E+00-2.73162486E-05-4.19029520E-09 4.95481845E-12-4.79553694E-16    2
 2.92260120E+04 4.92229457E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00 2.99687009E+04    4
O2 REF ELEMENT    RUS 89O   20   00   00   0G   200.000  6000.00  1000.0       1
 3.66096083E+00 6.56365523E-04-1.41149485E-07 2.05797658E-11-1.29913248E-15    2
-1.21597725E+03 3.41536184E+00 3.78245636E+00-2.99673415E-03 9.84730200E-06    3
-9.68129508E-09 3.24372836E-12-1.06394356E+03 3.65767573E+00 0.00000000E+00    4
OH                RUS 78O   1H   10   00   0G   200.000  6000.00  1000.0       1
 2.83864607E+00 1.10725586E-03-2.93914978E-07 4.20524247E-11-2.42169092E-15    2
 3.94395852E+03 5.84452662E+00 3.99201543E+00-2.40131752E-03 4.61793841E-06    3
-3.88113333E-09 1.36411470E-12 3.61508056E+03-1.03925458E-01 4.73234213E+03    4
HO2               L 5/89H   1O   20   00   0G   200.000  6000.00  1000.0       1
 0.41722659E+01 0.18812098E-02-0.34629297E-06 0.19468516E-10 0.17609153E-15    2
 0.61818851E+02 0.29577974E+01 0.43017880E+01-0.47490201E-02 0.21157953E-04    3
-0.24275961E-07 0.92920670E-11 0.29480876E+03 0.37167010E+01 0.15096500E+04    4
H2O               L 5/89H   2O   10   00   0G   200.000  6000.00  1000.0       1
 0.26770389E+01 0.29731816E-02-0.77376889E-06 0.94433514E-10-0.42689991E-14    2
-0.29885894E+05 0.68825500E+01 0.41986352E+01-0.20364017E-02 0.65203416E-05    3
-0.54879269E-08 0.17719680E-11-0.30293726E+05-0.84900901E+00-0.29084817E+05    4
H2O2              L 2/93H   2O   20   00   0G   200.000  6000.00  1000.0       1
 4.57333537E+00 4.04984070E-03-1.29479479E-06 1.97281710E-10-1.13402846E-14    2
-1.80548121E+04 7.04278488E-01 4.27611269E+00-5.42822417E-04 1.67335701E-05    3
-2.15770813E-08 8.62454363E-12-1.77542989E+04 3.43505074E+00-1.63942313E+04    4
CO                RUS 79C   1O   10   00   0G   200.000  6000.00  1000.0       1
 0.30484859E+01 0.13517281E-02-0.48579405E-06 0.78853644E-10-0.46980746E-14    2
-0.14266117E+05 0.60170977E+01 0.35795335E+01-0.61035369E-03 0.10168143E-05    3
 0.90700586E-09-0.90442449E-12-0.14344086E+05 0.35084093E+01-0.13293628E+05    4
CO2               L 7/88C   1O   20   00   0G   200.000  6000.00  1000.0       1
 0.46365111E+01 0.27414569E-02-0.99589759E-06 0.16038666E-09-0.91619857E-14    2
-0.49024904E+05-0.19348955E+01 0.23568130E+01 0.89841299E-02-0.71220632E-05    3
 0.24573008E-08-0.14288548E-12-0.48371971E+05 0.99009035E+01-0.47328105E+05    4
HCO               L12/89H   1C   1O   10   0G   200.000  6000.00  1000.0       1
 3.64896209E+00 3.08090819E-03-1.12429876E-06 1.86308085E-10-1.13951828E-14    2
 3.71209048E+03 5.06147406E+00 4.22118584E+00-3.24392532E-03 1.37799446E-05    3
-1.33144093E-08 4.33768865E-12 3.83956496E+03 3.39437243E+00 5.05141013E+03    4
CH3               L 8/88C   1H   30   00   0G   200.000  6000.00  1000.0       1
 0.29669735E+01 0.57936672E-02-0.19694809E-05 0.30556936E-09-0.17767843E-13    2
 0.16537116E+05 0.47918803E+01 0.36733375E+01 0.20020559E-02 0.57853135E-05    3
-0.69873054E-08 0.26055599E-11 0.16440180E+05 0.16018315E+01 0.17662905E+05    4
CH4   ANHARMONIC  L 8/88C   1H   40   00   0G   200.000  6000.00  1000.0       1
 0.16354256E+01 0.10084431E-01-0.33692369E-05 0.53497280E-09-0.31552817E-13    2
-0.10005603E+05 0.99936953E+01 0.51498792E+01-0.13671008E-01 0.49180130E-04    3
-0.48474403E-07 0.16669441E-10-0.10246648E+05-0.46413244E+01-0.89722666E+04    4
C2H6      ETHANE  L 5/84C   2H   60   00   0G   300.000  5000.00  1000.0       1
 0.47028847E+01 0.14042635E-01-0.46469377E-05 0.67473738E-09-0.35089312E-13    2
-0.12671988E+05-0.45433950E+01 0.15395260E+01 0.15040841E-01 0.66847115E-05    3
-0.13382948E-07 0.48561398E-11-0.11248766E+05 0.14107375E+02-0.10085085E+05    4
CH2O              L 8/88H   2C   1O   10   0G   200.000  6000.00  1000.0       1
 0.31694807E+01 0.61932742E-02-0.22505981E-05 0.36598245E-09-0.22015410E-13    2
-0.14478425E+05 0.60423533E+01 0.47937036E+01-0.99081518E-02 0.37321459E-04    3
-0.37927902E-07 0.13177015E-10-0.14308955E+05 0.60288702E+00-0.13059098E+05    4
C2H5              T12/91C   2H   5    0    0G   200.000  6000.00  1000.0       1
 0.42878814E+01 0.12433893E-01-0.44139119E-05 0.70654102E-09-0.42035136E-13    2
 0.12056455E+05 0.84602583E+00 0.43058580E+01-0.41833638E-02 0.49707270E-04    3
-0.59905874E-07 0.23048478E-10 0.12841714E+05 0.47100236E+01 0.14271225E+05    4
CH2      TRIPLET  RUS 79C   1H   2    0    0G   200.000  6000.00  1000.0       1
 0.25387122E+01 0.38225491E-02-0.12861304E-05 0.19800308E-09-0.11465743E-13    2
 0.46129253E+05 0.81054648E+01 0.41793655E+01-0.22178553E-02 0.79653602E-05    3
-0.69127339E-08 0.22475318E-11 0.45750857E+05-0.76113703E-02 0.46956131E+05    4
CH3O              L 8/88C   1H   3O   1    0G   200.000  6000.00  1000.0       1
 0.42653308E+01 0.78576406E-02-0.28410438E-05 0.46045190E-09-0.27631906E-13    2
 0.59351066E+02 0.39309947E+00 0.32652337E+01 0.33031665E-02 0.17048801E-04    3
-0.22709630E-07 0.88071768E-11 0.73229580E+03 0.74257357E+01 0.19625450E+04    4
CH2OH             T 9/92C   1H   3O   1    0G   200.000  6000.00  1000.0       1
 0.46776733E+01 0.65608772E-02-0.22636864E-05 0.35530469E-09-0.20842959E-13    2
-0.28928048E+04 0.48041879E+00 0.38644027E+01 0.55913161E-02 0.59514635E-05    3
-0.10477176E-07 0.43793199E-11-0.25050493E+04 0.54710652E+01-0.10704179E+04    4
CH                RUS 79C   1H   1    0    0G   200.000  6000.00  1000.0       1
 0.25209062E+01 0.17653726E-02-0.46147581E-06 0.59288567E-10-0.33473209E-14    2
 0.71134082E+05 0.74053223E+01 0.34898166E+01 0.32383554E-03-0.16889906E-05    3
 0.31621733E-08-0.14060907E-11 0.70799939E+05 0.20840111E+01 0.71845485E+05    4
C2H2              L 8/88C   2H   2    0    0G   200.000  6000.00  1000.0       1
 0.46587047E+01 0.48840949E-02-0.16083563E-05 0.24698787E-09-0.13861505E-13    2
 0.25663218E+05-0.39979074E+01 0.80869108E+00 0.23361395E-01-0.35516636E-04    3
 0.28014566E-07-0.85004459E-11 0.26332764E+05 0.13939671E+02 0.27349778E+05    4
C2H4    ETHYLENE  L 4/85C   2H   4    0    0G   300.000  5000.00  1000.0       1
 0.43985453E+01 0.96228607E-02-0.31663776E-05 0.45747628E-09-0.23659406E-13    2
 0.41153203E+04-0.24627438E+01 0.12176600E+01 0.13002675E-01 0.35037447E-05    3
-0.11155514E-07 0.47203222E-11 0.53373828E+04 0.15480169E+02 0.62902830E+04    4
C2H3              T06/93C   2H   3    0    0G   200.000  6000.00  1000.0       1
 0.47025310E+01 0.72642283E-02-0.25801992E-05 0.41319944E-09-0.24591492E-13    2
 0.34029675E+05-0.14293714E+01 0.30019602E+01 0.30304354E-02 0.24444315E-04    3
-0.35810242E-07 0.15108700E-10 0.34868173E+05 0.93304495E+01 0.36050230E+05    4
CH3OH             L 8/88C   1H   4O   1    0G   200.000  6000.00  1000.0       1
 0.36012593E+01 0.10243223E-01-0.35999217E-05 0.57251951E-09-0.33912719E-13    2
-0.25997155E+05 0.47056025E+01 0.57153948E+01-0.15230920E-01 0.65244182E-04    3
-0.71080873E-07 0.26135383E-10-0.25642765E+05-0.15040970E+01-0.24167389E+05    4
CH3HCO   ACETALDE L 8/88C   2H   4O   1    0G   200.000  6000.00  1000.0       1
 0.54041108E+01 0.11723059E-01-0.42263137E-05 0.68372451E-09-0.40984863E-13    2
-0.22593122E+05-0.34807917E+01 0.47294595E+01-0.31932858E-02 0.47534921E-04    3
-0.57458611E-07 0.21931112E-10-0.21572878E+05 0.41030159E+01-0.19987949E+05    4
C2H4O OXYRANE     L 8/88C   2H   4O   1    0G   200.000  6000.00  1000.0       1
 0.54887641E+01 0.12046190E-01-0.43336931E-05 0.70028311E-09-0.41949088E-13    2
-0.91804251E+04-0.70799605E+01 0.37590532E+01-0.94412180E-02 0.80309721E-04    3
-0.10080788E-06 0.40039921E-10-0.75608143E+04 0.78497475E+01-0.63304657E+04    4
C2H               T04/93C   2H   1    0    0G   200.000  6000.00  1000.0       1
 0.36646060E+01 0.38218694E-02-0.13650743E-05 0.21324828E-09-0.12309430E-13    2
 0.64597055E+05 0.39134973E+01 0.29018020E+01 0.13285982E-01-0.28050886E-04    3
 0.28930184E-07-0.10744742E-10 0.64490326E+05 0.61723506E+01 0.65750290E+05    4
CH2CO KETENE      T 6/94C   2H   2O   1    0G   200.000  6000.00  1000.0       1
 0.57577901E+01 0.63496507E-02-0.22584407E-05 0.36208462E-09-0.21569030E-13    2
-0.79786113E+04-0.61064037E+01 0.21401165E+01 0.18088368E-01-0.17324216E-04    3
 0.92767477E-08-0.19915011E-11-0.70430509E+04 0.12198699E+02-0.57366700E+04    4
HCCO              T 6/94C   2H   1O   1    0G   200.000  6000.00  1000.0       1
 0.58469006E+01 0.36405960E-02-0.12959007E-05 0.20796919E-09-0.12400022E-13    2
 0.19248496E+05-0.52916533E+01 0.23350118E+01 0.17010083E-01-0.22018867E-04    3
 0.15406447E-07-0.43455097E-11 0.20050299E+05 0.11976729E+02 0.21336387E+05    4
SCH2   SINGLET    C12/87C   1H   2    0     G   300.00   4000.00  1000.0       1
 0.03552888E+02 0.02066788E-01-0.01914116E-05-0.11046733E-09 0.02021349E-12    2
 0.04984975E+06 0.01686570E+02 0.03971265E+02-0.01699088E-02 0.10253689E-05    3
 0.02492550E-07-0.01981266E-10 0.04989367E+06 0.05753207E+00 0.51077663E+05    4
C2                RUS 79C   2    0    0    0G   200.000  6000.00  1000.0       1
 0.37913706E+01 0.51650473E-03-0.25486960E-07-0.82263554E-11 0.10086168E-14    2
 0.99023059E+05 0.28151802E+01 0.86470550E+00 0.39353120E-01-0.11981818E-03    3
 0.13908103E-06-0.55205503E-10 0.98731303E+05 0.11530141E+02 0.99928438E+05    4
C2O               RUS 79C   2O   1    0    0G   200.000  6000.00  1000.0       1
 0.51512722E+01 0.23726722E-02-0.76135971E-06 0.11706415E-09-0.70257804E-14    2
 0.33241888E+05-0.22183135E+01 0.28648610E+01 0.11990216E-01-0.18362448E-04    3
 0.15769739E-07-0.53897452E-11 0.33749932E+05 0.88867772E+01 0.35003406E+05    4
CH3CO             T 9/92C   2H   3O   1    0G   200.000  6000.00  1000.0       1
 0.59447731E+01 0.78667205E-02-0.28865882E-05 0.47270875E-09-0.28599861E-13    2
-0.37873075E+04-0.50136751E+01 0.41634257E+01-0.23261610E-03 0.34267820E-04    3
-0.44105227E-07 0.17275612E-10-0.26574529E+04 0.73468280E+01-0.12027167E+04    4
CH3CO3     4/24/92 THERMC   2H   3O   3    0G   300.000  5000.000 1575.000    21
 1.14020019E+01 8.59416721E-03-3.06838374E-06 4.90245722E-10-2.90175354E-14    2
-2.02745353E+04-3.30871657E+01 5.13747438E+00 1.59298658E-02-8.21627846E-07    3
-4.81982007E-09 1.57730795E-12-1.73772288E+04 3.13626395E+00                   4
CH3CO3H    4/24/92 THERMC   2H   4O   3    0G   300.000  5000.000 1374.000    21
 1.23794767E+01 1.00368883E-02-3.50324632E-06 5.51431811E-10-3.23046627E-14    2
-3.88133888E+04-3.71309181E+01 4.09002627E+00 2.54380593E-02-1.22828804E-05    3
 1.55365086E-09 2.77282260E-13-3.55663564E+04 8.74816008E+00                   4
CH3O2             L 1/84C   1H   3O   2    0G   300.000  5000.00   1000.0      1
 0.66812963E+01 0.80057271E-02-0.27188507E-05 0.40631365E-09-0.21927725E-13    2
 0.52621851E+03-0.99423847E+01 0.20986490E+01 0.15786357E-01 0.75683261E-07    3
-0.11274587E-07 0.56665133E-11 0.20695879E+04 0.15007068E+02 0.33715510E+04    4
CH3O2H            T11/96C   1H   4O   2    0G   200.000  6000.00  1000.0       1
 6.86907934E+00 1.00840883E-02-3.66515947E-06 5.96302681E-10-3.58894156E-14    2
-1.98402231E+04-1.24951986E+01 3.72654981E+00 7.51851847E-03 2.35970425E-05    3
-3.52694507E-08 1.42757614E-11-1.83982011E+04 6.52443362E+00-1.68074366E+04    4
C2H5O2H           T10/96C   2H   6O   2    0G   200.000  6000.00  1000.0       1
 9.99511555E+00 1.47311626E-02-5.30621235E-06 8.58442516E-10-5.14814807E-14    2
-2.53850722E+04-2.53504050E+01 4.37310002E+00 1.04422436E-02 4.63854723E-05    3
-7.02772770E-08 2.93034879E-11-2.29362227E+04 8.30134323E+00-2.08834916E+04    4
C2H5O2                  C   2H   5O   20   0G   300.00   5000.00  1000.00      1
 0.78275410E+01 0.15347060E-01-0.58108700E-05 0.10465070E-08-0.71520540E-13    2
-0.45527100E+04-0.12726800E+02 0.39390180E+01 0.13991490E-01 0.21720800E-04    3
-0.36137480E-07 0.14824910E-10-0.28469040E+04 0.10394470E+02                   4
CH3CO2     4/24/92 THERMC   2H   3O   2    0G   300.000  5000.000 1482.000    11
 8.37781397E+00 7.83020766E-03-2.53555405E-06 3.78870534E-10-2.13859882E-14    2
-2.38389350E+04-1.58386711E+01 4.30690287E+00 1.70945538E-02-1.02260887E-05    3
 3.08829567E-09-3.50372173E-13-2.24325681E+04 6.05640546E+00                   4
CH3CO2H    4/24/92 THERMC   2H   4O   2    0G   300.000  5000.000 1997.000    21
 7.34453556E+00 1.17718973E-02-4.13228619E-06 6.55097801E-10-3.86190560E-14    2
-4.51365468E+04-1.26872633E+01 3.34602242E+00 1.84778690E-02-6.73574398E-06    3
 1.07401993E-10 2.97848907E-13-4.35715193E+04 9.61466866E+00                   4
C2H5O             T11/82O   1C   2H   5    0G   300.000  5000.0   1000.0       1
 0.60114346E+01 0.12165219E-01-0.40449604E-05 0.59076588E-09-0.30969595E-13    2
-0.49366992E+04-0.67901798E+01 0.17302504E+01 0.16908489E-01 0.39996221E-05    3
-0.13711180E-07 0.57643603E-11-0.32922483E+04 0.17336115E+02-0.20138288E+04    4
C2H5OH            L 8/88C   2H   6O   1    0G   200.000  6000.00  1000.0       1
 0.65624365E+01 0.15204222E-01-0.53896795E-05 0.86225011E-09-0.51289787E-13    2
-0.31525621E+05-0.94730202E+01 0.48586957E+01-0.37401726E-02 0.69555378E-04    3
-0.88654796E-07 0.35168835E-10-0.29996132E+05 0.48018545E+01-0.28257829E+05    4
SC2H5O            T 4/83C   2O   1H   5    0G   300.000  5000.00  1000.0       1
 0.67665424E+01 0.11634436E-01-0.37790651E-05 0.53828875E-09-0.27315345E-13    2
-0.56092969E+04-0.93980442E+01 0.24813328E+01 0.16790036E-01 0.37755499E-05    3
-0.13923497E-07 0.60095193E-11-0.40120054E+04 0.14581622E+02-0.25172860E+04    4
PC2H5O         T 4/83   C   2H   5O   10   0G   300.00   5000.0   1000.00      1
 0.75944014E+01 0.93229339E-02-0.30303854E-05 0.43216319E-09-0.21970039E-13    2
-0.57727852E+04-0.13968735E+02 0.14019508E+01 0.21543175E-01-0.22326512E-05    3
-0.14464092E-07 0.80488420E-11-0.38464519E+04 0.19135818E+02                   4
CH2HCO           T04/83 C   2H   3O   10   0G   300.00   5000.00  1000.00      1
 0.59756699E+01 0.81305914E-02-0.27436245E-05 0.40703041E-09-0.21760171E-13    2
 0.49032178E+03-0.50452509E+01 0.34090624E+01 0.10738574E-01 0.18914925E-05    3
-0.71585831E-08 0.28673851E-11 0.15214766E+04 0.95582905E+01                   4
CN                T 6/94C   1N   1    0    0G   200.000  6000.00  1000.0       1
 0.37459804E+01 0.43450773E-04 0.29705984E-06-0.68651804E-10 0.44134174E-14    2
 0.52353188E+05 0.27867600E+01 0.36129350E+01-0.95551327E-03 0.21442976E-05    3
-0.31516324E-09-0.46430356E-12 0.52525340E+05 0.39804995E+01 0.52571034E+05    4
H2CN RADICAL      T05/97H   2C   1N   1    0G   200.000  6000.00  1000.0       1
 3.80315523E+00 5.47197456E-03-1.95314927E-06 3.13362513E-10-1.86249463E-14    2
 2.73218196E+04 3.31721893E+00 3.97799541E+00-3.43275678E-03 2.59134226E-05    3
-3.04692133E-08 1.16272702E-11 2.76769528E+04 4.43029598E+00 2.88846366E+04    4
N                 L 6/88N   1    0    0    0G   200.000  6000.00  1000.0       1
 0.24159429E+01 0.17489065E-03-0.11902369E-06 0.30226244E-10-0.20360983E-14    2
 0.56133775E+05 0.46496095E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.56104638E+05 0.41939088E+01 0.56850013E+05    4
NH                L11/89N   1H   1    0    0G   200.000  6000.00  1000.0       1
 0.27836929E+01 0.13298429E-02-0.42478047E-06 0.78348504E-10-0.55044470E-14    2
 0.42134514E+05 0.57407798E+01 0.34929084E+01 0.31179197E-03-0.14890484E-05    3
 0.24816442E-08-0.10356967E-11 0.41894294E+05 0.18483277E+01 0.42940822E+05    4
HCN               L 7/88H   1C   1N   1    0G   200.000  6000.00  1000.0       1
 0.38022392E+01 0.31464227E-02-0.10632185E-05 0.16619757E-09-0.97997567E-14    2
 0.14910512E+05 0.15754601E+01 0.22589885E+01 0.10051170E-01-0.13351763E-04    3
 0.10092349E-07-0.30089029E-11 0.15215853E+05 0.89164418E+01 0.16236675E+05    4
NO                RUS 89N   1O   1    0    0G   200.000  6000.00  1000.0       1
 3.26071234E+00 1.19101135E-03-4.29122646E-07 6.94481463E-11-4.03295681E-15    2
 9.92143132E+03 6.36900518E+00 4.21859896E+00-4.63988124E-03 1.10443049E-05    3
-9.34055507E-09 2.80554874E-12 9.84509964E+03 2.28061001E+00 1.09770882E+04    4
HCNO              120186H   1C   1N   1O   1G  0250.00   4000.00  1000.00      1
 0.06692412E+02 0.02368360E-01-0.02371510E-05-0.01275503E-08 0.02407137E-12    2
 0.01694737E+06-0.01245434E+03 0.03184859E+02 0.09752316E-01-0.01280203E-04    3
-0.06163104E-07 0.03226275E-10 0.01797907E+06 0.06123844E+02                   4
C                 L 7/88C   1    0    0    0G   200.000  6000.00  1000.0       1
 0.26055830E+01-0.19593434E-03 0.10673722E-06-0.16423940E-10 0.81870580E-15    2
 0.85411742E+05 0.41923868E+01 0.25542395E+01-0.32153772E-03 0.73379223E-06    3
-0.73223487E-09 0.26652144E-12 0.85442681E+05 0.45313085E+01 0.86195097E+05    4
HOCN              110193H   1C   1N   1O   1G   300.00   4000.00  1400.00      1
 0.06022112E+02 0.01929530E-01-0.01455029E-05-0.01045811E-08 0.01794814E-12    2
-0.04040321E+05-0.05866433E+02 0.03789424E+02 0.05387981E-01-0.06518270E-05    3
-0.01420164E-07 0.05367969E-11-0.03135335E+05 0.06667052E+02                   4
HNCO              T 6/94H   1N   1C   1O   1G   200.000  6000.00  1000.0       1
 0.52936894E+01 0.40307770E-02-0.14130589E-05 0.22445562E-09-0.13287683E-13    2
-0.15973489E+05-0.30864710E+01 0.22432188E+01 0.14491349E-01-0.15236174E-04    3
 0.83345851E-08-0.17104033E-11-0.15233708E+05 0.12157321E+02-0.14039745E+05    4
NCO               T 6/94C   1N   1O   1    0G   200.000  6000.00  1000.0       1
 0.51075979E+01 0.23371500E-02-0.88984637E-06 0.14920037E-09-0.91663122E-14    2
 0.14024945E+05-0.22908127E+01 0.27405490E+01 0.95089992E-02-0.10338762E-04    3
 0.68805052E-08-0.20963552E-11 0.14690320E+05 0.98908197E+01 0.15851325E+05    4
N2O               L 7/88N   2O   1    0    0G   200.000  6000.00  1000.0       1
 0.48230729E+01 0.26270251E-02-0.95850872E-06 0.16000712E-09-0.97752302E-14    2
 0.80734047E+04-0.22017208E+01 0.22571502E+01 0.11304728E-01-0.13671319E-04    3
 0.96819803E-08-0.29307182E-11 0.87417746E+04 0.10757992E+02 0.98141682E+04    4
NH2               L12/89N   1H   2    0    0G   200.000  6000.00  1000.0       1
 0.28476611E+01 0.31428453E-02-0.89866557E-06 0.13032357E-09-0.74885356E-14    2
 0.21823916E+05 0.64718133E+01 0.42055601E+01-0.21355282E-02 0.72682021E-05    3
-0.59302799E-08 0.18067218E-11 0.21535223E+05-0.14663231E+00 0.22747541E+05    4
HNO    WRA032498        H   1N   1O   1    0G   200.000  6000.00  1000.0       1
 3.16554762E+00 3.00005132E-03-3.94350282E-07-3.85787491E-11 7.08091931E-15    2
 1.18052184E+04 7.64764695E+00 4.53525882E+00-5.68546910E-03 1.85199976E-05    3
-1.71883674E-08 5.55833090E-12 1.16506820E+04 1.74314734E+00 1.28824477E+04    4
NO2               L 7/88N   1O   2    0    0G   200.000  6000.00  1000.0       1
 0.48847540E+01 0.21723955E-02-0.82806909E-06 0.15747510E-09-0.10510895E-13    2
 0.23164982E+04-0.11741695E+00 0.39440312E+01-0.15854290E-02 0.16657812E-04    3
-0.20475426E-07 0.78350564E-11 0.28966180E+04 0.63119919E+01 0.41124701E+04    4
C2N2              RUS 79C   2N   2    0    0G   200.000  6000.00  1000.0       1
 0.67055078E+01 0.36425829E-02-0.13094063E-05 0.21643797E-09-0.13121437E-13    2
 0.34860766E+05-0.10493904E+02 0.23292532E+01 0.26153785E-01-0.49000399E-04    3
 0.46191748E-07-0.16432385E-10 0.35668442E+05 0.98501993E+01 0.37175973E+05    4
NNH               T07/93N   2H   1    0    0G   200.000  6000.00  1000.0       1
 0.37667545E+01 0.28915081E-02-0.10416620E-05 0.16842594E-09-0.10091896E-13    2
 0.28650697E+05 0.44705068E+01 0.43446927E+01-0.48497072E-02 0.20059459E-04    3
-0.21726464E-07 0.79469538E-11 0.28791973E+05 0.29779411E+01 0.30009829E+05    4
NH3  AMONIA       RUS 89N   1H   3    0    0G   200.000  6000.00  1000.0       1
 2.71709692E+00 5.56856338E-03-1.76886396E-06 2.67417260E-10-1.52731419E-14    2
-6.58451989E+03 6.09289837E+00 4.30177808E+00-4.77127330E-03 2.19341619E-05    3
-2.29856489E-08 8.28992268E-12-6.74806394E+03-6.90644393E-01-5.52528050E+03    4
N2H2              L 5/90N   2H   2    0    0G   200.000  6000.00  1000.0       1
 0.13111509E+01 0.90018727E-02-0.31491187E-05 0.48144969E-09-0.27189798E-13    2
 0.24786417E+05 0.16409109E+02 0.49106602E+01-0.10779187E-01 0.38651644E-04    3
-0.38650163E-07 0.13485210E-10 0.24224273E+05 0.91027970E-01 0.25480756E+05    4
CH2CHOW      cc      cc C   2H   3O   1    0G   300.000  5000.000 1389.000    01
 6.44647136E+00 7.67418446E-03-2.55298889E-06 3.88933291E-10-2.22665944E-14    2
 1.95121257E+04-7.93000911E+00 2.76406282E+00 1.50800828E-02-7.88093083E-06    3
 1.97116311E-09-1.74094646E-13 2.09267565E+04 1.22979303E+01                   4
C2H3O             T 9/92C   2H   3O   1    0G   298.150  3000.000  1000.0      1
 0.48131470E+00 0.20711914E-01-0.12693155E-04 0.34579642E-08-0.35399703E-12    2
 0.15648642E+05 0.34629876E+02 0.10854772E+01 0.12845259E-01 0.24138660E-05    3
-0.44642672E-08-0.29381916E-12 0.15910655E+05 0.33395312E+02 0.16817588E+05    4
CH3HCOW      cc      cc C   2H   4O   1    0G   300.000  5000.000 1376.000    11
 6.89091884E+00 9.95180916E-03-3.45461153E-06 5.41621230E-10-3.16369494E-14    2
 1.91997658E+04-1.23847501E+01 1.66831767E+00 1.89949607E-02-8.35817476E-06    3
 1.20243883E-09 6.37236161E-14 2.13971005E+04 1.68980856E+01                   4
C3H6OH     4/24/92 THERMC   3H   7O   1    0G   300.000  5000.000 1403.000    31
 1.09340809E+01 1.49384611E-02-4.97490474E-06 7.58435125E-10-4.34420191E-14    2
-1.15876725E+04-2.89343754E+01 1.65867207E-01 4.13141037E-02-2.99664397E-05    3
 1.16234591E-08-1.86270510E-12-7.97925983E+03 2.84481654E+01                   4
O2C3H6OH   4/24/92 THERMC   3H   7O   3    0G   300.000  5000.000 1380.000    41
 1.55655181E+01 1.68822736E-02-5.91144639E-06 9.32274430E-10-5.46815496E-14    2
-2.95301139E+04-5.17312137E+01 2.20217327E+00 4.59913087E-02-3.03766163E-05    3
 1.05366009E-08-1.55928684E-12-2.45762626E+04 2.09594636E+01                   4
C3H5O2     4/24/92 THERMC   3H   5O   2    0G   300.000  5000.000 1370.000    21
 1.27832465E+01 1.25399349E-02-4.43433905E-06 7.03952842E-10-4.14804136E-14    2
 5.58979899E+03-4.03574985E+01 2.65293041E+00 3.19993514E-02-1.73652420E-05    3
 4.00148647E-09-2.51071126E-13 9.58815244E+03 1.56431929E+01                   4
C3H5O2H    4/24/92 THERMC   3H   6O   2    0G   300.000  5000.000 1386.000    31
 1.39009363E+01 1.34732295E-02-4.70091523E-06 7.39632892E-10-4.33128479E-14    2
-1.29274765E+04-4.49719922E+01 1.75161845E+00 4.16542055E-02-3.01268606E-05    3
 1.14209488E-08-1.80624472E-12-8.64152220E+03 2.04277976E+01                   4
C3H5O      4/24/92 THERMC   3H   5O   1    0G   300.000  5000.000 1375.000    21
 9.93380796E+00 1.21383493E-02-4.23299952E-06 6.65805507E-10-3.89817563E-14    2
 6.20904655E+03-2.68841926E+01 1.48979417E+00 2.77306345E-02-1.35044719E-05    3
 2.24225811E-09 8.28232368E-14 9.58557792E+03 1.99841620E+01                   4
NC3H7O2    4/24/92 THERMC   3H   7O   2    0G   300.000  5000.000 1376.000    31
 1.39295562E+01 1.62109297E-02-5.69426414E-06 8.99917649E-10-5.28613664E-14    2
-1.09646348E+04-4.81871702E+01 1.65182379E+00 4.10402887E-02-2.41422856E-05    3
 6.91955127E-09-7.89713852E-13-6.22330250E+03 1.92752981E+01                   4
NC3H7O2H   4/24/92 THERMC   3H   8O   2    0G   300.000  5000.000 1392.000    31
 1.45695935E+01 1.77104023E-02-6.09654929E-06 9.50484279E-10-5.53030086E-14    2
-2.92401303E+04-4.99984239E+01 6.27002549E-01 5.14448631E-02-3.83791359E-05    3
 1.54653552E-08-2.61561123E-12-2.44385400E+04 2.45962330E+01                   4
IC3H7O2    4/24/92 THERMC   3H   7O   2    0G   300.000  5000.000 1382.000    31
 1.45333419E+01 1.56502508E-02-5.48993478E-06 8.66859600E-10-5.08886097E-14    2
-1.33000915E+04-5.13639953E+01 1.03606539E+00 4.57999209E-02-3.15792877E-05    3
 1.14069765E-08-1.73998862E-12-8.39327437E+03 2.17529616E+01                   4
IC3H7O2H   4/24/92 THERMC   3H   8O   2    0G   300.000  5000.000 1395.000    31
 1.51552007E+01 1.71807437E-02-5.90648731E-06 9.20033944E-10-5.34977754E-14    2
-3.15766321E+04-5.30871350E+01 1.35497567E-01 5.55626488E-02-4.48202480E-05    3
 1.93393698E-08-3.43326450E-12-2.66246790E+04 2.65182104E+01                   4
IC3H7O       cc      cc C   3H   7O   1    0G   300.000  5000.000 1379.000    21
 1.24938036E+01 1.42024822E-02-4.76949687E-06 7.32349424E-10-4.21933677E-14    2
-1.25023025E+04-4.32454233E+01-1.45132586E-01 4.09328495E-02-2.46708437E-05    3
 6.60628776E-09-5.34569348E-13-7.91212088E+03 2.54949227E+01                   4
NC3H7O N-PROPOXY  T 3/96C   3H   7O   1    0G   298.150  5000.00  1000.0       1
 0.84124958E+01 0.19520193E-01-0.71317071E-05 0.12393621E-08-0.82483889E-13    2
-0.87750718E+04-0.18293360E+02 0.91452571E+00 0.33601264E-01-0.12282254E-04    3
-0.10739947E-08 0.72924952E-12-0.61847956E+04 0.22563171E+02-0.45289500E+04    4
IC4H5           I-A 8/83C   4H   5    0    0G   300.00   3000.0   1000.0       1
 0.72559180E+01 0.17141310E-01-0.74350450E-05 0.16410520E-08-0.15360780E-12    2
 0.36339860E+05-0.13420580E+02 0.10251400E+02-0.25627670E-01 0.10715090E-03    3
-0.11090190E-06 0.37550050E-10 0.37128550E+05-0.20574530E+02 0.39785870E+05    4
NC4H5           N-A 8/83C   4H   5    0    0G   300.0    3000.0   1000.0       1
 0.72559180E+01 0.17141310E-01-0.74350450E-05 0.16410520E-08-0.15360780E-12    2
 0.38227410E+05-0.12422830E+02-0.20329280E+01 0.53630530E-01-0.63382690E-04    3
 0.40641440E-07-0.10398020E-10 0.40219650E+05 0.32787830E+02 0.41505750E+05    4
C6H5  PHENYL RAD  L12/84C   6H   5    0    0G   300.000  5000.00  1000.0       1
 0.11431418E+02 0.17019045E-01-0.58387241E-05 0.88094687E-09-0.48050417E-13    2
 0.33942348E+05-0.38574219E+02-0.23405075E+01 0.42760305E-01-0.25518166E-05    3
-0.30668716E-07 0.16245519E-10 0.38376734E+05 0.35617355E+02 0.39502977E+05    4
C4H4 1-BUTEN-3YN  L 9/89C   4H   4    0    0G   200.000  6000.00  1000.0       1
 0.82948104E+01 0.11994381E-01-0.42624075E-05 0.68306978E-09-0.40680631E-13    2
 0.33550866E+05-0.18426826E+02 0.14049083E+01 0.29531073E-01-0.15596302E-04    3
-0.32142002E-08 0.45436937E-11 0.35507830E+05 0.17450183E+02 0.37097268E+05    4
C3H3 PROPARGYL    T 5/97C   3H   3    0    0G   200.000  6000.00  1000.0       1
 7.14221880E+00 7.61902005E-03-2.67459950E-06 4.24914801E-10-2.51475415E-14    2
 3.89087427E+04-6.82838123E+00 1.35110927E+00 3.27411223E-02-4.73827135E-05    3
 3.76309808E-08-1.18540923E-11 4.01057783E+04 2.09623551E+01 4.16139977E+04    4
C3H6   PROPYLENE  L 4/85C   3H   6    0    0G   300.000  5000.00  1000.0       1
 0.67213974E+01 0.14931757E-01-0.49652353E-05 0.72510753E-09-0.38001476E-13    2
-0.92453149E+03-0.12155617E+02 0.14575157E+01 0.21142263E-01 0.40468012E-05    3
-0.16319003E-07 0.70475153E-11 0.10740208E+04 0.17399460E+02 0.24557265E+04    4
C3H8     PROPANE  L 4/85C   3H   8    0    0G   300.000  5000.00  1000.0       1
 0.75341368E+01 0.18872239E-01-0.62718491E-05 0.91475649E-09-0.47838069E-13    2
-0.16467516E+05-0.17892349E+02 0.93355381E+00 0.26424579E-01 0.61059727E-05    3
-0.21977499E-07 0.95149253E-11-0.13958520E+05 0.19201691E+02-0.12489986E+05    4
IC3H7           I-L 9/84C   3H   7    0    0G   300.000  5000.00  1000.0       1
 0.65294638E+01 0.17193288E-01-0.57153220E-05 0.83408080E-09-0.43663532E-13    2
 0.77179102E+04-0.91399021E+01 0.14461584E+01 0.20988975E-01 0.77172672E-05    3
-0.18481391E-07 0.71269024E-11 0.98206094E+04 0.20108200E+02 0.11221480E+05    4
NC3H7           N-L 9/84C   3H   7    0    0G   300.000  5000.00  1000.0       1
 0.77026987E+01 0.16044203E-01-0.52833220E-05 0.76298590E-09-0.39392284E-13    2
 0.82984336E+04-0.15480180E+02 0.10515518E+01 0.25991980E-01 0.23800540E-05    3
-0.19609569E-07 0.93732470E-11 0.10631863E+05 0.21122559E+02 0.12087447E+05    4
SC3H5  CH3CH=CH*  T 6/96C   3H   5    0    0G   298.150  5000.00  1000.0       1
 0.61342175E+01 0.13325025E-01-0.48645230E-05 0.83773904E-09-0.55161424E-13    2
 0.28626198E+05-0.74112134E+01 0.22361012E+01 0.16387245E-01 0.76639733E-05    3
-0.19364818E-07 0.84155468E-11 0.30173536E+05 0.14806060E+02 0.31602007E+05    4
PC3H4     PROPYNE T 2/90H   4C   3    0    0G   200.000  6000.00  1000.0       1
 0.60252400E+01 0.11336542E-01-0.40223391E-05 0.64376063E-09-0.38299635E-13    2
 0.19620942E+05-0.86043785E+01 0.26803869E+01 0.15799651E-01 0.25070596E-05    3
-0.13657623E-07 0.66154285E-11 0.20802374E+05 0.98769351E+01 0.22302059E+05    4
TC3H5  CH3C*=CH2  T 6/96C   3H   5    0    0G   200.000  6000.00  1000.0       1
 0.61101805E+01 0.14673395E-01-0.53676822E-05 0.86904932E-09-0.51932006E-13    2
 0.25532442E+05-0.83555712E+01 0.25544033E+01 0.10986798E-01 0.30174305E-04    3
-0.47253568E-07 0.19771073E-10 0.27150242E+05 0.13207592E+02 0.28582707E+05    4
C3H6O PROPYLENE OX L9/85C   3H   6O   1    0G   300.000  5000.00  1000.0       1
 0.86900558E+01 0.16020987E-01-0.53971753E-05 0.79941542E-09-0.42656366E-13    2
-0.15420691E+05-0.22485016E+02 0.48733836E+00 0.28519690E-01 0.30096162E-05    3
-0.22652642E-07 0.10706728E-10-0.12556434E+05 0.22605270E+02-0.11156446E+05    4
C3H4CY            T12/81C   3H   4    0    0G   300.000  5000.00  1000.0       1
 0.66999931E+01 0.10357372E-01-0.34551167E-05 0.50652949E-09-0.26682276E-13    2
 0.30199051E+05-0.13378770E+02-0.24621047E-01 0.23197215E-01-0.18474357E-05    3
-0.15927593E-07 0.86846155E-11 0.32334137E+05 0.22729762E+02 0.3332728 E+05    4
C3H2 *HC=C=CH*    T 7/94C   3H   2    0    0G   298.150  5000.00  1000.0       1
 0.62476406E+01 0.56755154E-02-0.19541329E-05 0.31686787E-09-0.19722328E-13    2
 0.79187553E+05-0.63712327E+01 0.36613541E+01 0.89364130E-02 0.35404269E-05    3
-0.11072967E-07 0.52009416E-11 0.80115193E+05 0.79774693E+01 0.81616355E+05    4
C4H2              L 9/89C   4H   2    0    0G   200.000  6000.00  1000.0       1
 0.86670035E+01 0.67166371E-02-0.23544995E-05 0.37383079E-09-0.22118914E-13    2
 0.49856933E+05-0.21114205E+02-0.39518508E+00 0.51955813E-01-0.91786616E-04    3
 0.80523929E-07-0.26917088E-10 0.51451709E+05 0.20969101E+02 0.52978651E+05    4
IC4H7       CH2=CHCHCH3 C   4H   70   00   0G   300.00   3000.00  1000.00      1
  0.5521420E+01  0.2683684E-01 -0.1286430E-04  0.3088640E-08 -0.3030856E-12    2
  0.1198044E+05 -0.4495640E+01 -0.1080514E+01  0.4638686E-01 -0.3464697E-04    3
  0.1401374E-07 -0.2395040E-11  0.1375538E+05  0.2933150E+02                   4
C4H8    1-BUTENE  T 6/83C   4H   8    0    0G   300.000  5000.00  1000.0       1
 0.20535841E+01 0.34350507E-01-0.15883197E-04 0.33089662E-08-0.25361045E-12    2
-0.21397231E+04 0.15556360E+02 0.11811380E+01 0.30853380E-01 0.50865247E-05    3
-0.24654888E-07 0.11110193E-10-0.17904004E+04 0.21075639E+02-0.06494670E+03    4
T2C4H8                  C   4H   80   00   0G   300.00   5000.00  1000.00      1
 0.82797676E+00 0.35864539E-01-0.16634498E-04 0.34732759E-08-0.26657398E-12    2
-0.30521033E+04 0.21342545E+02 0.12594252E+01 0.27808424E-01 0.87013932E-05    3
-0.24402205E-07 0.98977710E-11-0.29647742E+04 0.20501129E+02                   4
C2C4H8  C4H8CIS   T 6/83C   4H   8    0    0G   300.000  5000.00  1000.0       1
 0.11097383E+01 0.35542578E-01-0.16481703E-04 0.34412202E-08-0.26411468E-12    2
-0.26507607E+04 0.19366680E+02 0.24108791E+01 0.25147773E-01 0.98473047E-05    3
-0.22716758E-07 0.86585895E-11-0.27758694E+04 0.14097700E+02-0.89121300E+03    4
C4H5                    C   4H   50   00   0G   300.00   5000.00  1000.00      1
 0.52872400E+01 0.19416250E-01-0.75259100E-05 0.13124900E-08-0.84820000E-13    2
 0.38354600E+05-0.15660000E+01 0.27518800E+01 0.29092550E-01-0.23393090E-04    3
 0.14159370E-07-0.42054700E-11 0.38953300E+05 0.10953000E+02                   4
C4H               T12/91C   4H   1    0    0G   200.000  6000.00  1000.0       1
 0.77680939E+01 0.49850386E-02-0.17648839E-05 0.28217408E-09-0.16779623E-13    2
 0.93912126E+05-0.14159577E+02 0.13210657E+01 0.38562824E-01-0.71343174E-04    3
 0.65319977E-07-0.22607050E-10 0.95021629E+05 0.15554575E+02 0.96617600E+05    4
C6H6    BENZENE   L12/84C   6H   6    0    0G   300.000  5000.00  1000.0       1
 0.11815166E+02 0.19169778E-01-0.65425238E-05 0.98228425E-09-0.53280361E-13    2
 0.40707441E+04-0.43973511E+02-0.32181215E+01 0.47168836E-01-0.21254918E-05    3
-0.34879005E-07 0.18425386E-10 0.89017773E+04 0.36999313E+02 0.99586128E+04    4
C3H5 SYMMETRIC    T 9/96C   3H   5    0    0G   200.000  6000.00  1000.0       1
 0.70094568E+01 0.13106629E-01-0.46533442E-05 0.74514323E-09-0.44350051E-13    2
 0.16412909E+05-0.13946114E+02 0.14698036E+01 0.19034365E-01 0.14480425E-04    3
-0.35468652E-07 0.16647594E-10 0.18325831E+05 0.16724114E+02 0.19675772E+05    4
C3H4      ALLENE  L 8/89C   3H   4    0    0G   200.000  6000.00  1000.0       1
 0.63168722E+01 0.11133728E-01-0.39629378E-05 0.63564238E-09-0.37875540E-13    2
 0.20117495E+05-0.10995766E+02 0.26130445E+01 0.12122575E-01 0.18539880E-04    3
-0.34525149E-07 0.15335079E-10 0.21541567E+05 0.10226139E+02 0.22962267E+05    4
C4H6    1,3-Bdiene      C   4H   60   00   0G   300.00   3000.00  1000.00      1
 0.7130943E+01  0.2064008E-01 -0.9323906E-05  0.2147307E-08 -0.2061361E-12     2
 0.9612654E+04 -0.1417222E+02 -0.2530900E+01  0.5387969E-01 -0.5173933E-04     3
 0.2510720E-07 -0.4324726E-11  0.1187691E+05  0.3391388E+02  0.1333562E+05     4
IC4H3                   C   4H   30   00   0G   300.00   5000.00  1000.00      1
 0.97570450E+01 0.88991060E-02-0.36321540E-05 0.67966640E-09-0.47498080E-13    2
 0.52220940E+05-0.22567330E+02 0.92063570E+01-0.68412580E-02 0.49426770E-04    3
-0.57756540E-07 0.21620850E-10 0.53230880E+05-0.15490760E+02                   4
NC4H3                   C   4H   30   00   0G   300.00   5000.00  1000.00      1
 0.91143560E+01 0.95471120E-02-0.38843590E-05 0.72623270E-09-0.50737070E-13    2
 0.61621520E+05-0.21296650E+02 0.56986800E+01 0.60107180E-02 0.28706380E-04    3
-0.41907200E-07 0.16944030E-10 0.63201220E+05-0.49851380E+00                   4
H2C4O             120189H   2C   4O   1     G  0300.00   4000.00  1000.00      1
 0.10268878E+02 0.04896164E-01-0.04885080E-05-0.02708566E-08 0.05107013E-12    2
 0.02346902E+06-0.02815985E+03 0.04810971E+02 0.13139988E-01 0.09865073E-05    3
-0.06120720E-07 0.16400028E-11 0.02545803E+06 0.02113424E+02                   4
C6H5O       PHENOXY RAD C   6H   5O   10   0G   300.00   5000.00  1000.00      1
 0.13833984E+02 0.17618403E-01-0.60696257E-05 0.91988173E-09-0.50449181E-13    2
-0.69212549E+03-0.50392990E+02-0.18219433E+01 0.48122510E-01-0.46792302E-05    3
-0.34018594E-07 0.18649637E-10 0.42429180E+04 0.33526199E+02                   4
CH3NO3     T05/98       C   1H   3O   3N   1G   200.000  6000.000  1000.0      1
 9.77845489E+00 1.10069541E-02-4.25928645E-06 7.18198185E-10-4.42041793E-14    2
-1.88804487E+04-2.39163197E+01 3.91363583E+00 1.52137945E-02 1.73479131E-05    3
-3.37074473E-08 1.44322204E-11-1.66103232E+04 9.44208392E+00-1.46737980E+04    4
CH2NO3     T05/98       C   1H   2O   3N   1G   200.000  6000.000  1000.0      1
 1.03913885E+01 7.66103917E-03-3.02728077E-06 5.16124915E-10-3.19767406E-14    2
 7.78486241E+03-2.54151556E+01 2.98654023E+00 2.47990510E-02-1.17175684E-05    3
-5.36820166E-09 4.80947389E-12 1.00202588E+04 1.36939353E+01 1.19010741E+04    4
HONO         HNO2 RUS 89H   1N   1O   2    0G   200.000  6000.000  1000.0      1
 0.57919018E+01 0.36515212E-02-0.12928936E-05 0.20688716E-09-0.12315254E-13    2
-0.11565589E+05-0.40558233E+01 0.32141709E+01 0.81276869E-02 0.16602559E-05    3
-0.95285182E-08 0.48715058E-11-0.10753237E+05 0.98219504E+01-0.94355439E+04    4
AR REF ELEMENT    L 6/88AR  1    0    0    0G   200.000  6000.00  1000.0       1
 0.25000000E+01 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.74537500E+03 0.43796749E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.74537500E+03 0.43796749E+01 0.00000000E+00    4
N2  REF ELEMENT   RUS 78N   2    0    0    0G   200.000  6000.00  1000.0       1
 0.29525407E+01 0.13968838E-02-0.49262577E-06 0.78600091E-10-0.46074978E-14    2
-0.92393753E+03 0.58718221E+01 0.35309628E+01-0.12365950E-03-0.50299339E-06    3
 0.24352768E-08-0.14087954E-11-0.10469637E+04 0.29674391E+01 0.00000000E+00    4
NO3               J12/64N   1O   3    0    0G   200.000  6000.00  1000.0       1
 7.48347734E+00 2.57772041E-03-1.00945831E-06 1.72314072E-10-1.07154015E-14    2
 5.70919428E+03-1.41618155E+01 2.17359310E+00 1.04902697E-02 1.10472650E-05    3
-2.81561854E-08 1.36583958E-11 7.39219877E+03 1.46022098E+01 8.55492386E+03    4
HNO3              L 4/90H   1N   1O   3    0G   200.000  6000.00  1000.0       1
 0.80037397E+01 0.44984461E-02-0.17365219E-05 0.29369198E-09-0.18148285E-13    2
-0.19256280E+05-0.16098258E+02 0.17449337E+01 0.18804057E-01-0.81595875E-05    3
-0.57859036E-08 0.44377077E-11-0.17380530E+05 0.16954532E+02-0.16105924E+05    4
N2H3              T 7/93H   3N   2    0    0G   200.000  6000.00  1000.0       1
 0.43414654E+01 0.75280979E-02-0.27478351E-05 0.44688178E-09-0.26846990E-13    2
 0.25176779E+05 0.98835045E+00 0.33151120E+01 0.21514763E-02 0.21849694E-04    3
-0.29813376E-07 0.12038856E-10 0.25844190E+05 0.82263324E+01 0.27068024E+05    4
N2H4 HYDRAZINE    L 5/90N   2H   4    0    0G   200.000  6000.00  1000.0       1
 4.93957357E+00 8.75017187E-03-2.99399058E-06 4.67278418E-10-2.73068599E-14    2
 9.28265548E+03-2.69439772E+00 3.83472149E+00-6.49129555E-04 3.76848463E-05    3
-5.00709182E-08 2.03362064E-11 1.00893925E+04 5.75272030E+00 1.14474575E+04    4
CNN               RUS 79C   1N   2    0    0G   200.000  6000.00  1000.0       1
 0.41398983E+01 0.38071002E-02-0.14753456E-05 0.24441991E-09-0.14746300E-13    2
 0.46790796E+05 0.32444306E+01 0.27584988E+01 0.12901042E-01-0.22802003E-04    3
 0.21393697E-07-0.75499090E-11 0.46953824E+05 0.91902188E+01 0.48186940E+05    4
HCNN              SRI/94C   1N   2H   10   0G   300.000  5000.000  1000.000    1
 0.58946362E+01 0.39895959E-02-0.15982380E-05 0.29249395E-09-0.20094686E-13    2
 0.53452941E+05-0.51030502E+01 0.25243194E+01 0.15960619E-01-0.18816354E-04    3
 0.12125540E-07-0.32357378E-11 0.54261984E+05 0.11675870E+02                   4
HCCOH             T 4/93C   2H   2O   1    0G   200.000  6000.000  1000.0      1
 0.63660255E+01 0.55038729E-02-0.18851901E-05 0.29446414E-09-0.17218598E-13    2
 0.89184965E+04-0.82504705E+01 0.19654173E+01 0.25585205E-01-0.38773334E-04    3
 0.31566335E-07-0.10081670E-10 0.97694090E+04 0.12602749E+02 0.11207642E+05    4
N2O4              RUS 89N   2O   4    0    0G   200.000  6000.000  1000.0      1
 1.15752899E+01 4.01616086E-03-1.57178323E-06 2.68274309E-10-1.66922019E-14    2
-2.92191226E+03-3.19488439E+01 3.02002308E+00 2.95904321E-02-3.01342458E-05    3
 1.42360407E-08-2.44100049E-12-6.40040162E+02 1.18059606E+01 1.33632866E+03    4
N2O3              L 4/90N   2O   3    0    0G   200.000  6000.000  1000.0      1
 9.08583845E+00 3.37756330E-03-1.31583890E-06 2.30762329E-10-1.47151267E-14    2
 7.27160146E+03-1.55361904E+01 5.81083964E+00 1.43330962E-02-1.96208597E-05    3
 1.73060735E-08-6.46553954E-12 8.19184453E+03 1.20461321E+00 1.04192062E+04    4
N2O5              L 4/90N   2O   5    0    0G   200.000  6000.000  1000.0      1
 1.31108082E+01 4.87435791E-03-1.87548389E-06 3.16374121E-10-1.95926845E-14    2
-3.11634700E+03-3.46877692E+01 3.68767444E+00 3.92120798E-02-5.53770029E-05    3
 4.20097833E-08-1.31260710E-11-8.30291184E+02 1.21967866E+01 1.59961321E+03    4
NH2OH  WRA032798        N   1H   3O   1    0G   200.000  6000.000  1000.0      1
 3.98241375E+00 7.99825642E-03-2.74883544E-06 4.22874218E-10-2.42498273E-14    2
-6.44279418E+03 3.22666600E+00 2.67285464E+00 1.13645347E-02-4.92179546E-06    3
-9.18041765E-11 6.06669407E-13-6.08956846E+03 1.00068112E+01-4.83091791E+00    4
HNOH              102290H   2N   1O   1     G  0300.00   4000.00  1500.00      1
 0.06396134E+02 0.01821067E-01-0.01870891E-05-0.07844471E-09 0.14448555E-13    2
 0.07859615E+05-0.10404785E+02 0.02125274E+02 0.10662818E-01-0.07602588E-04    3
 0.03081641E-07-0.05726498E-11 0.09553544E+05 0.13096718E+02                   4
H2NO              102290H   2N   1O   1     G  0300.00   4000.00  1500.00      1
 0.05673346E+02 0.02298836E-01-0.01774445E-05-0.11034818E-09 0.01859762E-12    2
 0.05569325E+05-0.06153540E+02 0.02530589E+02 0.08596035E-01-0.05471030E-04    3
 0.02276249E-07-0.04648073E-11 0.06868030E+05 0.11266506E+02                   4
HNNO              103190H   1N   2O   1     G  0300.00   4000.00  1500.00      1
 0.06991217E+02 0.01875970E-01-0.02124584E-05-0.06710472E-09 0.12305080E-13    2
 0.02497566E+06-0.11235229E+02 0.02238298E+02 0.13591997E-01-0.11798728E-04    3
 0.05392970E-07-0.10108589E-11 0.02660258E+06 0.14136789E+02                   4
C2H5CHO           T 9/92C   3H   6O   1    0G   273.150  5000.000 1500.00      1
 0.33137982E+01 0.26619606E-01-0.10475596E-04 0.18815334E-08-0.12761310E-12    2
-0.25459603E+05 0.96608447E+01 0.76044596E+01-0.86403564E-02 0.73930097E-04    3
-0.79687398E-07 0.28004927E-10-0.25489789E+05-0.67643691E+01-0.23097645E+05    4
C2H5CO            T 9/92C   3H   5O   1    0G   298.150  5000.000 1500.00      1
 0.30445698E+01 0.23236429E-01-0.86317936E-05 0.14799550E-08-0.96860829E-13    2
-0.61787211E+04 0.13122302E+02 0.67368294E+01-0.26945299E-02 0.49927017E-04    3
-0.50025808E-07 0.15011503E-10-0.65703366E+04-0.23398732E+01-0.43321855E+04    4
!TC4H7             A 8/83C   4H	 7    0    0G   300.     3000.	  1500.        1
! 0.4219753E+01  0.2882451E-01 -0.1399213E-04  0.3340718E-08 -0.3226427E-12     2
! 0.1266295E+05  0.3253111E+01 -0.2152314E+01  0.5547424E-01 -0.6226715E-04     3
! 0.4593056E-07 -0.1492297E-10  0.1407443E+05  0.3421103E+02  0.1543086E+05     4
H2CNO H2C*N=O     T 9/96H   2C   1N   1O   1G   200.000  6000.000 1500.        1
 0.54028152E+01 0.69057001E-02-0.25162977E-05 0.41014066E-09-0.24718300E-13    2
 0.24528690E+05-0.44574262E+01 0.38781858E+01-0.66530886E-02 0.53947610E-04    3
-0.68176813E-07 0.27181746E-10 0.25716857E+05 0.74618774E+01 0.26932156E+05    4
CH3NO             T12/92C   1H   3N   1O   1G   200.000  6000.000 1500.        1
 0.50677397E+01 0.93871079E-02-0.33958317E-05 0.55076729E-09-0.33095301E-13    2
 0.71852464E+04-0.10709779E+01 0.52463494E+01-0.68175691E-02 0.46713959E-04    3
-0.53482743E-07 0.19916692E-10 0.79241319E+04 0.18687355E+01 0.95017371E+04    4
HCNH cis          T05/97H   2C   1N   1    0G   200.000  6000.000 1500.        1
 4.21964804E+00 5.00385006E-03-1.76392053E-06 2.80725924E-10-1.65851919E-14    2
 3.67706419E+04 1.67138658E+00 3.68324269E+00-1.38553482E-03 2.40042191E-05    3
-3.11573905E-08 1.25791818E-11 3.72527355E+04 6.21248890E+00 3.84457533E+04    4
END
END
END
! Thermodynamic data in this file are mainly from :
!
! Alexander Burcat and Bonnie McBride
! "1997 Ideal Gas Thermodynamic Data for Combustion and Air-Pollution Use"
! Technion Aerospace Engineering (TAE) Report # 804 June 1997.
! or
!
! Alexander Burcat's Ideal Gas Thermochemical Database
! http://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics
! October 1999
!
! Anderson,W.R. Comb.Flame, 1999, v.117, p.394   for:
! HNO
! NH2OH
!
! For other species not found in the above reference data are taken from :
!
! CHEMKIN Thermodynamic Database (1997) for :
! H2C4O
! HCNO
! HOCN
! HNOH
! H2NO
! HNNO
! HCNN
!
! P.Dagaut and J.-C.Boettner (CNRS, LCSR, Orleans, France) for : 
! C2H5O2
! IC4H3
! NC4H3
!
! These species were added 27/1/94 by A.Konnov (Polynomial fits by 
! J.-C.Boettner)
!
! CH3CO3H
! CH3CO3
! CH3CO2H
! CH3CO2
! IC3H7O2H
! IC3H7O2
! NC3H7O2H
! NC3H7O2
! O2C3H6OH
! C3H6OH
! C3H5O2H
! C3H5O2
! C3H5O
! IC3H7O
! CH3HCOW
! CH2CHOW
!
!
! Updated 11/11/99 19/03/00 20/06/00
END

\\
\\
\\  This is the tran file
\\
\\
! edited 20/06/00
AR                 0   136.500     3.330     0.000     0.000     0.000
C                  0    71.400     3.298     0.000     0.000     0.000  !(*)
C2                 1    97.530     3.621     0.000     1.760     4.000
C2O                1   232.400     3.828     0.000     0.000     1.000 ! *
CN2                1   232.400     3.828     0.000     0.000     1.000 ! OIS
C2H                1   209.000     4.100     0.000     0.000     2.500
C2H2               1   209.000     4.100     0.000     0.000     2.500
C2H2OH             2   224.700     4.162     0.000     0.000     1.000  !(*)
CH2OH              2   417.000     3.690     1.700     0.000     2.000
C2H3               2   209.000     4.100     0.000     0.000     1.000  !(*)
C2H4               2   280.800     3.971     0.000     0.000     1.500
C2H5               2   252.300     4.302     0.000     0.000     1.500
C2H6               2   252.300     4.302     0.000     0.000     1.500
C2N                1   232.400     3.828     0.000     0.000     1.000 ! OIS
C2N2               1   349.000     4.361     0.000     0.000     1.000  !(OIS)
C3H2               2   209.000     4.100     0.000     0.000     1.000  !(*)
C3H3               1   252.000     4.760     0.000     0.000     1.000
C3H4               1   252.000     4.760     0.000     0.000     1.000
C3H4P              1   252.000     4.760     0.000     0.000     1.000 ! JAM
PC3H4              1   252.000     4.760     0.000     0.000     1.000 ! JAM
C3H6               2   266.800     4.982     0.000     0.000     1.000
C3H7               2   266.800     4.982     0.000     0.000     1.000
C4H6               2   357.000     5.180     0.000     0.000     1.000
IC3H7              2   266.800     4.982     0.000     0.000     1.000
NC3H7              2   266.800     4.982     0.000     0.000     1.000
C3H8               2   266.800     4.982     0.000     0.000     1.000
C4H                1   357.000     5.180     0.000     0.000     1.000
C4H2               1   357.000     5.180     0.000     0.000     1.000
C4H2OH             2   224.700     4.162     0.000     0.000     1.000  !(*)
C4H3               1   357.000     5.180     0.000     0.000     1.000
C4H4               1   357.000     5.180     0.000     0.000     1.000
C4H8               2   357.000     5.176     0.000     0.000     1.000
T2C4H8             2   357.000     5.176     0.000     0.000     1.000  !=C4H8
C2C4H8             2   357.000     5.176     0.000     0.000     1.000  !=C4H8
NC3H7O             2   357.000     5.176     0.000     0.000     1.000  !=C4H8
IC3H7O             2   357.000     5.176     0.000     0.000     1.000  !=C4H8
IC3H7O2H           2   357.000     5.176     0.000     0.000     1.000  !=C4H8
IC3H7O2            2   357.000     5.176     0.000     0.000     1.000  !=C4H8
NC3H7O2H           2   357.000     5.176     0.000     0.000     1.000  !=C4H8
NC3H7O2            2   357.000     5.176     0.000     0.000     1.000  !=C4H8
IC4H7              2   357.000     5.176     0.000     0.000     1.000  !=C4H8
TC4H7              2   357.000     5.176     0.000     0.000     1.000  !=C4H8
C3H5O2             2   357.000     5.176     0.000     0.000     1.000  !=C4H8
C3H5O2H            2   357.000     5.176     0.000     0.000     1.000  !=C4H8
C3H5O              2   357.000     5.176     0.000     0.000     1.000  !=C4H8
C3H6OH             2   357.000     5.176     0.000     0.000     1.000  !=C4H8
O2C3H6OH           2   357.000     5.176     0.000     0.000     1.000  !=C4H8
CH3CO2H            2   357.000     5.176     0.000     0.000     1.000  !=C4H8
CH3CO2             2   357.000     5.176     0.000     0.000     1.000  !=C4H8
CH3CO3             2   357.000     5.176     0.000     0.000     1.000  !=C4H8
CH3CO3H            2   357.000     5.176     0.000     0.000     1.000  !=C4H8
C2H5CHO            2   357.000     5.176     0.000     0.000     1.000  !=C4H8
C2H5CO             2   357.000     5.176     0.000     0.000     1.000  !=C4H8
C4H9               2   357.000     5.176     0.000     0.000     1.000
S*C4H9             2   357.000     5.176     0.000     0.000     1.000
C4H9               2   357.000     5.176     0.000     0.000     1.000
I*C4H9             2   357.000     5.176     0.000     0.000     1.000
C5H2               1   357.000     5.180     0.000     0.000     1.000
C5H3               1   357.000     5.180     0.000     0.000     1.000
C6H2               1   357.000     5.180     0.000     0.000     1.000
C6H5               2   412.300     5.349     0.000     0.000     1.000 ! JAM
C6H5(L)            2   412.300     5.349     0.000     0.000     1.000 ! JAM
C6H5O              2   450.000     5.500     0.000     0.000     1.000 ! JAM
C5H5OH             2   450.000     5.500     0.000     0.000     1.000 ! JAM
C6H6               2   412.300     5.349     0.000     0.000     1.000 ! SVE
C6H7               2   412.300     5.349     0.000     0.000     1.000 ! JAM
CH                 1    80.000     2.750     0.000     0.000     0.000
CH2                1   144.000     3.800     0.000     0.000     0.000
CH2(S)             1   144.000     3.800     0.000     0.000     0.000
CH2(SING)          1   144.000     3.800     0.000     0.000     0.000
SCH2               1   144.000     3.800     0.000     0.000     0.000
CH2CHCCH           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCCH2          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCH2           2   260.000     4.850     0.000     0.000     1.000 ! JAM
C3H5               2   260.000     4.850     0.000     0.000     1.000 ! JAM 
CH2CHCHCH          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCHCH2         2   357.000     5.180     0.000     0.000     1.000 ! JAM
IC4H5              2   357.000     5.180     0.000     0.000     1.000 
NC4H5              2   357.000     5.180     0.000     0.000     1.000 
CH2CO              2   436.000     3.970     0.000     0.000     2.000
CH2O               2   498.000     3.590     0.000     0.000     2.000
CH3                1   144.000     3.800     0.000     0.000     0.000
CH3CC              2   252.000     4.760     0.000     0.000     1.000 ! JAM
CH3CCCH2           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3CCCH3           2   357.000     5.180     0.000     0.000     1.000 ! JAM
IC4H3              2   357.000     5.180     0.000     0.000     1.000 
NC4H3              2   357.000     5.180     0.000     0.000     1.000 
CH3CCH2            2   260.000     4.850     0.000     0.000     1.000 ! JAM
TC3H5              2   260.000     4.850     0.000     0.000     1.000 ! JAM 
CH3CHCH            2   260.000     4.850     0.000     0.000     1.000 ! JAM
SC3H5              2   260.000     4.850     0.000     0.000     1.000 ! JAM   
CH3CH2CCH          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3HCO             2   436.000     3.970     0.000     0.000     2.000
CH3HCOW            2   436.000     3.970     0.000     0.000     2.000 !=AA
C2H4O              2   436.000     3.970     0.000     0.000     2.000 !=AA
C2H3O              2   436.000     3.970     0.000     0.000     2.000 !=AA
C3H6O              2   436.000     3.970     0.000     0.000     2.000 !=AA
CH3CO              2   436.000     3.970     0.000     0.000     2.000
CH3O               2   417.000     3.690     1.700     0.000     2.000
CH3OH              2   481.800     3.626     0.000     0.000     1.000 ! SVE
CH4                2   141.400     3.746     0.000     2.600    13.000
CH4O               2   417.000     3.690     1.700     0.000     2.000
CN                 1    75.000     3.856     0.000     0.000     1.000  !(OIS)
CNC                1   232.400     3.828     0.000     0.000     1.000 ! OIS
CNN                1   232.400     3.828     0.000     0.000     1.000 ! OIS
HCNN               1   232.400     3.828     0.000     0.000     1.000 !=CNN
CO                 1    98.100     3.650     0.000     1.950     1.800
CO2                1   244.000     3.763     0.000     2.650     2.100
H                  0   145.000     2.050     0.000     0.000     0.000
H2C4O              2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2                 1    38.000     2.920     0.000     0.790   280.000
H2CCCCH            2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2CCCCH2           2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2CCCH             2   252.000     4.760     0.000     0.000     1.000 ! JAM
H2CN               1   569.000     3.630     0.000     0.000     1.000 ! os/jm
H2NO               2   116.700     3.492     0.000     0.000     1.000 ! JAM
H2O                2   572.400     2.605     1.844     0.000     4.000
H2O2               2   107.400     3.458     0.000     0.000     3.800
H2S                2   301.000     3.600     0.000     0.000     1.000  !(OIS)
HC2N2              1   349.000     4.361     0.000     0.000     1.000 ! OIS
HCCHCCH            2   357.000     5.180     0.000     0.000     1.000 ! JAM
HCCO               2   150.000     2.500     0.000     0.000     1.000  !(*)
HCCOH              2   436.000     3.970     0.000     0.000     2.000
HCN                1   569.000     3.630     0.000     0.000     1.000  !(OIS)
HCNO               2   232.400     3.828     0.000     0.000     1.000 ! JAM
H2CNO              2   232.400     3.828     0.000     0.000     1.000 !=HCNO
CH3NO              2   232.400     3.828     0.000     0.000     1.000 !=HCNO
HCNH               2   232.400     3.828     0.000     0.000     1.000 !=HCNO
HCO                2   498.000     3.590     0.000     0.000     0.000
HCO+               1   498.000     3.590     0.000     0.000     0.000
CH2HCO             2   436.000     3.970     0.000     0.000     2.000
CH2CHOW            2   436.000     3.970     0.000     0.000     2.000 !=CH2HCO
HOCN               2   232.400     3.828     0.000     0.000     1.000  !(JAM)
HNCO               2   232.400     3.828     0.000     0.000     1.000  !(OIS)
HNNO               2   232.400     3.828     0.000     0.000     1.000  !(*)
HONO               2   232.400     3.828     0.000     0.000     1.000  !=HNNO
HNO                2   116.700     3.492     0.000     0.000     1.000  !(*)
HNOH               2   116.700     3.492     0.000     0.000     1.000 ! JAM
NH2OH              2   116.700     3.492     0.000     0.000     1.000 !=HNOH
HO2                2   107.400     3.458     0.000     0.000     1.000  !(*)
HSO2               2   252.000     4.290     0.000     0.000     1.000  !(OIS)
HE                 0    10.200     2.576     0.000     0.000     0.000  !(*)
N                  0    71.400     3.298     0.000     0.000     0.000  !(*)
N2                 1    97.530     3.621     0.000     1.760     4.000
N2H2               2    71.400     3.798     0.000     0.000     1.000  !(*)
N2H3               2   200.000     3.900     0.000     0.000     1.000  !(*)
N2H4               2   205.000     4.230     0.000     4.260     1.500
N2O                1   232.400     3.828     0.000     0.000     1.000  !(*)
NCN                1   232.400     3.828     0.000     0.000     1.000 ! OIS
NCO                1   232.400     3.828     0.000     0.000     1.000  !(OIS)
NH                 1    80.000     2.650     0.000     0.000     4.000
NH2                2    80.000     2.650     0.000     2.260     4.000
NH3                2   481.000     2.920     1.470     0.000    10.000
NNH                2    71.400     3.798     0.000     0.000     1.000  !(*)
NO                 1    97.530     3.621     0.000     1.760     4.000
NCNO               2   232.400     3.828     0.000     0.000     1.000 ! OIS
NO2                2   200.000     3.500     0.000     0.000     1.000  !(*)
CH3O2              2   200.000     3.500     0.000     0.000     1.000  !=NO2
CH3O2H             2   200.000     3.500     0.000     0.000     1.000  !=NO2
SO2                2   252.000     4.290     0.000     0.000     1.000  !(OIS)
SO3                2   378.400     4.175     0.000     0.000     1.000  !(OIS)
NO3                2   378.400     4.175     0.000     0.000     1.000  !=SO3
HNO3               2   378.400     4.175     0.000     0.000     1.000  !=SO3
O                  0    80.000     2.750     0.000     0.000     0.000
O2                 1   107.400     3.458     0.000     1.600     3.800
O3                 2   180.000     4.100     0.000     0.000     2.000
OH                 1    80.000     2.750     0.000     0.000     0.000
C2H5OH             2   362.6       4.53      0.000     0.000     1.000  !(sve)
PC2H5O             2   362.6       4.53      0.000     0.000     1.000  !=C2H5OH
SC2H5O             2   362.6       4.53      0.000     0.000     1.000  !=C2H5OH
C2H5O              2   362.6       4.53      0.000     0.000     1.000  !=C2H5OH
C2H5O2H            2   362.6       4.53      0.000     0.000     1.000  !=C2H5OH
C2H5O2             2   362.6       4.53      0.000     0.000     1.000  !=C2H5OH
C2F4               2   202.6       5.164     0.000     0.000     1.000  !(mec) 
N2O4               2   202.6       5.164     0.000     0.000     1.000  !=C2F4
N2O3               2   202.6       5.164     0.000     0.000     1.000  !=C2F4


#endif
