C $MEMBER=DILOG, DATE=76022516, USER=EAA
C $MEMBER=DILOG, FUNCTION DILOG/ENTRY DDILOG COMPUTE DILOGARITHM
C $MEMBER=DILOG, IN SINGLE/DOUBLE PRECISION FOR ARGUMENT IN
C $MEMBER=DILOG, SINGLE/DOUBLE PRECISION.
C $MEMBER=DILOG, AUTHOR: E. A. ALLTON   DATE FEB 21 L976
C **** OBTAINED FROM SLAC 28-JUN-82 RC YORK ****
C
      FUNCTION DDILOG(DARG)
      implicit none
      REAL*8 DDILOG,darg
      real*8 v(40),r(40)
      real*8 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10
      real*8 b11,b12,b13,b14,b15,b16,b17,b18,b19,b20
      real*8 b21,b22,b23,b24,b25,b26,b27,b28,b29,b30
      real*8 b31,b32,b33,b34,b35,b36,b37,b38,b39,b40
      real*8 a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      real*8 a11,a12,a13,a14,a15,a16,a17
      real*8 sp3,sp4,sp6,sp8,sp12
      real*8 dp3,dp4,dp6,dp8,dp12
      real*8 d,w,dpi,sph,dpisq,g,b,dsgn,dph,nsw
      integer*4 n,j,key
C
C   DDILOG(DARG) COMPUTES REAL PART OF SPENCE DILOGARITHM
C       FOR ANY REAL VALUE OF ARG.  REAL*8 DDILOG, DARG
C   AUTHOR:  E. A. ALLTON     DATE: APRIL 24, 1975
      EQUIVALENCE
     &  (V(1),B1),(V(2),B2),(V(3),B3),(V(4),B4),(V(5),B5),(V(6),B6),
     &(V(7),B7),(V(8),B8),(V(9),B9),(V(10),B10),(V(11),B11),(V(12),B12),
     &  (V(13),B13),(V(14),B14),(V(15),B15),(V(16),B16),(V(17),B17),
     &  (V(18),B18),(V(19),B19),(V(10),B20),(V(21),B21),(V(22),B22),
     &  (V(23),B23),(V(24),B24),(V(25),B25),(V(26),B26),(V(27),B27),
     &  (V(28),B28),(V(29),B29),(V(30),B30),(V(31),B31),(V(32),B32),
     &  (V(33),B33),(V(34),B34),(V(35),B35),(V(36),B36),(V(37),B37),
     &  (V(38),B38),(V(39),B39),(V(40),B40)
      EQUIVALENCE
     &  (R(1),A1),(R(2),A2),(R(3),A3),(R(4),A4),(R(5),A5),(R(6),A6),
     &(R(7),A7),(R(8),A8),(R(9),A9),(R(10),A10),(R(11),A11),(R(12),A12),
     &  (R(13),A13),(R(14),A14),(R(15),A15),(R(16),A16),(R(17),A17)
C
      LOGICAL LOCK
      DATA LOCK,KEY,N/.TRUE.,0,40/
      DATA DPI/3.1415926535897931D0/
      save
C
      IF(LOCK) GO TO 92
  400 D=DARG
      DSGN=1.D0
C  LOCATE DARG RELATIVE TO SUBINTERVAL BOUNDARY POINTS -1, 0, 1/2, 1, 2.
      IF(D-.5D0) 102, 101, 105
  101 DDILOG=DPH
      RETURN
  102 IF(D) 112, 103, 104
  103 DDILOG=0.D0
      RETURN
  104 W=D
      B=0.D0
      GO TO 120
  105 IF(D-1.D0) 107, 106, 108
  106 DDILOG=DP6
      RETURN
  107 W=1.D0-D
      B=DP6-DLOG(D)*DLOG(W)
      GO TO 119
  108 IF(D-2.D0) 110, 109, 111
  109 DDILOG=DP4
      RETURN
  110 W=1.D0-1.D0/D
      B=DP6-DLOG(D)*DLOG(W*W*D)/2.D0
      GO TO 120
  111 W=1.D0/D
      B=DP3-DLOG(D)**2/2.D0
      GO TO 119
  112 IF(D+1.D0) 114, 113, 115
  113 DDILOG=-DP12
      RETURN
  114 W=1.D0/(1.D0-D)
      B=-DP6+DLOG(W)*DLOG(D*D*W)/2.D0
      GO TO 120
  115 W=-D/(1.D0-D)
      B=-DLOG(-D/W)**2/2.D0
  119 DSGN=-1.D0
  120 NSW=20.D0*W+1.D0
      G=0.D0
      GO TO (125,125,124,123,123,123,123,122,122,121,121), NSW
C  BACKWARD SUMMATION OF SERIES IN FUNDAMENTAL INTERVAL 0<W<1/2.
  121 G=B33+W*(B34+W*(B35+W*(B36+W*(B37+W*(B38+W*(B39+W*B40))))))
  122 G=B25+W*(B26+W*(B27+W*(B28+W*(B29+W*(B30+W*(B31+W*(B32+W*G)))))))
  123 G=B18+W*(B19+W*(B20+W*(B21+W*(B22+W*(B23+W*(B24+W*G))))))
  124 G=B12+W*(B13+W*(B14+W*(B15+W*(B16+W*(B17+W*G)))))
  125 DDILOG=B+DSGN*W*(B1+W*(B2+W*(B3+W*(B4+W*(B5+W*(B6+W*(B7+W*(B8+
     & W*(B9+W*(B10+W*(B11+W*G)))))))))))
      RETURN
C
C     ONE-SHOT COMPUTATION AND STORAGE OF CONSTANTS.
C
   92 KEY=KEY+1
      LOCK=.FALSE.
      DPISQ=DPI**2
      DP3=DPISQ/3.D0
      DP4=DPISQ/4.D0
      DP6=DPISQ/6.D0
      DP8=DPISQ/8.D0
      DP12=DPISQ/12.D0
      DPH=(DP6-DLOG(2.D0)**2)/2.D0
      SP3=DP3
      SP4=DP4
      SP6=DP6
      SP8=DP8
      SP12=DP12
      SPH=DPH
      DO 98  J=1,N
      G=1.D0/DFLOAT(J*J)
      V(J)=G
   98 R(J)=G
      GO TO (400, 777), KEY
  777 RETURN
      END

