NAME    METALS   ( FROM IBM MANUAL) ( MIN)
ROWS
  N VALUE
  E YIELD
  L FE
  L MN
  L CU
  L MG
  G AL
  L SI
COLUMNS
    BIN1      VALUE       0.03         YIELD           1.
    BIN1      FE          0.15         CU             .03
    BIN1      MN          0.02         MG             .02
    BIN1      AL          0.7          SI             .02
    BIN2      VALUE       0.08         YIELD           1.
    BIN2      FE           .04         CU             .05
    BIN2      MN           .04         MG             .03
    BIN2      AL           .75         SI             .06
    BIN3      VALUE       0.17         YIELD         1.
    BIN3      FE           .02         CU             .08
    BIN3      MN           .01         AL             .8
    BIN3      SI           .08
    BIN4      VALUE       0.12         YIELD         1.
    BIN4      FE           .04         CU             .02
    BIN4      MN           .02         AL             .75
    BIN4      SI          0.12
    BIN5      VALUE       0.15         YIELD         1.
    BIN5      FE           .02         CU             .06
    BIN5      MN           .02         MG             .01
    BIN5      SI           .02         AL             .8
    ALUM      VALUE       0.21         YIELD         1.
    ALUM      FE           .01         CU             .01
    ALUM      AL           .97         SI             .01
    SILCON    VALUE       0.38         YIELD         1.
    SILCON    FE           .03         SI             .97
RHS
    ALOY1     YIELD        2000.       FE              60.
    ALOY1     CU            100.       MN              40.
    ALOY1     MG             30.       AL            1500.
    ALOY1     SI            300.
BOUNDS
 UP PROD1     BIN1            200.00
 UP PROD1     BIN2            750.00
 LO PROD1     BIN3            400.00
 UP PROD1     BIN3            800.00
 LO PROD1     BIN4            100.00
 UP PROD1     BIN4            700.00
 UP PROD1     BIN5           1500.00
RANGES
    AL1       SI               50.0
ENDATA
