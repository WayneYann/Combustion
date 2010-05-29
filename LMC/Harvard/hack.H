      integer HACKF_L1, HACKF_L2, HACKF_H1, HACKF_H2, HACK_NC
      integer HACKC_L1, HACKC_L2, HACKC_H1, HACKC_H2
      parameter (HACKF_L1=0, HACKF_L2=0, HACKF_H1=31, HACKF_H2=95)
      parameter (HACKC_L1=0, HACKC_L2=0, HACKC_H1=15, HACKC_H2=47)
      parameter (HACK_NC=32)
      integer hackflo(2), hackfhi(2)
      integer hackclo(2), hackchi(2)
      REAL_T hackfdat(HACKF_L1:HACKF_H1, HACKF_L2:HACKF_H2, HACK_NC)
      REAL_T hackcdat(HACKC_L1:HACKC_H1, HACKC_L2:HACKC_H2, HACK_NC)
      common / hackcomf / hackfdat, hackflo, hackfhi
      common / hackcomc / hackcdat, hackclo, hackchi
      save /hackcomf/, /hackcomc/
      
