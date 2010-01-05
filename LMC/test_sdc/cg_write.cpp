#include <Array.H>
#include <FArrayBox.H>
#include <fstream.h>

#include "cg_write.H"
#include "cg_write_F.H"

void cg_write(MultiFab& mf, char* fname, const int ix)	   
{
  std::ofstream outFile(fname);
  const int ncomp = mf.nComp();


  outFile<<"# choosing x-slice i = "<<ix<<std::endl;

  for (MFIter Smfi(mf); Smfi.isValid(); ++Smfi)
    {
      FArrayBox& fb = mf[Smfi];
      const Box& bx = Smfi.validbox();


      FORT_CG_WRITE(bx.loVect(), bx.hiVect(), &ncomp,
		    fb.dataPtr(), ARLIM(fb.loVect()), ARLIM(fb.hiVect()),
		    fname, &ix);
       
    }
  outFile.close();
}
