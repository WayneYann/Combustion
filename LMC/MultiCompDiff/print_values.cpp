void
HeatTransfer::print_values (const char* message,
			    const int i_coordinate,
			    const int j_coordinate,
			    const int first,
			    const int number,
			    MultiFab * multifab)
{
    // loop over fabs
    for (MFIter mfi(*multifab); mfi.isValid(); ++mfi)
    {
	// get index of the present box
	const int idx = mfi.index();
	const int* lo_vect = mfi.validbox().loVect();
	const int* hi_vect = mfi.validbox().hiVect();
	if (lo_vect[0] <= i_coordinate && i_coordinate <= hi_vect[0] && 
	    lo_vect[1] <= j_coordinate && j_coordinate <= hi_vect[1])
	{
	    std::cout << std::endl;
	    Real values[first+number];
	    (*multifab)[mfi].getVal(values,IntVect(i_coordinate,j_coordinate),first,number);
	    for (int n = first; n < first + number; ++n)
	    {
		std::cout << message
			  << " component " 
			  << n 
			  << " = " 
			  << values[n]
			  << std::endl;
	    }
	    std::cout << std::endl;
	}
    }
}
    
