#include <winstd.H>

#include "LevelBld.H"

#include "CNSReact.H"
#include "CNSReact_F.H"
#include "Derive_F.H"

using std::string;

static Box the_same_box (const Box& b) { return b; }
static Box grow_box_by_one (const Box& b) { return BoxLib::grow(b,1); }

typedef StateDescriptor::BndryFunc BndryFunc;

//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall
//
static int scalar_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static int norm_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD,  REFLECT_ODD,  REFLECT_ODD
};

static int tang_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
void
set_scalar_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i,scalar_bc[lo_bc[i]]);
        bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}

static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_vel_bc[lo_bc[0]]);
    bc.setHi(0,norm_vel_bc[hi_bc[0]]);
#if (BL_SPACEDIM >= 2)
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

#if (BL_SPACEDIM >= 2)
static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,norm_vel_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}
#endif

#if (BL_SPACEDIM == 3)
static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_bc[hi_bc[2]]);
}
#endif

void
CNSReact::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    // Initialize the network
    network_init();

    // Get options, set phys_bc
    read_params();
    //
    // Set number of state variables and pointers to components
    //
    int cnt = 0;
    Density = cnt++;
    Xmom = cnt++;
#if (BL_SPACEDIM >= 2)
    Ymom = cnt++;
#endif
#if (BL_SPACEDIM == 3)
    Zmom = cnt++;
#endif
    Eden = cnt++;
    Eint = cnt++;
    Temp = cnt++;
    NumAdv = 0;
    if (NumAdv > 0)
    {
        FirstAdv = cnt++;
        cnt += NumAdv - 2;
        LastAdv = cnt++;
    }

    int dm = BL_SPACEDIM;

    // Get the number of species from the network model.
    BL_FORT_PROC_CALL(GET_NUM_SPEC, get_num_spec)(&NumSpec);

    if (NumSpec > 0)
    {
        FirstSpec = cnt++;
        cnt += NumSpec - 2;
        LastSpec = cnt++;
    }

    // Get the number of species from the network model.
    BL_FORT_PROC_CALL(GET_NUM_AUX, get_num_aux)(&NumAux);

    if (NumAux > 0)
    {
        FirstAux = cnt++;
        cnt += NumAux - 2;
        LastAux = cnt++;
    }

    NUM_STATE = cnt;

    // Define NUM_GROW from the f90 module.
    BL_FORT_PROC_CALL(GET_METHOD_PARAMS, get_method_params)(&NUM_GROW);

    const Real run_strt = ParallelDescriptor::second() ; 

    BL_FORT_PROC_CALL(SET_METHOD_PARAMS, set_method_params)
        (dm, Density, Xmom, Eden, Eint, Temp, FirstAdv, FirstSpec, FirstAux, 
         NumAdv, small_dens, small_temp, small_pres, 
         allow_negative_energy,ppm_type,
         gamma,normalize_species);

    Real run_stop = ParallelDescriptor::second() - run_strt;
 
    ParallelDescriptor::ReduceRealMax(run_stop,ParallelDescriptor::IOProcessorNumber());
 
    if (ParallelDescriptor::IOProcessor())
        std::cout << "\nTime in set_method_params: " << run_stop << '\n' ;

    int coord_type = Geometry::Coord();
    BL_FORT_PROC_CALL(SET_PROBLEM_PARAMS, set_problem_params)
         (dm,phys_bc.lo(),phys_bc.hi(),Outflow,Symmetry,coord_type);

    Interpolater* interp = &cell_cons_interp;

    // Note that the default is state_data_extrap = false, store_in_checkpoint = true
    // We only need to put these explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

    store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,1,NUM_STATE,
                           interp,state_data_extrap,store_in_checkpoint);

#ifdef REACTIONS
    // Components 0:Numspec-1         are      omegadot_i
    // Component    NumSpec            is      enuc =      (eout-ein)
    // Component    NumSpec+1          is  rho_enuc= rho * (eout-ein)
    store_in_checkpoint = true;
    store_in_checkpoint = false;
    desc_lst.addDescriptor(Reactions_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,NumSpec+2,
                           &cell_cons_interp,state_data_extrap,store_in_checkpoint);
#endif

    Array<BCRec>       bcs(NUM_STATE);
    Array<std::string> name(NUM_STATE);
    
    BCRec bc;
    cnt = 0;
    set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "density";
    cnt++; set_x_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "xmom";
#if (BL_SPACEDIM >= 2)
    cnt++; set_y_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "ymom";
#endif
#if (BL_SPACEDIM == 3)
    cnt++; set_z_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "zmom";
#endif
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_E";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_e";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "Temp";

    for (int i=0; i<NumAdv; ++i)
    {
        char buf[64];
        sprintf(buf, "adv_%d", i);
        cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = string(buf);
    }

    // Get the species names from the network model.
    char* spec_names[NumSpec];
    for (int i = 0; i < NumSpec; i++) {
          int len = 20;
          Array<int> int_spec_names(len);
          // This call return the actual length of each string in "len" 
          BL_FORT_PROC_CALL(GET_SPEC_NAMES, get_spec_names)(int_spec_names.dataPtr(),&i,&len);
          spec_names[i] = new char[len+1];
          for (int j = 0; j < len; j++) 
             spec_names[i][j] = int_spec_names[j];
          spec_names[i][len] = '\0';
    }

    if ( ParallelDescriptor::IOProcessor())
    {
        std::cout << NumSpec << " Species: " << std::endl;
        for (int i = 0; i < NumSpec; i++)  
           std::cout << spec_names[i] << ' ' << ' ';
        std::cout << std::endl;
    } 

    for (int i=0; i<NumSpec; ++i)
    {
        cnt++; 
        set_scalar_bc(bc,phys_bc); 
        bcs[cnt] = bc; 
        string spec_string(spec_names[i]);
        name[cnt] = "rho_" + spec_string;
    }

    // Get the auxiliary names from the network model.
    char* aux_names[NumAux];
    for (int i = 0; i < NumAux; i++) {
          int len = 20;
          Array<int> int_aux_names(len);
          // This call return the actual length of each string in "len"
          BL_FORT_PROC_CALL(GET_AUX_NAMES, get_aux_names)(int_aux_names.dataPtr(),&i,&len);
          aux_names[i] = new char[len+1];
          for (int j = 0; j < len; j++)
             aux_names[i][j] = int_aux_names[j];
          aux_names[i][len] = '\0';
    }

    if ( ParallelDescriptor::IOProcessor())
    {
        std::cout << NumAux << " Auxiliary Variables: " << std::endl;
        for (int i = 0; i < NumAux; i++)
           std::cout << aux_names[i] << ' ' << ' ';
        std::cout << std::endl;
    }

    for (int i=0; i<NumAux; ++i)
    {
        cnt++;
        set_scalar_bc(bc,phys_bc);
        bcs[cnt] = bc;
        string aux_string(aux_names[i]);
        name[cnt] = "rho_" + aux_string;
    }

    desc_lst.setComponent(State_Type,
                          Density,
                          name,
                          bcs,
                          BndryFunc(BL_FORT_PROC_CALL(CA_DENFILL,ca_denfill),
                                    BL_FORT_PROC_CALL(CA_HYPFILL,ca_hypfill)));

#ifdef REACTIONS
    std::string name_react;
    for (int i=0; i<NumSpec; ++i)
    {
       set_scalar_bc(bc,phys_bc);
       string aux_string(spec_names[i]);
       name_react = "omegadot_" + aux_string;
       desc_lst.setComponent(Reactions_Type, i, name_react, bc,
                             BndryFunc(BL_FORT_PROC_CALL(CA_REACTFILL,ca_reactfill)));
    }
    desc_lst.setComponent(Reactions_Type, NumSpec  , "enuc", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_REACTFILL,ca_reactfill)));
    desc_lst.setComponent(Reactions_Type, NumSpec+1, "rho_enuc", bc,
                          BndryFunc(BL_FORT_PROC_CALL(CA_REACTFILL,ca_reactfill)));
#endif

    //
    // DEFINE DERIVED QUANTITIES
    //
    // Pressure
    //
    derive_lst.add("pressure",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERPRES,ca_derpres),the_same_box);
    derive_lst.addComponent("pressure",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Kinetic energy
    //
    derive_lst.add("kineng",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERKINENG,ca_derkineng),the_same_box);
    derive_lst.addComponent("kineng",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("kineng",desc_lst,State_Type,Xmom,BL_SPACEDIM);

    //
    // Sound speed (c)
    //
    derive_lst.add("soundspeed",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERSOUNDSPEED,ca_dersoundspeed),the_same_box);
    derive_lst.addComponent("soundspeed",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Mach number(M)
    //
    derive_lst.add("MachNumber",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERMACHNUMBER,ca_dermachnumber),the_same_box);
    derive_lst.addComponent("MachNumber",desc_lst,State_Type,Density,NUM_STATE);

#if (BL_SPACEDIM == 1)
    //
    // Wave speed u+c
    //
    derive_lst.add("uplusc",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERUPLUSC,ca_deruplusc),the_same_box);
    derive_lst.addComponent("uplusc",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Wave speed u-c
    //
    derive_lst.add("uminusc",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERUMINUSC,ca_deruminusc),the_same_box);
    derive_lst.addComponent("uminusc",desc_lst,State_Type,Density,NUM_STATE);
#endif

    //
    // Entropy (S)
    //
    derive_lst.add("entropy",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERENTROPY,ca_derentropy),the_same_box);
    derive_lst.addComponent("entropy",desc_lst,State_Type,Density,NUM_STATE);

#if (BL_SPACEDIM > 1)
    //
    // Vorticity
    //
    derive_lst.add("magvort",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERMAGVORT,ca_dermagvort),grow_box_by_one);
    // Here we exploit the fact that Xmom = Density + 1
    //   in order to use the correct interpolation.
    if (Xmom != Density+1)
       BoxLib::Error("We are assuming Xmom = Density + 1 in CNSReact_setup.cpp");
    derive_lst.addComponent("magvort",desc_lst,State_Type,Density,BL_SPACEDIM+1);
#endif

    //
    // Div(u)
    //
    derive_lst.add("divu",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERDIVU,ca_derdivu),grow_box_by_one);
    derive_lst.addComponent("divu",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("divu",desc_lst,State_Type,Xmom,BL_SPACEDIM);

    //
    // Internal energy as derived from rho*E, part of the state
    //
    derive_lst.add("eint_E",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DEREINT1,ca_dereint1),the_same_box);
    derive_lst.addComponent("eint_E",desc_lst,State_Type,Density,NUM_STATE);

    //
    // Internal energy as derived from rho*e, part of the state
    //
    derive_lst.add("eint_e",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DEREINT2,ca_dereint2),the_same_box);
    derive_lst.addComponent("eint_e",desc_lst,State_Type,Density,NUM_STATE);

#if 0
    //
    // Log(density)
    //
    derive_lst.add("logden",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERLOGDEN,ca_derlogden),the_same_box);
    derive_lst.addComponent("logden",desc_lst,State_Type,Density,NUM_STATE);
#endif

    derive_lst.add("StateErr",IndexType::TheCellType(),3,
                   BL_FORT_PROC_CALL(CA_DERSTATE,ca_derstate),grow_box_by_one);
    derive_lst.addComponent("StateErr",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("StateErr",desc_lst,State_Type,Temp,1);
    derive_lst.addComponent("StateErr",desc_lst,State_Type,FirstSpec,1);

    //
    // X from rhoX
    //
    for (int i = 0; i < NumSpec; i++){
      string spec_string(spec_names[i]);
      spec_string = "X("+spec_string+")";

      derive_lst.add(spec_string,IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(CA_DERSPEC,ca_derspec),the_same_box);
      derive_lst.addComponent(spec_string,desc_lst,State_Type,Density,1);
      derive_lst.addComponent(spec_string,desc_lst,State_Type,FirstSpec+i,1);
    }
    //
    // Velocities
    //
    derive_lst.add("x_velocity",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERVEL,ca_dervel),the_same_box);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,Xmom,1);

#if (BL_SPACEDIM >= 2)
    derive_lst.add("y_velocity",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERVEL,ca_dervel),the_same_box);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,Ymom,1);
#endif

#if (BL_SPACEDIM == 3)
    derive_lst.add("z_velocity",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERVEL,ca_dervel),the_same_box);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,Zmom,1);
#endif

    derive_lst.add("magvel",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERMAGVEL,ca_dermagvel),the_same_box);
    derive_lst.addComponent("magvel",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("magvel",desc_lst,State_Type,Xmom,BL_SPACEDIM);

#if (BL_SPACEDIM > 1)
    derive_lst.add("radvel",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERRADIALVEL,ca_derradialvel),the_same_box);
    derive_lst.addComponent("radvel",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("radvel",desc_lst,State_Type,Xmom,BL_SPACEDIM);
#endif

    derive_lst.add("magmom",IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERMAGMOM,ca_dermagmom),the_same_box);
    derive_lst.addComponent("magmom",desc_lst,State_Type,Xmom,BL_SPACEDIM);

    for (int i = 0; i < NumSpec; i++)  {
      string spec_string(spec_names[i]);
      derive_lst.add(spec_string,IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERSPEC,ca_derspec),the_same_box);
      derive_lst.addComponent(spec_string,desc_lst,State_Type,Density,1);
      derive_lst.addComponent(spec_string,desc_lst,State_Type,FirstSpec+i,1);
    }

    for (int i = 0; i < NumAux; i++)  {
      string aux_string(aux_names[i]);
      derive_lst.add(aux_string,IndexType::TheCellType(),1,
          BL_FORT_PROC_CALL(CA_DERSPEC,ca_derspec),the_same_box);
      derive_lst.addComponent(aux_string,desc_lst,State_Type,Density,1);
      derive_lst.addComponent(aux_string,desc_lst,State_Type,FirstAux+i,1);
    }

#if 0
    //
    // A derived quantity equal to all the state variables.
    //
    derive_lst.add("FULLSTATE",IndexType::TheCellType(),NUM_STATE,FORT_DERCOPY,the_same_box);
    derive_lst.addComponent("FULLSTATE",desc_lst,State_Type,Density,NUM_STATE);

#endif

    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //
    ErrorSetUp();
}
