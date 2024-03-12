#define _USE_MATH_DEFINES
#include "prob.H"
#include <cmath>

std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"') {
      count++;
    }
  }
  return (count % 2) == 0;
}

void
read_pmf(const std::string& myfile)
{
  std::string firstline;
  std::string secondline;
  std::string remaininglines;
  unsigned int pos1;
  unsigned int pos2;
  int variable_count;
  int line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline)) {
    amrex::Abort("PMF file variable quotes unbalanced");
  }
  std::getline(iss, secondline);
  pos1 = 0;
  pos2 = 0;
  variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  amrex::Vector<std::string> pmf_names;
  pmf_names.resize(variable_count);
  pos1 = 0;
  // pos2 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    pmf_names[i] = firstline.substr(pos1 + 1, pos2 - (pos1 + 1));
    pos1 = pos2 + 1;
  }

  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;
  for (int i = 0; i < variable_count; i++)
   amrex::Print() << "Variable found: " << pmf_names[i] <<
   std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  PeleC::h_prob_parm_device->pmf_N = line_count;
  PeleC::h_prob_parm_device->pmf_M = variable_count - 1;
  PeleC::prob_parm_host->h_pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->h_pmf_Y.resize(
    PeleC::h_prob_parm_device->pmf_N * PeleC::h_prob_parm_device->pmf_M);
  PeleC::prob_parm_host->pmf_Y.resize(
    PeleC::h_prob_parm_device->pmf_N * PeleC::h_prob_parm_device->pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (unsigned int i = 0; i < PeleC::h_prob_parm_device->pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> PeleC::prob_parm_host->h_pmf_X[i];
    for (unsigned int j = 0; j < PeleC::h_prob_parm_device->pmf_M; j++) {
      sinput >> PeleC::prob_parm_host
                  ->h_pmf_Y[j * PeleC::h_prob_parm_device->pmf_N + i];
    }
  }
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_X.begin(),
    PeleC::prob_parm_host->h_pmf_X.end(), PeleC::prob_parm_host->pmf_X.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_Y.begin(),
    PeleC::prob_parm_host->h_pmf_Y.end(), PeleC::prob_parm_host->pmf_Y.begin());
  PeleC::h_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::h_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
  PeleC::d_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::d_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
}

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("p_l", PeleC::h_prob_parm_device->p_l);
    pp.query("p_r", PeleC::h_prob_parm_device->p_r);
    pp.query("rho_l", PeleC::h_prob_parm_device->rho_l);
    pp.query("rho_r", PeleC::h_prob_parm_device->rho_r);
//    pp.query("u_l", PeleC::h_prob_parm_device->u_l);
    pp.query("u_r", PeleC::h_prob_parm_device->u_r);
//    pp.query("v_l", PeleC::h_prob_parm_device->v_l);
    pp.query("v_r", PeleC::h_prob_parm_device->v_r);
    pp.query("theta", PeleC::h_prob_parm_device->theta);
    
    PeleC::h_prob_parm_device->u_l = PeleC::h_prob_parm_device->V_l*sin(M_PI*PeleC::h_prob_parm_device->theta/180);
    PeleC::h_prob_parm_device->v_l = -1*PeleC::h_prob_parm_device->V_l*cos(M_PI*PeleC::h_prob_parm_device->theta/180);

    std::string znd_datafile;
    pp.query("znd_datafile", znd_datafile);
    read_pmf(znd_datafile);
  }

  auto eos = pele::physics::PhysicsType::eos();
  eos.RYP2T(
    PeleC::h_prob_parm_device->rho_l,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->p_l,
    PeleC::h_prob_parm_device->T_l);
  eos.RYP2E(
    PeleC::h_prob_parm_device->rho_l,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->p_l,
    PeleC::h_prob_parm_device->e_l);

  eos.RYP2T(
    PeleC::h_prob_parm_device->rho_r,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->p_r,
    PeleC::h_prob_parm_device->T_r);
  eos.RYP2E(
    PeleC::h_prob_parm_device->rho_r,
    PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->p_r,
    PeleC::h_prob_parm_device->e_r);
  
  /*amrex::Real cp = 0;
  eos.TY2Cp(
    PeleC::h_prob_parm_device->T_l, PeleC::h_prob_parm_device->massfrac.begin(),
    cp);
  
  auto& trans_parm = PeleC::trans_parms.host_trans_parm();
  amrex::Real Pr = 0.7;
  trans_parm.const_conductivity = trans_parm.const_viscosity * cp / Pr;
  PeleC::trans_parms.sync_to_device();*/
}
}

void
PeleC::problem_post_timestep()
{
 /* //Modify the checkpoint
  amrex::MultiFab& S_new = get_new_data(State_Type);
  amrex::MultiFab& S_old = get_old_data(State_Type);
  
  //MFIter to iterate across the domain
  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();++mfi) 
  {
    const amrex::Box& box = mfi.tilebox();
    auto sfab_new = S_new.array(mfi);
    auto sfab_old = S_old.array(mfi);
    
    const auto geomdata = geom.data();
    
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    
    const ProbParmDevice* lprobparm = d_prob_parm_device;

    if(lprobparm->PostStep>0)
    {
      amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
         //Defined at prob.H, just like pc_initdata
         pc_poststep(i, j, k, sfab_new, sfab_old, geomdata, *lprobparm);
        });
    }
  }*/
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
