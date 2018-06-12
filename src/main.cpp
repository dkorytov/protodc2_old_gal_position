
#include "dtk/all.hpp"
#include <mpi.h>
#include "GenericIO.h"



std::string output;
std::string gio_output,gio_pos_output;
int step;
float *host_halo_mass, *infall_halo_mass, *host_halo_x, *host_halo_y, *host_halo_z;
float *x, *y, *z, *vx, *vy,*vz;
double* sdss_u, *sdss_g, *sdss_r, *sdss_i, *sdss_z,*stellar_mass;
int64_t* host_halo_tag, *infall_halo_tag, *nodeIndex, *hostIndex;
int* pos_type;

void read_param(char* param_loc){
  std::cout<<"reading params..."<<std::endl;
  dtk::AutoTimer t;
  dtk::Param param(param_loc);
  step = param.get<int>("step");
  output = param.get<std::string>("output");
  output = dtk::rep_str(output,"${step}",step);
  gio_output = dtk::rep_str(output,".hdf5",".gio");
  gio_pos_output = dtk::rep_str(output,".hdf5","_pos.gio");
  std::cout<<"\t"<<gio_output<<" "<<gio_pos_output<<std::endl;
  std::cout<<"\tdone. "<<t<<std::endl;
}
void alloc_arrays(size_t size){
  std::cout<<"allocing arrays..."<<std::endl;
  dtk::AutoTimer t;
  host_halo_tag = new int64_t[size];
  host_halo_mass = new float[size];
  host_halo_x = new float[size];
  host_halo_y = new float[size];
  host_halo_z = new float[size];
  infall_halo_tag = new int64_t[size];
  infall_halo_mass = new float[size];
  x = new float[size];
  y = new float[size];
  z = new float[size];
  vx = new float[size];
  vy = new float[size];
  vz = new float[size];
  sdss_u = new double[size];
  sdss_g = new double[size];
  sdss_r = new double[size];
  sdss_i = new double[size];
  sdss_z = new double[size];
  stellar_mass = new double[size];
  nodeIndex = new int64_t[size];
  hostIndex = new int64_t[size];
  pos_type = new int[size];
  std::cout<<"\tdone. "<<t<<std::endl;
}
void load_data(H5::H5File& hfile,std::string base_name,size_t* gal_size, int num){
  std::cout<<"loading data for "<<base_name<<" "<<num<<std::endl;
  dtk::AutoTimer t;
  size_t offset = 0;
  for(int i =0;i<num;++i)
    offset +=gal_size[i];
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/host_halo_tag"),host_halo_tag+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/host_halo_mass"),host_halo_mass+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/host_halo_x"),host_halo_x+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/host_halo_y"),host_halo_y+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/host_halo_z"),host_halo_z+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/infall_halo_tag"),infall_halo_tag+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/infall_halo_mass"),infall_halo_mass+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/nodeIndex"),nodeIndex+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/hostIndex"),hostIndex+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/x"),x+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/y"),y+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/z"),z+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/vx"),vx+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/vy"),vy+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/vz"),vz+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/sdss_u"),sdss_u+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/sdss_g"),sdss_g+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/sdss_r"),sdss_r+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/sdss_i"),sdss_i+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/sdss_z"),sdss_z+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/sdss_u"),sdss_u+offset);
  dtk::read_hdf5(hfile,dtk::make_str(base_name,"/stellar_mass"),stellar_mass+offset);
  std::fill(pos_type+offset,pos_type+offset+gal_size[num],num);
  std::cout<<"\tdone. "<<t<<std::endl;
}
void write(std::string loc,size_t size,bool all_data=true){
  std::cout<<"writing out data"<<std::endl;
  dtk::AutoTimer t;
  gio::GenericIO gio(MPI_COMM_WORLD, loc);
  gio.setNumElems(size);
  gio.addVariable("x",x,gio::GenericIO::VarIsPhysCoordX);
  gio.addVariable("y",y,gio::GenericIO::VarIsPhysCoordY);
  gio.addVariable("z",z,gio::GenericIO::VarIsPhysCoordZ);
  gio.addVariable("vx",vx);
  gio.addVariable("vy",vy);
  gio.addVariable("vz",vz);
  gio.addVariable("nodeIndex",nodeIndex);
  if(all_data){
    gio.addVariable("host_halo_mass",host_halo_mass);
    gio.addVariable("host_halo_tag",host_halo_tag);
    gio.addVariable("host_halo_x",host_halo_x);
    gio.addVariable("host_halo_y",host_halo_y);
    gio.addVariable("host_halo_z",host_halo_z);
    gio.addVariable("infall_halo_mass",infall_halo_mass);
    gio.addVariable("infall_halo_tag",infall_halo_tag);
    gio.addVariable("sdss_u",sdss_u);
    gio.addVariable("sdss_g",sdss_g);
    gio.addVariable("sdss_r",sdss_r);
    gio.addVariable("sdss_i",sdss_i);
    gio.addVariable("sdss_z",sdss_z);
    gio.addVariable("hostIndex",hostIndex);
    gio.addVariable("positionType",pos_type);
    gio.addVariable("stellar_mass",stellar_mass);
  }
  gio.write();
  std::cout<<"\tdone. "<<t<<std::endl;
}
int main(int argc, char** argv){
  MPI_Init(&argc,&argv);
  if(argc<=1){
    std::cout<<"need to give param file"<<std::endl;
    exit(-1);
  }
  
  read_param(argv[1]);
  H5::H5File hfile(output, H5F_ACC_RDONLY);
  size_t gal_size[6];
  size_t gal_total_size;

  gal_size[0] = dtk::hdf5_num_elems(hfile,"gal_bhp/host_halo_mass");
  gal_size[1] = dtk::hdf5_num_elems(hfile,"gal_rnd/host_halo_mass");
  gal_size[2] = dtk::hdf5_num_elems(hfile,"gal_core/host_halo_mass");
  gal_size[3] = dtk::hdf5_num_elems(hfile,"gal_blank/host_halo_mass");
  gal_size[4] = dtk::hdf5_num_elems(hfile,"gal_nfw/host_halo_mass");
  gal_size[5] = dtk::hdf5_num_elems(hfile,"gal_central/host_halo_mass");
  gal_total_size = gal_size[0]+gal_size[1]+gal_size[2]+gal_size[3]+gal_size[4]+gal_size[5];
  alloc_arrays(gal_total_size);
  load_data(hfile,"gal_bhp",gal_size,0);
  load_data(hfile,"gal_rnd",gal_size,1);
  load_data(hfile,"gal_core",gal_size,2);
  load_data(hfile,"gal_blank",gal_size,3);
  load_data(hfile,"gal_nfw",gal_size,4);
  load_data(hfile,"gal_central",gal_size,5);
  write(gio_output,gal_total_size,true);
  write(gio_pos_output,gal_total_size,false);
  MPI_Finalize();
}
