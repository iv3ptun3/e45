#include "getenv.hh"
#include "globals.hh"

void getenv_read(){

  if(getenv("Truncated_mean_cut")){
    G4String env_truncated_mean = getenv("Truncated_mean_cut");
    env.truncated_mean = atof( env_truncated_mean.c_str() );
    std::cout<<"!======================================!"<<std::endl;
    std::cout<<"!======================================!"<<std::endl;
    std::cout<<"!======================================!"<<std::endl;
    std::cout<<"!======================================!"<<std::endl;
    std::cout<<"read getenv from shell variables"<<std::endl;
    std::cout<<env.truncated_mean<<std::endl;
    std::cout<<"!======================================!"<<std::endl;
    std::cout<<"!======================================!"<<std::endl;
    std::cout<<"!======================================!"<<std::endl;
    std::cout<<"!======================================!"<<std::endl;
  }else{
    std::cout<<"!======================================!"<<std::endl;
    std::cout<<"getenv error: Truncated_mean_cut"<<std::endl;
    std::cout<<"!======================================!"<<std::endl;
    //    break;
  }
}
