#ifndef INPUTPROCESSOR_H
#define INPUTPROCESSOR_H

#include <string>
#include <vector>
#include <array>
#include "computation.h"
#include <fstream>

using namespace std;

struct AorticInput{
    string aorta_file = "aorta.stl";
    string leaflet_boundaries = "aorta_rec_bnd.bnd";
    array<double, 3> main_flow_direction = {0, 0, 1}; //middle direction of blood flow in aorta
};

struct LeafletInput{
    bool generative = false;
    //if true then leaflet will be created according to model
    //else leaflet will be loaded from leaf_file
    string model;
    std::vector<double> model_params;
    string leaflet_file = "leaflet.txt";
    array<double, 3> main_axis_direction = {0, 1, 0};
    //this information is necessary for sewing with right orientation
    array<double, 3> second_line_symmetry= {0, 0, 1};
    /*
     * The plane formed by the main_axis_direction and second_line_symmetry
     * should divide the pattern so that the extreme points of the
     * pattern are extreme relative to this plane
     */
};

struct SewEnergyParams{
    /*The weights with which these restrictions are taken in the energy functional*/
    double plane_w = 1.0;
    double sqr_length_w = 1.0;
    double digedral_angle_w = 0.7;
    //Additional factor increasing energy from concavities
    //The more factor, the less concavities
    double convexity_w = 10;
    double force_w = 0.0;
    double shift = 0.5;
    struct EnergyMinimizerParams {
        //solver params - sp
        int freq = 50;
        double step_sz = 1.e-4;
        double tol = 1.e-4;
        double epsabs = 1.e-2;
        int maxits = 301;
        double time/*secs*/ = 3;
    };
    EnergyMinimizerParams sp;
};

struct CameraParams{
    array<double, 3> origin = {0.800262, -2.22123, 3.87111};
    double pitch = 33.2;
    double yaw = 9;
};

struct SolverParams{
    double delta = 5.5e-7;                      //coefficient converting force into shift
    double eps = 0.001;                         //stop condition
    double max_possible_shift_scale = 200;
    double max_possible_recommend_shift_scale = 200;
    double max_time = 7200;
};

struct ObjectTraits{
    int elastic_model_type = EMOD_MSM;          //number of used elastic model
    std::vector<double> elastic_model_params;   //coefficients of model
    double thickness = 0.5;
    double pressure = 80.0;                     //pressure acted on object
    double collision_margin = 0.05;             //The distance at which collision detection is turned on
    int renderLabel = 0;
    void add_model_param(double param){ elastic_model_params.push_back(param); }
    ObjectTraits(){
        set_EMOD_MSM_defaults();
    }
private:
    void set_EMOD_MSM_defaults(){
        add_model_param(1350.0   /*Pa*/); //E1_0
        add_model_param(265000.0 /*Pa*/); //E1_1
        add_model_param(0.93           ); //Eps_lim1
        add_model_param(10800.0  /*Pa*/); //E2_0
        add_model_param(1570000.0/*Pa*/); //E2_1
        add_model_param(0.185          ); //Eps_lim2
    }
};

struct RequiredOutput{
    bool use = true;
    string postfix = "";
    string res_dir = "";
    string aorta_to_nts = "";
    bool divide_leaflets = true;
    vector<string> leaf_names;
    string log_name = "";
    //в какие файла и что сохранять
    //частоты записи
};

class InputProcessor{
public:
    AorticInput a_in;
    string aorticInputFile = "";

    vector<LeafletInput> l_ins;
    string leafletInputFile = "";

    SewEnergyParams sep;
    string sewEnergyParamsInputFile = "";

    CameraParams cp;
    string cameraParamsInputFile = "";

    std::vector<ObjectTraits> ots;
    string objectTraitsInputFile = "";

    SolverParams sp;
    string solverParamsInputFile = "";

    RequiredOutput ro;
    string outputTemplateInputFile = "";

    InputProcessor() {}
    InputProcessor(int argc, char* argv[]) { InputProcessorInit(argc, argv); } ;
    void InputProcessorInit(int argc, char* argv[]);

    friend string to_string(const InputProcessor& ip) { return ip.__input_text; }

private:
    void processConfigFile(const char* filename);
    void processConfigFile_aorta(string filename);
    void processConfigFile_leaflets(string filename);
    void processConfigFile_sew_energy(string filename);
    void processConfigFile_solver(string filename);
    void processConfigFile_out_template(string filename);

    string __input_text;
};


#endif