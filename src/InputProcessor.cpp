#include "InputProcessor.h"
#include "nlohmann/json.hpp"
#include <stdexcept>
#include <iostream>
#include <map>
#include <algorithm>
#include <cctype>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <experimental/filesystem>


using json = nlohmann::json;

json get_json(string filename){
    std::ifstream js;
    js.open(filename);
    if (! js.good()) {
        throw std::runtime_error("Unsuccessful attempt to open file \"" + filename + "\"");
    }
    json conf;
    js >> conf;
    js.close();

    return std::move(conf);
}

void no_field(string field){
    std::cout << "Warning: Can't find the field \"" + field + "\". For this will be used defaults.\n";
}

array<double, 3> from_vec(vector<double> v){
    array<double, 3> res;
    for (int i = 0; i < 3 && i < v.size(); ++i)
        res[i] = v[i];
    return std::move(res);
}

#define SET(X, Y, Z) if (Y.count(Z)) X = Y[Z]

#define SETW(X, Y, Z, W) if (Y.count(Z)) X = Y[Z]; else no_field(W)

void InputProcessor::processConfigFile_aorta(string filename) {
    if (aorticInputFile.size() && filename != aorticInputFile)
        cout << "Warning: The field \"aorta\" in file " + aorticInputFile + " will be ignored.\n";
    aorticInputFile = filename;
    cameraParamsInputFile = filename;
    json conf = get_json(filename);
    __input_text += "aorta: " + to_string(conf) + "\n";
    if (!conf.count("aorta")){
        no_field("aorta");
        return;
    }
    if (conf["aorta"].count("mesh_file"))
        a_in.aorta_file = conf["aorta"]["mesh_file"];
    else
        no_field("aorta:mesh_file");

    if (conf["aorta"].count("boundary_file"))
        a_in.leaflet_boundaries = conf["aorta"]["boundary_file"];
    else
        no_field("aorta:boundary_file");

    if (conf["aorta"].count("blood_flow"))
        a_in.main_flow_direction = from_vec(conf["aorta"]["blood_flow"]);

    if (conf["aorta"].count("camera_params")){
        json cam = conf["aorta"]["camera_params"];
        if (cam["use_camera"]){
            if (cam.count("origin"))
                cp.origin = from_vec(cam["origin"]);
            SET(cp.pitch, cam, "pitch");
            SET(cp.yaw, cam, "yaw");
        }
    }
}

int getModel(string ModelName){
    transform(ModelName.begin(), ModelName.end(), ModelName.begin(), ::tolower);
    map<string, ElasticModels> models = {
            {"msm", EMOD_MSM},
            {"mass-spring", EMOD_MSM},
            {"mass-spring model", EMOD_MSM},
            {"trqs", EMOD_TRQS},
            {"st venan", EMOD_TRQS},
            {"neogook", EMOD_NEOGOOK},
            {"neogokean", EMOD_NEOGOOK},
            {"neo-gookean", EMOD_NEOGOOK},
            {"reinforcing", EMOD_REINFORCING},
            {"compressible neogook", EMOD_REINFORCING}
    };
    if (!models.count(ModelName)) return -1;
    else return models[ModelName];
}

void InputProcessor::processConfigFile_leaflets(string filename) {
    if (leafletInputFile.size() && filename != leafletInputFile)
        cout << "Warning: The field \"leaflets\" in file " + leafletInputFile + " will be ignored.\n";
    leafletInputFile = filename;
    objectTraitsInputFile = filename;
    json conf = get_json(filename);
    __input_text += "leaflets: " + to_string(conf) + "\n";
    if (!conf.count("leaflets")){
        no_field("leaflets");
        return;
    }
    conf.size();
    int n = conf["leaflets"].size();
    l_ins.resize(n);
    ots.resize(n);
    for (int i = 0; i < n; ++i){
        json leaf = conf["leaflets"][i];
        if (!leaf.count("mesh")){
            no_field("leaflets[" + to_string(i) + "]:mesh");
        } else {
            json mesh = leaf["mesh"];
            if (!mesh.count("generative")) {
                no_field("leaflets[" + to_string(i) + "]:mesh:generative");
            }
            else
                l_ins[i].generative = mesh["generative"];
            SET(l_ins[i].model, mesh, "model");
            if (mesh.count("model_params"))
                l_ins[i].model_params = mesh["model_params"].get<vector<double >>();
            if (mesh.count("mesh_file"))
                l_ins[i].model = mesh["mesh_file"];
            else
                if (!l_ins[i].generative)
                    no_field("leaflets[" + to_string(i) + "]:mesh:mesh_file");
            if (mesh.count("main_axis_direction"))
                l_ins[i].main_axis_direction = from_vec(mesh["main_axis_direction"]);
            if (mesh.count("second_line_symmetry"))
                l_ins[i].second_line_symmetry = from_vec(mesh["second_line_symmetry"]);
        }
        if (!leaf.count("traits")){
            no_field("leaflets[" + to_string(i) + "]:traits");
        } else {
            json traits = leaf["traits"];
            bool has_model = true;
            if (!traits.count("elastic_model")) {
                no_field("leaflets[" + to_string(i) + "]:traits:elastic_model");
                has_model = false;
            } else ots[i].elastic_model_type = getModel(traits["elastic_model"]);
            if (traits.count("elastic_params")){
                if (!has_model)
                    cout << "Warning: in leaflets[\" + to_string(i) + \"]:traits: \"elastic_model\" not setted but \"elastic_params\" setted.\n";
                ots[i].elastic_model_params = traits["elastic_params"].get<vector<double>>();
            }
            SET(ots[i].pressure, traits, "pressure");
            SET(ots[i].collision_margin, traits, "collision_margin");
            SET(ots[i].renderLabel, traits, "render_label");
        }
    }
}
void InputProcessor::processConfigFile_sew_energy(string filename) {
    if (sewEnergyParamsInputFile.size() && filename != sewEnergyParamsInputFile)
        cout << "Warning: The field \"sew_energy_params\" in file " + sewEnergyParamsInputFile + " will be ignored.\n";
    sewEnergyParamsInputFile = filename;
    json conf = get_json(filename);
    __input_text += "sew_energy: " + to_string(conf) + "\n";
    if (!conf.count("sew_energy_params")){
        no_field("sew_energy_params");
        return;
    }
    conf = conf["sew_energy_params"];
    SET(sep.plane_w, conf, "plane_constraint_weight");
    SET(sep.sqr_length_w, conf, "mesh_edge_constraint_weight");
    SET(sep.digedral_angle_w, conf, "mesh_digedral_element_weight");
    SET(sep.convexity_w, conf, "mesh_convexity_factor");
    if (conf.count("energy_minimizer_params")){
        json emp = conf["energy_minimizer_params"];
        SET(sep.sp.freq, emp, "freq");
        SET(sep.sp.step_sz, emp, "step_size");
        SET(sep.sp.tol, emp, "tol");
        SET(sep.sp.epsabs, emp, "epsabs");
        SET(sep.sp.maxits, emp, "maxits");
        SET(sep.sp.time, emp, "max_time_per_leaflet");
    }
}
void InputProcessor::processConfigFile_solver(string filename) {
    if (solverParamsInputFile.size() && filename != solverParamsInputFile)
        cout << "Warning: The field \"solver_params\" in file " + solverParamsInputFile + " will be ignored.\n";
    solverParamsInputFile = filename;
    json conf = get_json(filename);
    __input_text += "solver: " + to_string(conf) + "\n";
    if (!conf.count("solver_params")){
        no_field("solver_params");
        return;
    }
    conf = conf["solver_params"];
    SETW(sp.delta, conf, "delta", "solver_params:delta");
    SETW(sp.eps, conf, "eps", "solver_params:eps");
    SET(sp.max_possible_shift_scale, conf, "max_shift_scale");
    SET(sp.max_possible_recommend_shift_scale, conf, "max_recommend_shift_scale");
}

string getPostfix(){
    return to_string(time(NULL));
}

void InputProcessor::processConfigFile_out_template(string filename){
    if (outputTemplateInputFile.size() && filename != outputTemplateInputFile)
        cout << "Warning: The field \"result_template\" in file " + outputTemplateInputFile + " will be ignored.\n";
    outputTemplateInputFile = filename;
    json conf = get_json(filename);
    __input_text += "out_template: " + to_string(conf) + "\n";
    if (!conf.count("result_template")){
        no_field("result_template");
        return;
    }
    conf = conf["result_template"];
    SET(ro.use, conf, "save_results");
    if (!ro.use) return;
    ro.postfix = getPostfix();
    __input_text += "postfix: " + ro.postfix + "\n";
    bool use_dir = false;
    bool create_dir = false;
    bool use_dir_postfix = false;
    SET(use_dir, conf, "use_common_result_directory");
    if (use_dir) {
        SET(create_dir, conf, "create_directory");
        SET(ro.res_dir, conf, "common_result_directory");
        if (create_dir) {
            SET(use_dir_postfix, conf, "use_result_directory_postfix");
            if (use_dir_postfix) {
                while (ro.res_dir.length() && ro.res_dir[ro.res_dir.length() - 1] == '/')
                    ro.res_dir.pop_back();
                ro.res_dir += ro.postfix;
            }
            if (ro.res_dir.length() && ro.res_dir[ro.res_dir.length() - 1] != '/') ro.res_dir += "/";
            if (!experimental::filesystem::create_directory(ro.res_dir)) {
                std::cout << "Failed to create directory " << ro.res_dir << "\n";
                exit(-1);
            }
        }
    }
    bool save_aorta = false;
    SET(save_aorta, conf, "save_aorta_to_nts");
    if (save_aorta) SET(ro.aorta_to_nts, conf, "aorta_to_nts_name");
    bool save_leafs = false;
    SET(save_aorta, conf, "save_leaflets");
    if (save_aorta) {
        SET(ro.divide_leaflets, conf, "divide_leaflets_per_files");
        if (conf.count("leaflets_names"))
            ro.leaf_names = conf["leaflets_names"].get<vector<string>>();
    }
    SET(ro.log_name, conf, "logfile_name");
    bool use_res_postfix = false;
    SET(use_res_postfix, conf, "use_results_postfix");
    struct AddPostfix{
        static void addPosix(string& name, const string& posix){
            if (name.length())
                name.insert(name.find_last_of("."), posix);
        }
    };
    AddPostfix::addPosix(ro.aorta_to_nts, ro.postfix);
    ro.aorta_to_nts = ro.res_dir + ro.aorta_to_nts;
    AddPostfix::addPosix(ro.log_name, ro.postfix);
    ro.log_name = ro.res_dir + ro.log_name;
    for (int i = 0; i < ro.leaf_names.size(); ++i) {
        AddPostfix::addPosix(ro.leaf_names[i], ro.postfix);
        ro.leaf_names[i] = ro.res_dir + ro.leaf_names[i];
    }
}

void InputProcessor::processConfigFile(const char* filename){
    json conf = get_json(filename);
    __input_text += "main: " + to_string(conf) + "\n";
//    string mesh = conf["aorta"]["mesh_file"];
//    std::vector<double> test = conf["aorta"]["blood_flow"];
//    std::cout << test.size() << test[0] << test[1] << test[2];
    using Processor = void;
    std::vector<std::pair<string, string>> conf_types = {
            {"aorta_file", "aorta"},
            {"leaflets_file", "leaflets"},
            {"sew_energy_file", "sew_energy_params"},
            {"solver_file", "solver_params"}
    };
    string thisfile = filename;
    int __pos = thisfile.find_last_of('/');
    __pos = (__pos == thisfile.size() - 1) ? 0 : __pos;
    string dir = thisfile.substr(0, __pos + ((thisfile.size() > 1) ? 1 : 0));
//    cout << "array size = " << conf["leaflets"].size() << std::endl;
//    exit(-1);
    if (conf.count("aorta_file")) {
        processConfigFile_aorta(dir + string(conf["aorta_file"]));
        if (conf.count("aorta"))
            cout << "Warning: The field \"aorta\" in file " + thisfile + " will be ignored.\n";
    }
    else processConfigFile_aorta(thisfile);
    if (conf.count("leaflets_file")) {
        processConfigFile_leaflets(dir + string(conf["leaflets_file"]));
        if (conf.count("leaflets"))
            cout << "Warning: The field \"leaflets\" in file " + thisfile + " will be ignored.\n";
    }
    else processConfigFile_leaflets(thisfile);
    if (conf.count("sew_energy_file")) {
        processConfigFile_sew_energy(dir + string(conf["sew_energy_file"]));
        if (conf.count("sew_energy_params"))
            cout << "Warning: The field \"sew_energy_params\" in file " + thisfile + " will be ignored.\n";
    }
    else processConfigFile_sew_energy(thisfile);
    if (conf.count("solver_file")) {
        processConfigFile_solver(dir + string(conf["solver_file"]));
        if (conf.count("solver_params"))
            cout << "Warning: The field \"solver_params\" in file " + thisfile + " will be ignored.\n";
    }
    else processConfigFile_solver(thisfile);
    if (conf.count("result_template_file")) {
        processConfigFile_out_template(dir + string(conf["result_template_file"]));
        if (conf.count("result_template"))
            cout << "Warning: The field \"result_template\" in file " + thisfile + " will be ignored.\n";
    }
    else processConfigFile_out_template(thisfile);
}

void InputProcessor::InputProcessorInit(int argc, char **argv) {
    int i = 0;
    if (argc == 1) goto helpMessage;
    for (i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            helpMessage:
            std::cout << "Help message: " << std::endl;
            std::cout << "Command line options: " << std::endl;
            std::cout << "-c, --config <Main configuration parameters file name>" << std::endl;
            std::cout << "-ac, --aortacnfg <Aortic configuration parameters file name>" << std::endl;
            std::cout << "-lc, --leafscnfg <Leaflets configuration parameters file name>" << std::endl;
            std::cout << "-sc, --solvercnfg <Solver configuration parameters file name>" << std::endl;
            exit(0);
        }
        //Parameters file name found with -d or --database options
        if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--config") == 0) {
            std::cout << "Main configuration file found: " << argv[i + 1] << std::endl;
            processConfigFile(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-ac") == 0 || strcmp(argv[i], "--aortacnfg") == 0) {
            std::cout << "Aortic configuration file externally set to : " << argv[i + 1] << std::endl;
            processConfigFile_aorta(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-lc") == 0 || strcmp(argv[i], "--leafscnfg") == 0) {
            std::cout << "Leaflets configuration file externally set to : " << argv[i + 1] << std::endl;
            processConfigFile_leaflets(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-sc") == 0 || strcmp(argv[i], "--solvercnfg") == 0) {
            std::cout << "Solver configuration file externally set to : " << argv[i + 1] << std::endl;
            processConfigFile_solver(argv[i + 1]);
            i++;
            continue;
        }
    }
}

//TODO: абсолютные и относительные пути до папки результатов

#undef SET
#undef SETW