#include "Simulator.h"
#include "util_sim.h"


/**Function*************************************************************

  Synopsis    [Initailize simulator]

  Description [This function will set #qubits n, construct initial state, and enable dynamic reordering]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::init_simulator(int nQubits)
{
    n = nQubits; // set the number n here
    manager = Cudd_Init(n, n, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    int *constants = new int[n];
    for (int i = 0; i < n; i++)
        constants[i] = 0; // TODO: costom initial state
    measured_qubits_to_clbits = std::vector<std::vector<int>>(n, std::vector<int>(0));
    init_state(constants);
    delete[] constants;
    if (isReorder) Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);
}

bool Simulator::check_equ()
{
    if (Single_Bdd == Initial_BDD){
        std::cout<<"Test "<<current_bit<<" succeeded."<<std::endl;
        return true;
    }
    else{
        std::cout<<"Test "<<current_bit<<" failed."<<std::endl;
        return false;
    }
}
/**Function*************************************************************

  Synopsis    [parse and simulate the qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void Simulator::sim_qasm_eqc(int cur_bit)
{
    // std::cout<<n<<std::endl;
    // std::cout<<cur_bit<<std::endl;
    std::string inStr;
    std::stringstream inFile_ss(qasmfile);
    while (getline(inFile_ss, inStr))
    {
        inStr = inStr.substr(0, inStr.find("//"));
        if (inStr.find_first_not_of("\t\n ") != std::string::npos)
        {
            std::stringstream inStr_ss(inStr);
            getline(inStr_ss, inStr, ' ');
            if (inStr == "qreg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                // init_simulator(stoi(inStr));
            }
            else if (inStr == "creg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                nClbits = stoi(inStr);                
            }
            else if (inStr == "OPENQASM"){;}
            else if (inStr == "include"){;}
            else
            {
                if (inStr == "x")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    PauliX(stoi(inStr), cur_bit);
                }
                else if (inStr == "cx")
                {
                    std::vector<int> cont(1);
                    std::vector<int> ncont(0);
                    int targ;
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    cont[0] = stoi(inStr);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    targ = stoi(inStr);
                    Toffoli(targ, cont, ncont, cur_bit);
                    cont.clear();
                    ncont.clear();
                }
                else if (inStr == "swap")
                {
                    int swapA, swapB;
                    std::vector<int> cont(0);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        if (i == 0)
                            swapA = stoi(inStr);
                        else
                            swapB = stoi(inStr);
                    }
                    Fredkin(swapA, swapB, cont, cur_bit);
                    cont.clear();
                }
                else if (inStr == "cswap")
                {
                    int swapA, swapB;
                    std::vector<int> cont(1);
                    for (int i = 0; i < 3; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        if (i == 0)
                            cont[i] = stoi(inStr);
                        else if (i == 1)
                            swapA = stoi(inStr);
                        else
                            swapB = stoi(inStr);
                    }
                    Fredkin(swapA, swapB, cont, cur_bit);
                    cont.clear();
                }
                // else if (inStr == "mcswap")
                // {
                //     std::vector<int> cont(0);
                //     std::vector<int> ncont(0);
                //     int targ;
                //     getline(inStr_ss, inStr, '[');
                //     while(getline(inStr_ss, inStr, ']'))
                //     {
                //         cont.push_back(stoi(inStr));
                //         getline(inStr_ss, inStr, '[');
                //     }
                //     targ = cont.back();
                //     cont.pop_back();
                // }
                else if (inStr == "ccx" || inStr == "mcx")
                {
                    std::vector<int> cont(0);
                    std::vector<int> ncont(0);
                    int targ;
                    getline(inStr_ss, inStr, '[');
                    while(getline(inStr_ss, inStr, ']'))
                    {
                        cont.push_back(stoi(inStr));
                        getline(inStr_ss, inStr, '[');
                    }
                    targ = cont.back();
                    cont.pop_back();
                    Toffoli(targ, cont, ncont, cur_bit);
                    cont.clear();
                    ncont.clear();
                }
                else
                {
                    std::cerr << std::endl
                            // << "[warning]: Gate \'" << inStr << "\' is not supported in this simulator. The gate is ignored ..." << std::endl;
                            << "[warning]: Syntax \'" << inStr << "\' is not supported in this simulator. The line is ignored ..." << std::endl;
                }
            }
        }
    }
    // std::cout<<cur_bit<<std::endl;
    // std::cout<<i<<std::endl;
    // if (!check_equ()){
    //     flag = false;
    //     std::cout<<"The circuits are not equivalent."<<std::endl;
    // }
}
/**Function*************************************************************

  Synopsis    [simulate the circuit described by a qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm(std::string qasm)
{
    std::string inStr0;
    std::stringstream inFile_ss0(qasm);
    bool flag = true;
    while (getline(inFile_ss0, inStr0))
    {
        inStr0 = inStr0.substr(0, inStr0.find("//"));
        if (inStr0.find_first_not_of("\t\n ") != std::string::npos)
        {
            std::stringstream inStr_ss0(inStr0);
            getline(inStr_ss0, inStr0, ' ');
            if (inStr0 == "qreg")
            {
                getline(inStr_ss0, inStr0, '[');
                getline(inStr_ss0, inStr0, ']');
                init_simulator(stoi(inStr0));
                break;
            }
            else {;}
        }
    }
    qasmfile = qasm;
    bool equivalent = true;
    for(int i=0;i<n;i++){
        sim_qasm_eqc(i);
        if(Initial_BDD[i] == Single_Bdd[i]){
            std::cout<<"Test "<<i<<" succeeded."<<std::endl;
        }
        else{
            std::cout<<"Test "<<i<<" failed."<<std::endl;
            equivalent = false;
            break;
        }
        Cudd_RecursiveDeref(manager, Single_Bdd[i]);
    }
    if(equivalent) std::cout<<"Equivalent"<<std::endl;
    else std::cout<<"not Equivalent"<<std::endl;
}



/**Function*************************************************************

  Synopsis    [print state vector and distribution of sampled outcomes]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::print_results()
{
    // write output string based on state_count and statevector
    std::unordered_map<std::string, int>::iterator it;
    
    run_output = "{";
    if (state_count.begin() != state_count.end()){
        run_output += "\"counts\": { ";
        for (it = state_count.begin(); it != state_count.end(); it++)
        {
            if (std::next(it) == state_count.end())
                run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second);
            else
                run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second) + ", ";
        }
        run_output += " }";
        run_output += (statevector != "null") ? ", " : ""; 
    }    

    run_output += (statevector != "null") ? "\"statevector\": " + statevector + " }" : " }";
    //return;
    std::cout << run_output << std::endl;
}