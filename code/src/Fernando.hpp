#include "Instance.hpp"
#include "lp.hpp"
#include <unordered_set>
#include <deque>
extern "C"{
#include "cgraph/cgraph.h"
#include "cgraph/clique_separation.h"
}

#include "Hash.cpp"


class Fernando{
public:
    Fernando( Instance &_inst );


    virtual ~Fernando();
private:
    Instance &inst_;

    std::vector< int > fim;
    std::vector< std::string > names;

    std::vector< std::vector< std::vector< int > > > xIdx_;
    std::vector< std::vector< std::vector< int > > > eIdx_;
    std::vector< std::vector< int > > fIdx_;

    std::vector<std::vector<std::vector<std::vector<int>>>> process;
    std::vector<std::vector<std::vector<std::vector<int>>>> enter_flow;

    double lifting(double c);
    void lifting_linear();
    void lifting_binario();
    int manual_cuts();
    int qtd_manual_cuts = 0;
    //std::vector<std::vector<int>> variables_pack;
    std::unordered_set<std::vector<int>> variables_pack;

    int totalCliques = 0;
    int totalFenchel = 0;
    int qtdLimitesEnumeracao = 0;
    int qtdLimiteVarsEnumeracao = 0;
    int qtdJanelaModificada = 0;
    int cIdx_;
    LinearProgram *mip;
    double teto(double v);
    CGraph *cgraph;
    bool clique = true;
    bool continuo = true;
    bool binario = false;
    void optimize();

    int cliques(int *idxs,double *coefs);

    void cgraph_creation();

    int oddHoles();

    void buildProblem();
    void buildCliqueCuts();

    /* melhorar o lst */
    void combinacao(int job, unsigned int tam, std::vector<int> &vec, std::vector<std::vector<int> > &combinacoes);
    void reduz_lst_kondili(int k_max);

    /*corte de fenchel*/
    bool insertVar(std::vector<S> sol, S var);
    bool backtrack(int j, int op, int ti, std::vector<S> sol);
    template <typename T> bool isSubset(std::vector<T> &A, std::vector<T> &B);
    bool limite_enumeracao;
    void enumeracao_fenchel(const std::deque<S> &vars, int index, std::vector<std::vector<S>> &solutions, std::vector<S> solution);
    template <typename T> bool dominancia(std::vector<T> &vec, std::vector<std::vector<T>> &set);
    int fenchel_tempo(int ti, int tf); // pega as variáveis por um intervalo de tempo para fazer o corte
    int fenchel_vars(std::deque<S> &vars, std::deque<S> &vars1); // pega as variáveis pelo intervalo de variáveis maiores que zero
    int executeFenchel();
    std::vector< std::vector< int > > enum_time; // armazena o tempo das variáveis a serem alocadas
    std::vector< std::vector< int > > novoEst; // armazena o novo est do backtracking
    std::vector< std::vector< int > > novoLst; // armazena o novo lst do backtracking
    std::vector< std::vector< int > > modificadoresEst; // armazena o modificador aplicado no novo est do backtracking
    std::vector< std::vector< int > > modificadoresLst; // armazena o modificador aplicado no novo lst do backtracking

    int runFenchel(std::deque<S> &vars, std::deque<S> &vars1);

    void setAuxiliaresBacktrack(const S &var);
    void desetAuxiliaresBacktrack(const S &var);

    bool canInsert(S var);
};
